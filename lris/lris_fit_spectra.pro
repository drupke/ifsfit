;+
; NAME:
;     LRIS_FIT_SPECTRA
;
; PURPOSE:
;     Wrapper script for fitting LRIS spectra.
;
; EXPLANATION:
;
; CALLING SEQUENCE
;     lris_fit_spectra,mask,slitrange=,aprange=,/time,/verbose,/onecomp
;
; INPUTS:
;     mask - String name of mask, as specified in LRIS_INITFIT_SPECTRA
;     slitrange - scalar or two-element array specifying range of slits to fit
;     aprange - scalar or two-element array specifying range of apertures to fit
;     time - selects timing of continuum and line fits for each spectrum
;     verbose - selects verbose output from fitting procedures
;     onecomp - fit only a single component to each spectrum
;
; OUTPUT:
;     A binary (.genx) file containing the results of the fit.
;
; METHOD:
;
; REVISION HISTORY:
;     09may13  David Rupke  created
;     09jun07  David Rupke  added multiple components
;     10jan28  DSNR         re-written to fit in observed frame
;
;-

pro lris_fit_spectra,lab,mask,slitrange=slitrange,aprange=aprange,$
                     time=time,verbose=verbose,onecomp=onecomp

starttime = systime(1)

; Set fit range in rest frame
fitran_rest = [3700d,6800d]
; Intercept of blue tilt, in rest frame
blueint = 4500d
; Continuum fitting function
fcncontfit='lris_fit_continuum'
fcnz = 'lris_redshift_spec'

if ~keyword_set(time) then time=0
if keyword_set(verbose) then quiet=0 else quiet=1

; Get fit initialization
initdat = lris_initfit_spectra(lab,mask)
infile = initdat.infile
errfile = initdat.errfile
outfile = initdat.outfile
compfile = initdat.compfile
initfile = initdat.initfile
zguess = initdat.zinit
zfix = initdat.zfix
defaultvdisp = initdat.vdisp
stitchwave = initdat.stitchwave
startempfile = initdat.startempfile
strong = initdat.stronglines
addexp = initdat.addexp
fcninitpar = initdat.initparinfo

; Set up multiple components
if ~ keyword_set(onecomp) AND file_test(compfile) then begin
   readcol300,compfile,c_slit,c_ap,c_icent,format='(I,I,D)',/silent,skip=1
   mc_refwave = 6562.8d
endif
if file_test(initfile) then $
   readcol300,initfile,in_slit,in_ap,in_zinit,in_dzstel,$
              in_vdisp,in_masksig,in_strong,in_noem,$
              format='(I,I,D,D,D,D,I,I)',/silent,skip=1

; Set up emission-line masking parameters
defaultmasksig = 5d

; Get number of slits
fits_info,infile+'.fits',/silent,N_ext=next
next+=1
print,'Found ',next,' slits in ',infile,format='(/,A,I0,A,A,/)'
if n_elements(zfix) eq 1 then zfix = intarr(next)+zfix
if n_elements(zguess) eq 1 then zguess = intarr(next)+zguess
if n_elements(strong) eq 1 then strong = intarr(next)+strong

; Cycle through slits
if ~ keyword_set(slitrange) then slitrange=[1,next] $
else if n_elements(slitrange) eq 1 then slitrange = [slitrange,slitrange]
slitrange = fix(slitrange)
for i=slitrange[0]-1,slitrange[1]-1 do begin

   print,'Processing slit ',i+1,' of ',next,'...',format='(A,I0,A,I0,A)'

   spec = readspec(infile+'.fits',extension=i)
   specerr = readspec(errfile+'.fits',extension=i)
   wave = spec[*,0]
   naps = n_elements(spec[0,*])-1

;  Cycle through apertures
   if ~ keyword_set(aprange) then loc_aprange=[1,naps] $
   else if n_elements(aprange) eq 1 then loc_aprange=[aprange,aprange] $
   else loc_aprange = aprange
   loc_aprange = fix(loc_aprange)
   for j=loc_aprange[0],loc_aprange[1] do begin

;     Common initializations
      flux = spec[*,j]
      err = specerr[*,j]

      print,'  Processing aperture ',j,' of ',naps,'...',$
            format='(A,I0,A,I0,A)'
      slaplab='sl'+string(i+1,format='(I02)')+'_ap'+$
              string(j,format='(I02)')
      outfileslap = outfile + '_' + slaplab

;     Total number of velocity components.  Default is 1 component.
      ncomp=1
      if file_test(initfile) then begin
         in_ind = where(in_slit eq i+1 AND in_ap eq j,count)
         if count eq 1 then $
            if in_noem[in_ind] eq 1 then ncomp=0
      endif
;     Initial LRIS set of redshifts (blue, red, blue tilt, stellar),
;     and estimated shifts in redshift space for extra components.
      dzstel = 0d
      vdisp = defaultvdisp
      masksig = defaultmasksig
      tmp = zguess[i]
      zinit = double(tmp[0])
      dostrong = strong[i]
      disperse = 0
      if ncomp gt 0 then begin
         zinit = refinez(zinit,wave/(1+zinit),flux)
         if ~ quiet then print,'Redshift adjusted from ',zguess[i],$
                               ' to ',zinit,$
                               format='(A0,D0.5,A0,D0.5)'
      endif
      if file_test(initfile) then begin
         in_ind = where(in_slit eq i+1 AND in_ap eq j,count)
         if count eq 1 then begin
            if in_zinit[in_ind] ne -1 then begin
               tmp = in_zinit[in_ind]
               zinit = double(tmp[0])
            endif
            if in_dzstel[in_ind] ne -1 then begin
               tmp = in_dzstel[in_ind]
               dzstel = double(tmp[0])
            endif
            if in_vdisp[in_ind] ne -1 then begin
               disperse = 1
               tmp = in_vdisp[in_ind]
               vdisp = double(tmp[0])
            endif
            if in_masksig[in_ind] ne -1 then begin
               tmp = in_masksig[in_ind]
               masksig = double(tmp[0])
            endif
            if in_strong[in_ind] ne -1 then begin
               tmp = in_strong[in_ind]
               dostrong = double(tmp[0])
            endif
         endif
      endif
      zred = zinit
      zblue = zinit
;   Case of multiple components
      if ~ keyword_set(onecomp) AND file_test(compfile) then begin
         mc_ind = where(c_slit eq i+1 AND c_ap eq j,count)
         ncomp+=count
         if count ge 1 then begin
            dzcomp = c_icent[mc_ind]/mc_refwave
            zred = [zinit,zinit + dzcomp]
            zblue = [zinit,zinit + dzcomp]
         endif
      endif
;     Initialize redshift structure
      z = {star:zinit+dzstel,$
           gas: zred,$
           blue: zblue,$
           bluetilt:0d,$
           blueint:blueint}


;     Get linelist
      linelist = lris_initlinelist(strong=dostrong)
      
;     Fit the First

      argsinitpar={stitchwave:stitchwave,zfix:zfix[i]}
      argscntpar={stitchwave:stitchwave}
      argspoly={addexp:addexp}

;     Remove telluric line
      gd_indx = where(wave lt 6855 OR wave gt 6945)
      wave = wave[gd_indx]
      flux = flux[gd_indx]
      err = err[gd_indx]

;     Redshift fit range
      fitran = lris_redshift_spec(fitran_rest,z)

      structinit = fit_spectrum(wave,flux,err,startempfile,z,ncomp,$
                                /subtract,obj_id=slaplab,$
                                linelist=linelist,quiet=quiet,$
                                time=time,fitran=fitran,$
                                fcninitpar=fcninitpar,$
                                argsinitpar=argsinitpar,$
                                disperse=disperse,fcncontfit=fcncontfit,$
                                argscontfit=argscntpar,vdisp=vdisp,$
                                fcnz=fcnz,argspoly=argspoly)

      if ~ quiet then print,'FIT STATUS: ',structinit.fitstatus

      if ncomp gt 0 then begin

;     Repeat fit with new z estimate and better estimate for line widths

;     Re-order components in struct.param based on Halpha flux if
;     there is more than 1 component.
         if ncomp gt 1 then lris_orderlines,structinit
;     Update z array to reflect new emission-line redshifts
         lris_updatez,structinit.param,z
;     New stellar redshift estimate
         z.star = z.gas[0] + dzstel
;     Make emission-line masks based on best-fit parameters
         linepars = sepfitpars(structinit.param,structinit.perror)
         masklines = reform(linepars.wave[*,0:ncomp-1],$
                            n_elements(linepars.wave[*,0:ncomp-1]))
         maskwidths = masksig*reform(linepars.sigma[*,0:ncomp-1],$
                                     n_elements(linepars.sigma[*,0:ncomp-1]))
;     Redshift fit range
         fitran = lris_redshift_spec(fitran_rest,z,stitchwave=stitchwave)

;     This time, when we redshift everything include the blue/red
;     differences.
         argsz={stitchwave:stitchwave}


;     Fit the Second
         struct = fit_spectrum(wave,flux,err,startempfile,z,ncomp,$
                               /subtract,obj_id=slaplab,$
                               linelist=linelist,quiet=quiet,$
                               time=time,fitran=fitran,$
                               fcninitpar=fcninitpar,$
                               argsinitpar=argsinitpar,$
                               maskwidths=maskwidths,masklines=masklines,$
                               sigguess=linepars.sigma,$
                               peakguess=linepars.fluxpk,$
                               disperse=disperse,fcncontfit=fcncontfit,$
                               argscontfit=argscntpar,vdisp=vdisp,$
                               fcnz=fcnz,argsz=argsz,argspoly=argspoly)

         if ~ quiet then print,'FIT STATUS: ',struct.fitstatus
      
;     Update z array to reflect new emission-line redshifts.  Do *not*
;     update stellar redshift in z array, since that hasn't changed.
         lris_updatez,struct.param,z

      endif else begin
         
         struct = structinit
      
      endelse
         
;     Update structure with z array and fit range.
      struct = jjadd_tag(struct,'z',z)
      struct = jjadd_tag(struct,'fitrange',fitran,/array_tag)
      
;     Save result to binary file
;      savegen,struct=struct,file=outfileslap+'.genx'
      save,struct,file=outfileslap+'.xdr'

   endfor
   
endfor

print,'Total time for calculation: ',systime(1)-starttime,' s.',$
      format='(/,A0,I0,A0,/)'

end
