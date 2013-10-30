;
; NAME:
;     LRIS_FEB09_FIT_SPECTRA
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
;     09jun07  DSNR         added multiple components
;     10apr21  DSNR         converted to observed-frame fits
;-

pro lris_feb09_fit_spectra,lab,mask,slitrange=slitrange,aprange=aprange,$
                           time=time,verbose=verbose

starttime = systime(1)

; Set fit range in rest frame
fitran_rest = [3700d,6800d]
; Intercept of blue tilt, in rest frame
blueint = 4500d

if ~keyword_set(time) then time=0
if keyword_set(verbose) then quiet=0 else quiet=1

; Get fit initialization
initdat = lris_feb09_initfit_spectra(lab,mask)
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
fcninitpar = initdat.initparinfo
fcnlinefit = initdat.fcnlinefit
fcncontfit = initdat.fcncontfit
redord = initdat.redord
fcnz = 'lris_redshift_spec'

; Set up multiple components
if file_test(compfile) then begin
   readcol300,compfile,c_slit,c_ap,c_icent,format='(I,I,D)',/silent,skip=1
   mc_refwave = 6562.8d
endif
if file_test(initfile) then $
   readcol300,initfile,in_slit,in_ap,in_zinit,in_dzstel,in_vdisp,in_masksig,in_floatf,in_floatl,in_floats,in_noem,format='(I,I,D,D,D,D,D,D,D,I)',/silent,skip=1

; Set up emission-line masking parameters
defaultmasksig = 5d

; Get number of slits
fits_info,infile+'.fits',/silent,N_ext=next
next+=1
print,'Found ',next,' slits in ',infile,format='(/,A,I0,A,A,/)'
if n_elements(zguess) eq 1 then zguess = intarr(next)+zguess
if n_elements(zfix) eq 1 then zfix = intarr(next)+zfix
if n_elements(strong) eq 1 then strong = intarr(next)+strong
if n_elements(fcncontfit) eq 1 then fcncontfit = $
   make_array(next,/string,value=fcncontfit)
if n_elements(redord) eq 1 then redord = intarr(next)+redord
;if n_elements(stitchwave) eq 1 then stitchwave = intarr(next)+stitchwave
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

      print,'  Processing aperture ',j,' of ',naps,'...',$
            format='(A,I0,A,I0,A)'
      slaplab='sl'+string(i+1,format='(I02)')+'_ap'+$
              string(j,format='(I02)')
      outfileslap = outfile + '_' + slaplab

;   Common initializations
      flux = spec[*,j]
      err = specerr[*,j]

;   Total number of velocity components.  Default is 1 component.
      ncomp=1
      if file_test(initfile) then begin
         in_ind = where(in_slit eq i+1 AND in_ap eq j,count)
         if count eq 1 then $
            if in_noem[in_ind] eq 1 then ncomp=0
      endif

;   Initial LRIS set of redshifts (blue, red, blue tilt, stellar), and
;   estimated shifts in redshift space for extra components.
      dzstel = 0d
      vdisp = defaultvdisp
      masksig = defaultmasksig
      tmp = zguess[i]
      zinit = double(tmp[0])
      disperse = 0
      floatcomp = 0
      dostrong = strong[i]
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
           if in_floatf[in_ind] ne -1 AND $
              in_floatl[in_ind] ne -1 AND $
              in_floats[in_ind] ne -1 then begin
              floatcomp = dblarr(3)
              floatcomp[0] = in_floatf[in_ind]
              floatcomp[1] = in_floatl[in_ind]
              floatcomp[2] = in_floats[in_ind]
           endif
         endif
      endif
      zred = zinit
      zblue = zinit
;   Case of multiple components
      if file_test(compfile) then begin
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

;     Initialize "broad" and "red" components for parameterizing Feb
;     09 red-side PSF.
;
;     Flux relative to narrow comp and sigma in A
      initbcomp = [0.03d,10d]
;     Wavelength offset, relative flux, sigma
      initrcomp = [3d,0.1d,1d]

;     Get linelist
      linelist = lris_initlinelist(strong=dostrong)
      
      
;     Fit the First


      argsinitpar={stitchwave:stitchwave,$
                   zfix:zfix[i],$
                   floatcomp:floatcomp,$
                   initbcomp:initbcomp,$
                   initrcomp:initrcomp}
      if fcncontfit[i] eq 'lris_fit_continuum' then $
         argscontfit={stitchwave:stitchwave,$
                      redord:redord[i]}
      if fcncontfit[i] eq 'fit_two_continua' then $
         argscontfit={stitchwave:stitchwave}
      argscntpar={stitchwave:stitchwave}
      argslinefit={floatcomp:floatcomp,$
                   stitchwave:stitchwave}

;     Remove telluric line
      gd_indx = where(wave lt 6855 OR wave gt 6945)
      wave = wave[gd_indx]
      flux = flux[gd_indx]
      err = err[gd_indx]

;     Set masking of floating component properly
      if ncomp gt 0 then begin
;     Initialize emission line list
         linewave = linelist.wave
         linelabel = linelist.label
         nlines = n_elements(linewave)
;     Redshift emission line list
         linewavez = dblarr(nlines,ncomp)
         for ic=0,ncomp-1 do $
            linewavez[*,ic] = $
            call_function(fcnz,linewave,z,$
                          _extra={icomp:ic,gas:1})
         masklines = reform(linewavez,ncomp*nlines)
         meansig = vdisp/299792d*6500d
         maskwidths = masksig*replicate(meansig,nlines*ncomp)
         if floatcomp[0] ne 0 then begin
            masklines = [masklines,replicate(floatcomp[1],ncomp)]
            maskwidths = [maskwidths,replicate(2d*floatcomp[2],ncomp)]
         endif
      endif else begin
         masklines = 0
         maskwidths = 0
      endelse

;     Redshift fit range
      fitran = lris_redshift_spec(fitran_rest,z)

      structinit = fit_spectrum(wave,flux,err,startempfile,z,ncomp,$
                                /subtract,obj_id=slaplab,$
                                linelist=linelist,quiet=quiet,$
                                time=time,fitran=fitran,$
                                fcninitpar=fcninitpar,$
                                argsinitpar=argsinitpar,$
                                disperse=disperse,$
                                fcncontfit=fcncontfit[i],$
                                argscontfit=argscntpar,vdisp=vdisp,$
                                fcnz=fcnz,argspoly=argspoly,$
                                fcnlinefit=fcnlinefit,$
                                argslinefit=argslinefit,$
                                masklines=masklines,maskwidths=maskwidths)

      if ~ quiet then print,'FIT STATUS: ',structinit.fitstatus

;   Repeat fit with new z estimate and better estimate for line widths

      if ncomp gt 0 then begin

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

;     Add broad Halpha component, and any floating component, to
;     masking algorithm
;      masklines = [masklines,replicate(6562.8d,ncomp)]
;      maskwidths = [maskwidths,replicate(100d,ncomp)]
         if floatcomp[0] ne 0 then begin
            floatcomp[0] = structinit.param[13]
            floatcomp[1] = structinit.param[14]
            floatcomp[2] = structinit.param[15]
            argsinitpar.floatcomp = floatcomp
            argslinefit.floatcomp = floatcomp
            masklines = [masklines,replicate(structinit.param[14],ncomp)]
            maskwidths = [maskwidths,replicate(2d*structinit.param[15],ncomp)]
         endif

;        This time, when we redshift everything include the blue/red
;        differences.
         argsz={stitchwave:stitchwave}
;        Redshift fit range
         fitran = lris_redshift_spec(fitran_rest,z,stitchwave=stitchwave)


         struct = fit_spectrum(wave,flux,err,startempfile,z,ncomp,$
                               linelist=linelist,quiet=quiet,$
                               obj_id=slaplab,$
                               time=time,fitran=fitran,/subtract,$
                               fcncontfit=fcncontfit[i],$
                               argscontfit=argscontfit,$
                               fcninitpar=fcninitpar,$
                               argsinitpar=argsinitpar,$
                               sigguess=linepars.sigma,$
                               peakguess=linepars.fluxpk,$
                               disperse=disperse,vdisp=vdisp,$
                               fcnz=fcnz,argsz=argsz,$
                               fcnlinefit=fcnlinefit,$
                               argslinefit=argslinefit,$
                               masklines=masklines,$
                               maskwidths=maskwidths)
         
         if ~ quiet then print,'FIT STATUS: ',struct.fitstatus

;        Update z array and fitted line wavelengths
         lris_updatez,struct.param,z

      endif else begin
         
         struct = structinit
      
      endelse
         
;     Update structure with z array and fit range.
      struct = jjadd_tag(struct,'z',z)
      struct = jjadd_tag(struct,'fitrange',fitran,/array_tag)

;     Save result to binary file
      save,struct,file=outfileslap+'.xdr'

   endfor

endfor

print,'Total time for calculation: ',systime(1)-starttime,' s.',$
      format='(/,A0,I0,A0,/)'

end
