;+
; NAME:
;     LRISBR_FIT_SPECTRA
;
; PURPOSE:
;     Wrapper script for fitting LRIS spectra.
;
; EXPLANATION:
;
; CALLING SEQUENCE
;     lrisbr_fit_spectra,mask,slitrange=,aprange=,/time,/verbose
;
; INPUTS:
;     mask - String name of mask, as specified in LRISBR_INITFIT_SPECTRA
;     slitrange - scalar or two-element array specifying range of slits to fit
;     aprange - scalar or two-element array specifying range of apertures to fit
;     time - selects timing of continuum and line fits for each spectrum
;     verbose - selects verbose output from fitting procedures
;
; OUTPUT:
;     A binary (.genx) file containing the results of the fit.
;
; METHOD:
;
; REVISION HISTORY:
;     10mar18  David Rupke  created
;
;-

pro lrisbr_fit_spectra,lab,mask,slitrange=slitrange,aprange=aprange,$
                       time=time,verbose=verbose

starttime = systime(1)

; Set fit range in rest frame
fitranBrest = [3700d,5500d]
fitranRrest = [5850d,6800d]
; Intercept of tilts, in rest frame
blueint = 4500d
redint = 6600d
; Continuum fitting function
fcnz = 'lrisbr_redshift_spec'

if ~keyword_set(time) then time=0
if  keyword_set(verbose) then quiet=0 else quiet=1

; Get fit initialization
initdat = lrisbr_initfit_spectra(lab,mask)
infile = initdat.infile
errfile = initdat.errfile
outfile = initdat.outfile
compfile = initdat.compfile
initfile = initdat.initfile
zguess = initdat.zinit
zfixB = initdat.zfixB
zfixR = initdat.zfixR
defaultdzstel = initdat.dzstel
defaultBdisp = initdat.disperseB
defaultRdisp = initdat.disperseR
defaultBvdisp = initdat.vdispB
defaultRvdisp = initdat.vdispR
startempfile = initdat.startempfile
strong = initdat.stronglines
fcninitpar = initdat.initparinfo
fcnlinefit = initdat.linefit

defaultncompB = 1
defaultncompR = 1

; Set up multiple components
if file_test(compfile) then begin
   readcol300,compfile,c_slit,c_ap,c_dv,c_br,$
              format='(I,I,D,A)',/silent,skip=1
endif
if file_test(initfile) then $
   readcol300,initfile,in_slit,in_ap,in_zinit,in_dzstel,$
              in_vdispB,in_vdispR,in_noemB,in_noemR,$
              in_masksig,in_floatfR,in_floatlR,in_floatsR,$
              in_floatfB,in_floatlB,in_floatsB,in_addmwstar,$
              format='(I,I,D,D,D,D,I,I,D,D,D,D,D,D,D,D)',$
              /silent,skip=1

; Set up emission-line masking parameters
defaultmasksig = 5d

; Get number of slits
fits_info,infile+'.fits',/silent,N_ext=next
next+=1
print,'Found ',next,' slits in ',infile,format='(/,A,I0,A,A,/)'
if n_elements(strong) eq 1 then strong = intarr(next)+strong
if n_elements(zguess) eq 1 then zguess = intarr(next)+zguess
if n_elements(zfixB) eq 1 then zfixB = intarr(next)+zfixB
if n_elements(zfixR) eq 1 then zfixR = intarr(next)+zfixR

; Cycle through slits
if ~ keyword_set(slitrange) then slitrange=[1,next] $
else if n_elements(slitrange) eq 1 then slitrange = [slitrange,slitrange]
slitrange = fix(slitrange)
for i=slitrange[0]-1,slitrange[1]-1 do begin
   
   print,'Processing slit ',i+1,' of ',next,'...',format='(A,I0,A,I0,A)'

;  Get linelist
   linelist = lris_initlinelist(strong=strong[i])

   spec = readspec(infile+'.fits',extension=i)
   specerr = readspec(errfile+'.fits',extension=i)
   naps = n_elements(spec[0,*])-1

;  Cycle through apertures
   if ~ keyword_set(aprange) then loc_aprange=[1,naps] $
   else if n_elements(aprange) eq 1 then loc_aprange=[aprange,aprange] $
   else loc_aprange = aprange
   loc_aprange = fix(loc_aprange)
   for j=loc_aprange[0],loc_aprange[1] do begin

;     Common initializations
      wave = spec[*,0]
      flux = spec[*,j]
      err = specerr[*,j]

      print,'  Processing aperture ',j,' of ',naps,'...',$
            format='(A,I0,A,I0,A)'
      slaplab='sl'+string(i+1,format='(I02)')+'_ap'+$
              string(j,format='(I02)')
      outfileslap = outfile + '_' + slaplab

;     Total number of velocity components.  Default is 1 component.
      ncompB=defaultncompB
      ncompR=defaultncompR
;     Initial set of redshifts and estimated shifts in redshift space
;     for extra components.
      vdispB = defaultBvdisp
      vdispR = defaultRvdisp
      masksig = defaultmasksig
      tmp = zguess[i]
      zinit = double(tmp[0])
      dzstel = defaultdzstel
      disperseB = defaultBdisp
      disperseR = defaultRdisp
      floatcompB = 0
      floatcompR = 0
      addmwstar = 0
      if file_test(initfile) then begin
         in_ind = where(in_slit eq i+1 AND in_ap eq j,count)
         if count eq 1 then begin
            if in_noemB[in_ind] ne -1 then ncompB = 0
            if in_noemR[in_ind] ne -1 then ncompR = 0
            if in_zinit[in_ind] ne -1 then begin
               tmp = in_zinit[in_ind]
               zinit = double(tmp[0])
            endif
            if in_dzstel[in_ind] ne -1 then begin
               tmp = in_dzstel[in_ind]
               dzstel = double(tmp[0])
            endif
            if in_vdispB[in_ind] ne -1 then begin
               disperseB = 1
               tmp = in_vdispB[in_ind]
               vdispB = double(tmp[0])
            endif
            if in_vdispR[in_ind] ne -1 then begin
               disperseR = 1
               tmp = in_vdispR[in_ind]
               vdispR = double(tmp[0])
            endif
            if in_masksig[in_ind] ne -1 then begin
               tmp = in_masksig[in_ind]
               masksig = double(tmp[0])
            endif
            if in_floatfR[in_ind] ne -1 AND $
               in_floatlR[in_ind] ne -1 AND $
               in_floatsR[in_ind] ne -1 then begin
               floatcompR = dblarr(3)
               floatcompR[0] = in_floatfR[in_ind]
               floatcompR[1] = in_floatlR[in_ind]
               floatcompR[2] = in_floatsR[in_ind]
            endif
            if in_floatfB[in_ind] ne -1 AND $
               in_floatlB[in_ind] ne -1 AND $
               in_floatsB[in_ind] ne -1 then begin
               floatcompB = dblarr(3)
               floatcompB[0] = in_floatfB[in_ind]
               floatcompB[1] = in_floatlB[in_ind]
               floatcompB[2] = in_floatsB[in_ind]
            endif
            if in_addmwstar[in_ind] ne -1 then begin
               tmp = in_addmwstar[in_ind]
               addmwstar = double(tmp[0])
            endif
         endif
      endif
      if ncompR gt 0 then begin
         zinit = refinez(zinit,wave/(1+zinit),flux)
         if ~ quiet then print,'Redshift adjusted from ',zguess[i],' to ',zinit,$
                               format='(A0,D0.5,A0,D0.5)'
      endif
      zred = zinit
      zblue = zinit
;   Case of multiple components
      if file_test(compfile) then begin
         mc_ind = where(c_slit eq i+1 AND c_ap eq j,count)
         if count ge 1 then begin
            dzcomp = c_dv[mc_ind]/299792d
            brtmp = c_br[mc_ind]
            if brtmp[0] eq 'b' or brtmp[0] eq 'br' then begin
               ncompB+=count
               zblue = [zinit,zinit + dzcomp]
            endif
            if brtmp[0] eq 'r' or brtmp[0] eq 'br' then begin
               ncompR+=count
               zred = [zinit,zinit + dzcomp]
            endif
         endif
      endif
;     Initialize redshift structure
      z = {star:zinit+dzstel,$
           gas: zred,$
           blue: zblue,$
           bluetilt:0d,$
           redtilt:0d,$
           blueint:blueint,$
           redint:redint}
           

;=======================      
;     Fit blue side
;=======================

      fcncontfit='fit_continuum'
      argsinitpar={tielambda:'Hbeta',$
                   bluered:'blue',$
                   floatcomp:floatcompB,$
                   zfix:zfixB[i]}
      argsz={bluered:'blue'}
      argspoly={addexp:1}
      argslinefit={floatcomp:floatcompB}
      fitran = lrisbr_redshift_spec(fitranBrest,z,bluered='blue')


;     Set masking of floating component properly
      if ncompB gt 0 then begin
;     Initialize emission line list
         linewave = linelist.wave
         linelabel = linelist.label
         nlines = n_elements(linewave)
;     Redshift emission line list
         linewavez = dblarr(nlines,ncompB)
         for ic=0,ncompB-1 do $
            linewavez[*,ic] = $
            call_function(fcnz,linewave,z,$
                          _extra={bluered:'blue',icomp:ic,gas:1})
         masklines = reform(linewavez,ncompB*nlines)
         meansig = vdispB/299792d*6500d
         maskwidths = masksig*replicate(meansig,nlines*ncompB)
         if floatcompB[0] ne 0 then begin
            masklines = [masklines,replicate(floatcompB[1],ncompB)]
            maskwidths = [maskwidths,replicate(2d*floatcompB[2],ncompB)]
         endif
      endif else begin
         masklines = 0
         maskwidths = 0
      endelse

      structinit = fit_spectrum(wave,flux,err,startempfile,z,ncompB,$
                                /subtract,$
                                linelist=linelist,quiet=quiet,$
                                obj_id=slaplab,$
                                time=time,fitran=fitran,$
                                fcninitpar=fcninitpar,$
                                argsinitpar=argsinitpar,$
                                disperse=disperseB,fcncontfit=fcncontfit,$
                                vdisp=vdispB,fcnz=fcnz,argsz=argsz,$
                                fcnlinefit=fcnlinefit,$
                                argslinefit=argslinefit,$
                                argspoly=argspoly,masklines=masklines,$
                                maskwidths=maskwidths,$
                                addmwstar=addmwstar)


      if ~ quiet then print,'FIT STATUS: ',structinit.fitstatus

      if ncompB gt 0 then begin

         lrisbr_updatez,structinit.param,z,'blue'
         z.star = z.blue[0] + dzstel
         linepars = sepfitpars(structinit.param,structinit.perror)
         masklines = reform(linepars.wave[*,0:ncompB-1],$
                            n_elements(linepars.wave[*,0:ncompB-1]))
         maskwidths = masksig*reform(linepars.sigma[*,0:ncompB-1],$
                                     n_elements(linepars.sigma[*,0:ncompB-1]))
         fitran = lrisbr_redshift_spec(fitranBrest,z,bluered='blue')
         argsz={bluered:'blue'}
         if floatcompB[0] ne 0 then begin
            floatcompB[0] = structinit.param[6]
            floatcompB[1] = structinit.param[7]
            floatcompB[2] = structinit.param[8]
            argsinitpar.floatcomp = floatcompB
            argslinefit.floatcomp = floatcompB
            masklines = [masklines,replicate(structinit.param[7],ncompB)]
            maskwidths = [maskwidths,replicate(2d*structinit.param[8],ncompB)]
         endif


         structblue = fit_spectrum(wave,flux,err,startempfile,z,ncompB,$
                                   /subtract,$
                                   linelist=linelist,quiet=quiet,$
                                   obj_id=slaplab,$
                                   time=time,fitran=fitran,$
                                   fcninitpar=fcninitpar,$
                                   argsinitpar=argsinitpar,$
                                   maskwidths=maskwidths,$
                                   masklines=masklines,$
                                   sigguess=linepars.sigma,$
                                   peakguess=linepars.fluxpk,$
                                   disperse=disperseB,$
                                   fcncontfit=fcncontfit,$
                                   vdisp=vdispB,fcnz=fcnz,argsz=argsz,$
                                   fcnlinefit=fcnlinefit,$
                                   argslinefit=argslinefit,$
                                   argspoly=argspoly,$
                                   addmwstar=addmwstar)
         

         if ~ quiet then print,'FIT STATUS: ',structblue.fitstatus
         lrisbr_updatez,structblue.param,z,'blue'

      endif else begin

         structblue = structinit

      endelse

      structblue = jjadd_tag(structblue,'z',z)
      structblue = jjadd_tag(structblue,'fitrange',fitran,/array)
      save,structblue,file=outfileslap+'_B.xdr'

;=======================      
;     Fit red side
;=======================

      fitran = lrisbr_redshift_spec(fitranRrest,z,bluered='red')
      if wave[n_elements(wave)-1] gt fitran[0] then begin

         fcncontfit='lrisbr_fit_redcont'
         argscontfit={ctBcoeff:structblue.ct_coeff}
         argsinitpar={tielambda:'Halpha',$
                      bluered:'red',$
                      floatcomp:floatcompR,$
                      zfix:zfixR[i]}
         argslinefit={floatcomp:floatcompR}
         argsz={bluered:'red'}
;      gd_indx = where(wave lt 6855 OR wave gt 6945)
;      wave = wave[gd_indx]
;      flux = flux[gd_indx]
;      err = err[gd_indx]

;     Set masking of floating component properly
         if ncompR gt 0 then begin
;     Initialize emission line list
            linewave = linelist.wave
            linelabel = linelist.label
            nlines = n_elements(linewave)
;     Redshift emission line list
            linewavez = dblarr(nlines,ncompR)
            for ic=0,ncompR-1 do $
               linewavez[*,ic] = $
               call_function(fcnz,linewave,z,$
                             _extra={bluered:'red',icomp:ic,gas:1})
            masklines = reform(linewavez,ncompR*nlines)
            meansig = vdispR/299792d*6500d
            maskwidths = masksig*replicate(meansig,nlines*ncompR)
            if floatcompR[0] ne 0 then begin
               masklines = [masklines,replicate(floatcompR[1],ncompR)]
               maskwidths = [maskwidths,replicate(2d*floatcompR[2],ncompR)]
            endif
         endif else begin
            masklines = 0
            maskwidths = 0
         endelse

         
         structinit = fit_spectrum(wave,flux,err,startempfile,z,ncompR,$
                                   linelist=linelist,quiet=quiet,$
                                   obj_id=slaplab,$
                                   time=time,fitran=fitran,/subtract,$
                                   fcncontfit=fcncontfit,$
                                   argscontfit=argscontfit,$
                                   fcninitpar=fcninitpar,$
                                   argsinitpar=argsinitpar,$
                                   disperse=disperseR,vdisp=vdispR,$
                                   fcnz=fcnz,argsz=argsz,$
                                   fcnlinefit=fcnlinefit,$
                                   argslinefit=argslinefit,$
                                   masklines=masklines,$
                                   maskwidths=maskwidths,$
                                   argspoly=argspoly,$
                                   addmwstar=addmwstar)
         


         if ~ quiet then print,'FIT STATUS: ',structinit.fitstatus

         if ncompR gt 0 then begin
            
            lrisbr_updatez,structinit.param,z,'red'
            z.star = z.gas[0] + dzstel
            linepars = sepfitpars(structinit.param,structinit.perror)
            masklines = reform(linepars.wave[*,0:ncompR-1],$
                               n_elements(linepars.wave[*,0:ncompR-1]))
            maskwidths = masksig*reform(linepars.sigma[*,0:ncompR-1],$
                                        n_elements(linepars.sigma[*,0:ncompR-1]))
            fitran = lrisbr_redshift_spec(fitranRrest,z,bluered='red')
            argsz={bluered:'red'}
            if floatcompR[0] ne 0 then begin
               floatcompR[0] = structinit.param[6]
               floatcompR[1] = structinit.param[7]
               floatcompR[2] = structinit.param[8]
               argsinitpar.floatcomp = floatcompR
               argslinefit.floatcomp = floatcompR
               masklines = [masklines,replicate(structinit.param[7],ncompR)]
               maskwidths = [maskwidths,replicate(2d*structinit.param[8],ncompR)]
            endif


            structred = fit_spectrum(wave,flux,err,startempfile,z,ncompR,$
                                     linelist=linelist,quiet=quiet,$
                                     obj_id=slaplab,$
                                     time=time,fitran=fitran,/subtract,$
                                     fcncontfit=fcncontfit,$
                                     argscontfit=argscontfit,$
                                     fcninitpar=fcninitpar,$
                                     argsinitpar=argsinitpar,$
                                     sigguess=linepars.sigma,$
                                     peakguess=linepars.fluxpk,$
                                     disperse=disperseR,vdisp=vdispR,$
                                     fcnz=fcnz,argsz=argsz,$
                                     fcnlinefit=fcnlinefit,$
                                     argslinefit=argslinefit,$
                                     masklines=masklines,$
                                     maskwidths=maskwidths,$
                                     argspoly=argspoly,$
                                     addmwstar=addmwstar)


            if ~ quiet then print,'FIT STATUS: ',structred.fitstatus
            lrisbr_updatez,structred.param,z,'red'
            
         endif else begin

            structred = structinit

         endelse

         structred = jjadd_tag(structred,'z',z)
         structred = jjadd_tag(structred,'fitrange',fitran,/array)
         save,structred,file=outfileslap+'_R.xdr'

      endif else begin

         print,'No red-side fit: redshift too high.'

      endelse

   endfor
   
endfor

print,'Total time for calculation: ',systime(1)-starttime,' s.',$
      format='(/,A0,I0,A0,/)'

end
