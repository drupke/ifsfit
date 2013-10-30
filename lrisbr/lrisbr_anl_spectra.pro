pro lrisbr_anl_spectra,lab,mask,slitrange=slitrange,aprange=aprange,$
                       nobcor=nobcor,showplot=showplot
;+
; NAME:
;     LRISBR_ANL_SPECTRA
;
; PURPOSE:
;     Wrapper script for analyzing fitted LRIS spectra.
;
; EXPLANATION:
;
; CALLING SEQUENCE
;     lrisbr_anl_spectra,mask,slitrange=,aprange=,/dobcor,/showplot
;
; INPUTS:
;     mask - String name of mask, as specified in LRIS_INITFIT_SPECTRA
;     slitrange - scalar or two-element array specifying range of slits to fit
;     aprange - scalar or two-element array specifying range of apertures to fit
;     nobcor - do not correct blue line fluxes based on Balmer lines
;     showplot - display to screen each fit as it is plotted
;
; INPUT:
;     A binary (.genx) file containing the results of the fit.
;
; OUTPUT:
;     JPG plots of stellar, emission-line, and continuum fits, and data files
;     summarizing outputs.
;
; METHOD:
;
; REVISION HISTORY:
;     09may13  David Rupke  created
;     09jun07  David Rupke  added multiple components
;-

  if keyword_set(showplot) then zbuf=0 else zbuf=1
  if keyword_set(nobcor) then dobcor=0 else dobcor=1

  fwhmtosig = 2d*sqrt(2d*alog(2d))

; Get parameters to initialize fit
  initdat = lrisbr_initfit_spectra(lab,mask)
  infile = initdat.infile
  errfile = initdat.errfile
  outfile = initdat.outfile
  compfile = initdat.compfile

; Get linelist
  linelist = lris_initlinelist()
; which lines to output
  outlines=['[OII]3726','[OII]3729','Hbeta','[OIII]5007','[OI]6300','Halpha',$
            '[NII]6583','[SII]6716','[SII]6731']
  lris_printlinepars,[''],[0],[0],outfile+'.lines.dat','',/init,whichlines=outlines

; Get number of slits
  fits_info,infile+'.fits',/silent,N_ext=next
  next+=1
  print,'Found ',next,' slits in ',infile,format='(/,A,I0,A,A,/)'
  
; Open output file for fit results
  openw,fitunit,outfile+'.fit.dat',/get_lun
  printf,fitunit,'# Slit/Ap','comp','RchiB','RchiR','fwhmB','fwhmR',$
         'zB','zR','zstelB','zstelR','ztiltB','ztiltR',$
         format='(A-10,A6,A8,A8,A8,A8,A9,A9,A9,A9,A9,A9)'
; Open output file for Balmer line results
  openw,balmunit,outfile+'.balmer.dat',/get_lun
  printf,balmunit,'# ','Blue side multiplier','E(B-V)',$ 
         format='(A-10,A-32,A-32)'
  printf,balmunit,'# ','-----------------------------  ',$
         '-----------------------------  ',$
         format='(A-10,A-32,A-32)'
  printf,balmunit,'# Slit/Ap','Ha/Hb','Hg/Hb','Hd/Hb','He/Hb',$
         'Ha/Hb','Hg/Hb','Hd/Hb','He/Hb',$
         format='(A-10,A8,A8,A8,A8,A8,A8,A8,A8)'
; 
  openw,ebalmunit,outfile+'.balmerr.dat',/get_lun
  printf,ebalmunit,'# ','Blue side multiplier','E(B-V)',$ 
         format='(A-10,A-32,A-32)'
  printf,ebalmunit,'# ','-----------------------------  ',$
         '-----------------------------  ',$
         format='(A-10,A-32,A-32)'
  printf,ebalmunit,'# Slit/Ap','Ha/Hb','Hg/Hb','Hd/Hb','He/Hb',$
         'Ha/Hb','Hg/Hb','Hd/Hb','He/Hb',$
         format='(A-10,A8,A8,A8,A8,A8,A8,A8,A8)'

; Cycle through slits
  if ~ keyword_set(slitrange) then slitrange=[1,next] $
  else if n_elements(slitrange) eq 1 then slitrange = [slitrange,slitrange]
  slitrange = fix(slitrange)
  for i=slitrange[0]-1,slitrange[1]-1 do begin

     print,'Processing slit ',i+1,' of ',next,'...',format='(A,I0,A,I0,A)'

;   Read in slit
     spec = readspec(infile+'.fits',extension=i)
     naps = n_elements(spec[0,*])-1

;   Cycle through apertures
     if ~ keyword_set(aprange) then loc_aprange=[1,naps] $
     else if n_elements(aprange) eq 1 then loc_aprange=[aprange,aprange] $
     else loc_aprange = aprange
     loc_aprange = fix(loc_aprange)
     for j=loc_aprange[0],loc_aprange[1] do begin

        print,'  Processing aperture ',j,' of ',naps,'...',format='(A,I0,A,I0,A)'
       
        outfileslap = outfile + '_sl' + string(i+1,format='(I02)') + '_ap' + $
                      string(j,format='(I02)')
        
        if (i eq 0 AND j eq 1) then append=0 else append=1
        
;      Restore fit from binary file
        restore,file=outfileslap+'_B.xdr'
        dored=1
        if file_test(outfileslap+'_R.xdr') $
        then restore,file=outfileslap+'_R.xdr' $
        else dored=0

        if dored then wave = [structblue.wave,structred.wave] $
        else wave = structblue.wave
        
;      Wavelength range in data
        waveranB = [structblue.wave[0],$
                    structblue.wave[n_elements(structblue.wave)-1]]
        if dored then $
           waveranR = [structred.wave[0],$
                       structred.wave[n_elements(structred.wave)-1]]
        
;      Number of velocity components
        if structblue.param[0] eq 0 then ncompB = 0 $
        else ncompB = structblue.param[1]
        ncompR = 0
        if dored then $
           if structred.param[0] ne 0 then $
              ncompR = structred.param[1]
        
;      Separate best-fit parameters
        nlines=0
        if ncompB gt 0 then begin
           fluxpkerrB = lris_fluxerrors(structblue)
           lineparsB = sepfitpars(structblue.param,structblue.perror,$
                                  waveran=waveranB,fluxpkerr=fluxpkerrB)
           nlines +=  n_elements(lineparsB.flux[*,0])
        endif
        if ncompR gt 0 then begin
           fluxpkerrR = lris_fluxerrors(structred)
           lineparsR = sepfitpars(structred.param,structred.perror,$
                                  waveran=waveranR,fluxpkerr=fluxpkerrR)
           nlines += n_elements(lineparsR.flux[*,0])
        endif
        
;      Plot results
        if ~ keyword_set(noplot) then begin
           if dored then $
              structBR = {wave: [structblue.wave,structred.wave], $
                          spec: [structblue.spec,structred.spec], $
                          specfit: [structblue.specfit,$
                                    structred.specfit], $
                          spec_nocnt: [structblue.spec_nocnt,$
                                       structred.spec_nocnt], $
                          fitrange: [structblue.fitrange[0],structred.fitrange[1]], $
                          z: structred.z $
                         } $
           else $
              structBR = {wave: structblue.wave, $
                          spec: structblue.spec, $
                          specfit: structblue.specfit, $
                          spec_nocnt: structblue.spec_nocnt, $
                          fitrange: structblue.fitrange, $
                          z: structblue.z $
                         }
           plotstelfit,structBR,outfileslap+'_stel',zbuf=zbuf
           if ncompB gt 0 OR ncompR gt 0 then begin
              plotstronglines,structBR,outfileslap+'_strg',zbuf=zbuf
              plotbalmlines,structBR,outfileslap+'_balm',zbuf=zbuf
           endif
        endif
        
;      Print fit parameters to a text file
        fwhmB=0d
        fwhmR=0d
        if ncompB gt 0 then $
           fwhmB = fwhmtosig*lineparsB.sigma[0]
        if ncompR gt 0 then $
           fwhmR = fwhmtosig*lineparsR.sigma[0]
        if dored then $
           printf,fitunit,structblue.obj_id,$
                  1,$
                  structblue.redchisq,$
                  structred.redchisq,$
                  fwhmB,$
                  fwhmR,$
                  structred.z.blue[0],$
                  structred.z.gas[0],$
                  structblue.z.star,$
                  structred.z.star,$
                  structblue.z.bluetilt,$
                  structred.z.redtilt,$
                  format='(A-10,I6.0,D8.2,D8.2,D8.2,D8.2,D9.5,D9.5,D9.5,D9.5,E9.1,E9.1)' $
        else $
           printf,fitunit,structblue.obj_id,$
                  1,$
                  structblue.redchisq,$
                  0d,$
                  fwhmB,$
                  fwhmR,$
                  structblue.z.blue[0],$
                  0d,$
                  structblue.z.star,$
                  0d,$
                  structblue.z.bluetilt,$
                  0d,$
                  format='(A-10,I6.0,D8.2,D8.2,D8.2,D8.2,D9.5,D9.5,D9.5,D9.5,E9.1,E9.1)'
        if ncompB gt 1 OR ncompR gt 1 then begin
           ncomp = max([ncompB,ncompR])
           for k=1,ncomp-1 do begin
              if ncompB ge k+1 then begin
                 redchisqB = structblue.redchisq
                 fwhmB = fwhmtosig*lineparsB.sigma[0,k]
                 zB = structblue.z.blue[k]
                 zBstar = structblue.z.star
                 zBtilt = structblue.z.bluetilt
              endif else begin
                 redchisqB = 0d
                 fwhmB = 0d
                 zB = 0d
                 zBstar = 0d
                 zBtilt = 0d
              endelse
              if ncompR ge k+1 then begin
                 redchisqR = structred.redchisq
                 fwhmR = fwhmtosig*lineparsR.sigma[0,k]
                 zR = structred.z.gas[k]
                 zRstar = structred.z.star
                 zRtilt = structred.z.redtilt
              endif else begin
                 redchisqR = 0d
                 fwhmR = 0d
                 zR = 0d
                 zRstar = 0d
                 zRtilt = 0d
              endelse
              printf,fitunit,structred.obj_id,$
                     k+1,$
                     redchisqB,$
                     redchisqR,$
                     fwhmB,$
                     fwhmR,$
                     zB,$
                     zR,$
                     zBstar,$
                     zRstar,$
                     zBtilt,$
                     zRtilt,$
                     format='(A-10,I6.0,D8.2,D8.2,D8.2,D8.2,D9.5,D9.5,D9.5,D9.5,E9.1,E9.1)'
           endfor
        endif

        if ncompB gt 0 AND ncompR gt 0 then begin

;     Use Balmer lines to check flux calibration
           if ncompB eq 1 AND ncompR eq 2 then begin
              lineparsBflux = [[lineparsB.flux],$
                               [lineparsB.flux-lineparsB.flux]]
              lineparsBfluxerr = [[lineparsB.fluxerr],$
                                  [lineparsB.fluxerr-lineparsB.fluxerr]]
           endif else begin
              lineparsBflux = lineparsB.flux
              lineparsBfluxerr = lineparsB.fluxerr
           endelse
           bfluxesblue = getbalmfluxes(structblue.linelabel,$
                                       lineparsBflux,$
                                       errors=lineparsBfluxerr)
           if dored then begin
              if ncompB eq 2 AND ncompR eq 1 then begin
                 lineparsRflux = [[lineparsR.flux],$
                                  [lineparsR.flux-lineparsR.flux]]
                 lineparsRfluxerr = [[lineparsR.fluxerr],$
                                     [lineparsR.fluxerr-lineparsR.fluxerr]]
              endif else begin
                 lineparsRflux = lineparsR.flux
                 lineparsRfluxerr = lineparsR.fluxerr
              endelse
              bfluxesred = getbalmfluxes(structred.linelabel,$
                                         lineparsRflux,$
                                         errors=lineparsRfluxerr)
              zbfluxesblue = where(bfluxesblue.value eq 0)
              bfluxesblue.value[zbfluxesblue] = $
                 bfluxesred.value[zbfluxesblue]
              bfluxesblue.error[zbfluxesblue] = $
                 bfluxesred.error[zbfluxesblue]
           endif

           if abs(ncompB-ncompR) ge 2 then begin
              print,'Balmer correction not performed properly for '+$
                    'this spectrum; components not treated right.'
           endif

           bcorr = lris_chkbalmlines(bfluxesblue,ebmv=ebmv)
           printf,balmunit,structblue.obj_id,bcorr[0,0,0],bcorr[1,0,0],$
                  bcorr[2,0,0],bcorr[3,0,0],$
                  ebmv[0,0,0],ebmv[1,0,0],ebmv[2,0,0],ebmv[3,0,0],$
                  format='(A-10,D8.2,D8.2,D8.2,D8.2,D8.2,D8.2,D8.2,D8.2)'
           printf,ebalmunit,structblue.obj_id,bcorr[0,1,0],bcorr[1,1,0],$
                  bcorr[2,1,0],bcorr[3,1,0],$
                  ebmv[0,1,0],ebmv[1,1,0],ebmv[2,1,0],ebmv[3,1,0],$
                  format='(A-10,D8.2,D8.2,D8.2,D8.2,D8.2,D8.2,D8.2,D8.2)'
        
;     Correct line fluxes using Balmer decrement
           if dobcor then begin
              if bcorr[1,0,0] lt 0.2d OR $
                 bcorr[1,0,0] gt 5d then begin
                 usebcorr = 1d
                 usebcorrerr = 0d
                 nobcor = 1
              endif else if bcorr[0,0,0] lt 1 AND $
                 bcorr[1,0,0] gt bcorr[0,0,0] then begin
                 usebcorr = bcorr[0,0,0]
                 usebcorrerr = bcorr[0,1,0]
                 nobcor=0
              endif else begin
                 usebcorr = bcorr[1,0,0]
                 usebcorrerr = bcorr[1,1,0]
                 nobcor=0
              endelse
              lineparsB.flux *= usebcorr
;;               lineparsB.fluxerr = lineparsB.flux * $
;;                                   sqrt((lineparsB.fluxerr $
;;                                         / lineparsB.flux)^2d + $
;;                                        (usebcorrerr/usebcorr)^2d)
;;               nans = where(finite(lineparsB.fluxerr) eq 0,ctnan)
;;               if ctnan gt 0 then lineparsB.fluxerr[nans] = 0d
           endif else nobcor=1

        endif
        
;     Print line fluxes to a text file
;     Initialize output linelist
        label = string('sl',i+1,'_ap',j,format='(A0,I0,A0,I0)')
        if ncompB gt 0 AND ncompR gt 0 then begin
           lrisbr_printlinepars,structblue.linelabel,$
                                lineparsB.flux,lineparsB.fluxerr,$
                                lineparsR.flux,lineparsR.fluxerr,$
                                outfile+'.lines.dat',label,/append,$
                                whichlines=outlines,nobcor=nobcor
        endif else if ncompB gt 0 then begin
           lrisbr_printlinepars,structblue.linelabel,$
                                lineparsB.flux,lineparsB.fluxerr,$
                                0,0,$
                                outfile+'.lines.dat',label,/append,$
                                whichlines=outlines,nobcor=nobcor
        endif else if ncompR gt 0 then begin
           lrisbr_printlinepars,structblue.linelabel,$
                                0,0,$
                                lineparsR.flux,lineparsR.fluxerr,$
                                outfile+'.lines.dat',label,/append,$
                                whichlines=outlines,nobcor=nobcor
        endif else begin
           lrisbr_printlinepars,structblue.linelabel,0,0,0,0,$
                                outfile+'.lines.dat',label,/append,$
                                whichlines=outlines,nobcor=nobcor
        endelse
        
     endfor

  endfor
  
  free_lun,fitunit
  free_lun,balmunit
  free_lun,ebalmunit
  
end
