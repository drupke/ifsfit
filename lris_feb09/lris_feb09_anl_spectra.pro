pro lris_feb09_anl_spectra,lab,mask,slitrange=slitrange,$
                           aprange=aprange,$
                           nobcor=nobcor,showplot=showplot
;+
; NAME:
;     LRIS_FEB09_ANL_SPECTRA
;
; PURPOSE:
;     Wrapper script for analyzing fitted LRIS spectra.
;
; EXPLANATION:
;
; CALLING SEQUENCE
;     lris_anl_spectra,mask,slitrange=,aprange=,/dobcor,/showplot
;
; INPUTS:
;     mask - String name of mask, as specified in LRIS_INITFIT_SPECTRA
;     slitrange - scalar or two-element array specifying 
;                 range of slits to fit
;     aprange - scalar or two-element array specifying 
;               range of apertures to fit
;     nobcor - do not correct blue line fluxes based on Balmer lines
;     showplot - display to screen each fit as it is plotted
;
; INPUT:
;     A binary (.genx) file containing the results of the fit.
;
; OUTPUT: 
;     JPG plots of stellar, emission-line, and continuum fits, and
;     data files summarizing outputs.
;
; METHOD:
;
; REVISION HISTORY:
;     09may13  David Rupke  created
;     09jun07  DSNR         added multiple components
;     10apr21  DSNR         converted to observed-frame fits
;-

  if keyword_set(showplot) then zbuf=0 else zbuf=1

  fwhmtosig = 2d*sqrt(2d*alog(2d))

; Get parameters to initialize fit
  initdat = lris_feb09_initfit_spectra(lab,mask)
  infile = initdat.infile
  errfile = initdat.errfile
  outfile = initdat.outfile
  compfile = initdat.compfile
  stitchwave = initdat.stitchwave

; Get linelist
  linelist = lris_initlinelist()
; which lines to output
  outlines=['[OII]3726','[OII]3729','Hbeta','[OIII]5007','[OI]6300',$
            'Halpha','[NII]6583','[SII]6716','[SII]6731']
  lris_printlinepars,[''],[0],[0],outfile+'.lines.dat','',$
                     /init,whichlines=outlines

; Get number of slits
  fits_info,infile+'.fits',/silent,N_ext=next
  next+=1
  print,'Found ',next,' slits in ',infile,format='(/,A,I0,A,A,/)'
  
; Open output file for fit results
  openw,fitunit,outfile+'.fit.dat',/get_lun
  printf,fitunit,'# Slit/Ap','comp','Rchi2','Niter','fwhmB','fwhmR',$
         'zR','zstel','zR-zB','dzdlB','fwhmRbroad','fwhmRred',$
         'dlRred(A)','fRbroad(%)','fRred(%)','fRfloat','lRfloat',$
         'fwhmRfloat',$
         format='(A-10,A6,A8,A6,A8,A8,A9,A9,A9,A9,A11,A11,A11,A11,'+$
         'A11,A11,A11,A11)'
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
  else if n_elements(slitrange) eq 1 then slitrange=[slitrange,slitrange]
  slitrange = fix(slitrange)
  for i=slitrange[0]-1,slitrange[1]-1 do begin

     print,'Processing slit ',i+1,' of ',next,'...',$
           format='(A,I0,A,I0,A)'

;   Read in slit
     spec = readspec(infile+'.fits',extension=i)
     naps = n_elements(spec[0,*])-1

;   Cycle through apertures
     if ~ keyword_set(aprange) then loc_aprange=[1,naps] $
     else if n_elements(aprange) eq 1 then loc_aprange=[aprange,aprange] $
     else loc_aprange = aprange
     loc_aprange = fix(loc_aprange)
     for j=loc_aprange[0],loc_aprange[1] do begin

        print,'  Processing aperture ',j,' of ',naps,'...',$
              format='(A,I0,A,I0,A)'
        outfileslap = outfile + '_sl' + string(i+1,format='(I02)') + $
                      '_ap' + string(j,format='(I02)')

        if (i eq 0 AND j eq 1) then append=0 else append=1

;     Restore fit from binary file
        restore,file=outfileslap+'.xdr'

;     Rest wavelength range in data
        waveran = [struct.wave[0],struct.wave[n_elements(struct.wave)-1]]

;     Number of velocity components
        ndpar = size(struct.param)
        if ndpar[0] eq 0 then ncomp = 0 $
        else ncomp = struct.param[1]


        if ncomp gt 0 then begin

           fluxpkerr = lris_fluxerrors(struct)

;     Separate best-fit parameters
           linepars = lris_feb09_sepfitpars(struct.param,struct.perror,$
                                            stitchwave,waveran=waveran,$
                                            fluxpkerr=fluxpkerr)
           nlines = n_elements(linepars.flux[*,0])

;     Plot results
           if ~ keyword_set(noplot) then begin
              plotstelfit,struct,outfileslap+'_stel',zbuf=zbuf
              plotstronglines,struct,outfileslap+'_strg',zbuf=zbuf
              plotbalmlines,struct,outfileslap+'_balm',zbuf=zbuf
           endif

;     Print fit parameters to a text file
           ihalpha = where(struct.linelabel eq 'Halpha')
           ishalpha = struct.param[0] + ihalpha*3 + 2
           gint_br = gaussarea(1d/(2d*struct.param[9]^2d))
           gint_rd = gaussarea(1d/(2d*struct.param[12]^2d))
           gint_ha = gaussarea(1d/(2d*struct.param[ishalpha]^2d))
           pct_br = 100d*struct.param[8]*gint_br.area/gint_ha.area
           pct_rd = 100d*struct.param[11]*gint_rd.area/gint_ha.area
           printf,fitunit,struct.obj_id,$
                  1,$
                  struct.redchisq,$
                  struct.niter,$
                  fwhmtosig*linepars.sigma[0],$
                  fwhmtosig*linepars.sigma[nlines-1],$
                  struct.z.gas[0],$
                  struct.z.star,$
                  struct.z.gas[0]-struct.z.blue[0],$
                  struct.z.bluetilt,$
                  fwhmtosig*struct.param[9],$
                  fwhmtosig*struct.param[12],$
                  struct.param[10],$
                  pct_br,pct_rd,$
                  struct.param[13],$
                  struct.param[14],$
                  fwhmtosig*struct.param[15],$
                  format='(A-10,I6.0,D8.2,I6.0,D8.2,D8.2,D9.5,D9.5'+$
                  ',D9.5,E9.1,D11.2,D11.2,D11.1,D11.1,D11.1,E11.2,'+$
                  'D11.1,D11.1)'
           if ncomp gt 1 then $
              for k=1,ncomp-1 do $
                 printf,fitunit,struct.obj_id,$
                        k+1,$
                        struct.redchisq,$
                        struct.niter,$
                        fwhmtosig*linepars.sigma[0,k],$
                        fwhmtosig*linepars.sigma[nlines-1,k],$
                        struct.z.gas[k],$
                        struct.z.star,$
                        struct.z.gas[0]-struct.z.blue[0],$
                        struct.z.bluetilt,$
                        fwhmtosig*struct.param[9],$
                        fwhmtosig*struct.param[12],$
                        struct.param[10],$
                        pct_br,pct_rd,$
                        struct.param[13],$
                        struct.param[14],$
                        fwhmtosig*struct.param[15],$
                        format='(A-10,I6.0,D8.2,I6.0,D8.2,D8.2,'+$
                        'D9.5,D9.5,D9.5,E9.1,D11.2,D11.2,D11.1,'+$
                        'D11.1,D11.1,E11.2,D11.1,D11.1)'
         
;     Use Balmer lines to check flux calibration
           bfluxes = getbalmfluxes(struct.linelabel,linepars.flux,$
                                   errors=linepars.fluxerr)
           bcorr = lris_chkbalmlines(bfluxes,ebmv=ebmv)
           printf,balmunit,struct.obj_id,bcorr[0,0,0],bcorr[1,0,0],$
                  bcorr[2,0,0],bcorr[3,0,0],$
                  ebmv[0,0,0],ebmv[1,0,0],ebmv[2,0,0],ebmv[3,0,0],$
                  format='(A-10,D8.2,D8.2,D8.2,D8.2,D8.2,D8.2,D8.2,D8.2)'
           printf,ebalmunit,struct.obj_id,bcorr[0,1,0],bcorr[1,1,0],$
                  bcorr[2,1,0],bcorr[3,1,0],$
                  ebmv[0,1,0],ebmv[1,1,0],ebmv[2,1,0],ebmv[3,1,0],$
                  format='(A-10,D8.2,D8.2,D8.2,D8.2,D8.2,D8.2,D8.2,D8.2)'
         
;     Correct line fluxes using Balmer decrement
           if ~ keyword_set(nobcor) then begin
              nobcor=0
              bluelines = where(struct.linewave*(1+struct.z.star) lt $
                                stitchwave)
              if bcorr[1,0,0] lt 0 then begin
                 usebcorr = 1d
                 usebcorrerr = 0d
              endif else begin
                 usebcorr = bcorr[1,0,0]
                 usebcorrerr = bcorr[1,1,0]
              endelse
              linepars.flux[bluelines,*] *= usebcorr
;;               linepars.fluxerr[bluelines,*] = $
;;                  linepars.flux[bluelines,*] * $
;;                  sqrt((linepars.fluxerr[bluelines,*]/ $
;;                        linepars.flux[bluelines,*])^2d + $
;;                       (usebcorrerr/usebcorr)^2d)
;;               nans = where(finite(linepars.fluxerr) eq 0,ctnan)
;;               if ctnan gt 0 then linepars.fluxerr[nans] = 0d
           endif else nobcor=1
         
;     Print line fluxes to a text file
;     Initialize output linelist
           label = string('sl',i+1,'_ap',j,format='(A0,I0,A0,I0)')
           lris_printlinepars,struct.linelabel,linepars.flux,$
                              linepars.fluxerr,$
                              outfile+'.lines.dat',label,/append,$
                              whichlines=outlines,nobcor=nobcor

        endif else begin

           if ~ keyword_set(noplot) then $
              plotstelfit,struct,outfileslap+'_stel',zbuf=zbuf

           printf,fitunit,struct.obj_id,$
                  1,1d,1,0d,0d,$
                  struct.z.gas[0],$
                  struct.z.star,$
                  struct.z.gas[0]-struct.z.blue[0],$
                  struct.z.bluetilt,$
                  0d,0d,0d,0d,0d,0d,0d,0d,$
                  format='(A-10,I6.0,D8.2,I6.0,D8.2,D8.2,'+$
                  'D9.5,D9.5,D9.5,E9.1,'+$
                  'D11.2,D11.2,D11.1,D11.1,D11.1,E11.2,D11.1,D11.1)'

           label = string('sl',i+1,'_ap',j,format='(A0,I0,A0,I0)')
           lris_printlinepars,0d,0d,0d,outfile+'.lines.dat',$
                              label,/append,whichlines=outlines,$
                              /zeroflux
           
        endelse

     endfor

  endfor

  free_lun,fitunit
  free_lun,balmunit
  free_lun,ebalmunit

end
