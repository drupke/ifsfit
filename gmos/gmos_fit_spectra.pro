pro gmos_fit_spectra,time=time,verbose=verbose
;
; History
;  09jul08  DSNR  copied from LRIS routine to GMOS
;  10may27  DSNR  started re-write for new data
;
  
  starttime = systime(1)

; Set rest-frame fit range
  fitran = [6350d,7100d] ; KPNO
  vdisp=100
  fcninitpar='gmos_initparinfo'
  fcnlinefit='manygauss_slow'
  fcncontfit='gmos_fit_flat_continuum'
  mc_refwave = 6562.8d
  checkcomp_sigcut_init=5d
;  checkcomp_sigcut_init=1d ; for corner spectra
; Total number of velocity components.  Default is 1 component.
  ncomp_init=2
;  ncomp_init=1 ; for corner spectra

  if ~ keyword_set(time) then time=0
  if keyword_set(verbose) then quiet=0 else quiet=1

;;   startempfile = '/Users/drupke/src/idl/uhspecfit/stellar_models/'+$
;;                  'gonzalezdelgado/SSPGeneva_z008+020+040.sav'
  startempfile = '/Users/drupke/src/idl/uhspecfit/stellar_models/'+$
                 'bruzualcharlot/bc03_tremonti.sav'
; For GMOS data
;;   rootdir = '/Users/drupke/winds/gmos/f10565/'
;;   specpre = 'spec_'
;;   errpre = 'err_'
;;   outpre = 'f10565p2448_'
;;   startrow=1
;;   nrow=11
;;   startcol=1
;;   ncol=13
;;   startrow=6
;;   nrow=1
;;   startcol=9
;;   ncol=1
; for corner spectra
;;   startrow=101
;;   nrow=2
;;   startcol=101
;;   ncol=2
; For KPNO longslit spectra
  rootdir = '/Users/drupke/winds/gmos/f10565/kpnolong/'
  specpre = 'f10565p2448.02dec28.'
  outpre = 'f10565p2448_'
  startrow=1
  nrow=9
  startcol=1
  ncol=1

; Get linelist
  linelist = gmos_initlinelist()

  for i=startrow,startrow+nrow-1 do begin

     for j=startcol,startcol+ncol-1 do begin

        ncomp = ncomp_init
        checkcomp_sigcut = checkcomp_sigcut_init
      
        print,'Fitting spectrum ',i,', ',j,format='(A,I0,A,I0)'

        zguess = 0.043d
;       GMOS
;;        lab = string(i,'_',j,format='(I03,A,I03)')
;;         specfile = rootdir+'spec/'+specpre+lab+'.fits'
;;         errfile = rootdir+'err/'+errpre+lab+'.fits'
;;         outfile = rootdir+'specfit/'+outpre+lab+'.genx'
;       KPNO
        lab = string(i,format='(I0)')
        specfile = rootdir+'red/'+specpre+lab+'.fits'
        errfile = rootdir+'red/'+specpre+lab+'.err.fits'
        outfile = rootdir+'specfit/'+outpre+lab+'.genx'
        spec = readspec(specfile)
        specerr  = readspec(errfile)

        if (i eq 1 AND j eq 1) then append=0 else append=1

fit:

;       Estimated shifts in redshift space for components.  0d is
;       default for a single component.
        if ncomp eq 2 then mcomp=[0d,-5d/mc_refwave] $
        else if ncomp eq 1 then mcomp=[0d] $
        else mcomp=-1
;       Initial set of redshifts
        if ncomp gt 0 then zarr = zguess + mcomp else zarr = zguess

        wave = spec[*,0]
        flux = spec[*,1]
        err = specerr[*,1]
        
;       Create de-redshifted wavelengths
        waverest = wave / (1d + zarr[0])
;       Remove NaI D line
        nad_indx = where(waverest lt 5870 OR waverest gt 5900)
;       Remove telluric line
        tel_indx = where(wave gt 6855 AND wave lt 6945)
        gd_indx = cmset_op(nad_indx,'AND',/not2,tel_indx)

        waverest = waverest[gd_indx]
        flux = flux[gd_indx]
        err = err[gd_indx]

;       First fit
        structinit = fit_spectrum(waverest,flux,err,startempfile,/subtract,$
                                  vdisp=vdisp,linelist=linelist,quiet=quiet,$
                                  time=time,mcomp=mcomp,fitran=fitran,$
                                  obj_id=lab,fcninitpar=fcninitpar,$
                                  fcnlinefit=fcnlinefit,$
                                  fcncontfit=fcncontfit,/disperse)

;       Repeat fit with new z estimate and better estimate for line
;       widths, if we're fitting emission lines.

        if mcomp[0] ne -1 then begin

;          Make sure components are ordered such that component with
;          highest peak flux is the first ("reference") component
           gmos_orderlines,structinit,zarr

;          Update de-redshifted wavelengths and z array
           gmos_updatez,structinit.param,waverest,zarr
           mcomp = zarr - zarr[0]

;          Update linelist so that maskwidths array works OK
           newlinelist = {wave:structinit.linewave,label:structinit.linelabel}
;          Estimate linewidths
;          Mask at MASKSIG # of sigma away from line in either direction
           masksig = 5d
           linepars = sepfitpars(structinit.param,structinit.perror)
           masklines = linepars.wave[*,0]
           maskwidths = masksig*linepars.sigma[*,0]
           for k=1,ncomp-1 do begin
              masklines = [masklines,linepars.wave[*,k]]
              maskwidths = [maskwidths,masksig*linepars.sigma[*,k]]
           endfor

           struct = fit_spectrum(waverest,flux,err,startempfile,$
                                 vdisp=vdisp,/subtract,$
                                 linelist=newlinelist,quiet=quiet,obj_id=lab,$
                                 time=time,mcomp=mcomp,fitran=fitran,$
                                 fcninitpar=fcninitpar,maskwidths=maskwidths,$
                                 masklines=masklines,sigguess=linepars.sigma,$
                                 peakguess=linepars.fluxpk,/disperse,$
                                 fcnlinefit=fcnlinefit,fcncontfit=fcncontfit)

           gmos_orderlines,struct,zarr

;          Update z array and fitted line wavelengths, and add to structure
           gmos_updatez,struct.param,waverest,zarr

           goodcomp = gmos_checkcomp(struct,sigcut=checkcomp_sigcut)
           ngood = n_elements(goodcomp)

;;            if ngood lt ncomp OR goodcomp[0] eq -1 then begin
;;               if goodcomp[0] ne -1 then begin
;;                  ncomp = ngood
;;                  zarr = zguess
;;                  mcomp = zarr - zarr[0]
;;                  checkcomp_sigcut = 1d
;;               endif else begin
;;                  ncomp = 0
;;                  zarr = zguess
;;                  mcomp = -1d
;;               endelse
;;               print,'Repeating the fit with ',ncomp,' component(s).',format='(A,I0,A)'
;;               goto,fit
;;            endif

           if ngood lt ncomp then begin
              ncomp--
              if ncomp eq 1 then checkcomp_sigcut = 1d
              print,'Repeating the fit with ',ncomp,' component(s).',format='(A,I0,A)'
              goto,fit
           endif

        endif else begin

           struct = structinit

        endelse

;       Save result to binary file
        struct = add_tag(struct,zarr,'z')
        struct = add_tag(struct,fitran,'fitrange')
        savegen,struct=struct,file=outfile

     endfor

  endfor

  print,'Total time for calculation: ',systime(1)-starttime,' s.',$
        format='(/,A0,I0,A0,/)'

end
