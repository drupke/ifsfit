pro gmos_fit_spectra,gal,bin,time=time,verbose=verbose,$
                     rows=rows,cols=cols,fibers=fibers,$
                     onefit=onefit,keepnad=keepnad
;
; History
;  09jul08  DSNR  copied from LRIS routine to GMOS
;  10may27  DSNR  started re-write for new data
;
  
  starttime = systime(1)
  if ~ keyword_set(time) then time=0
  if keyword_set(verbose) then quiet=0 else quiet=1
  if keyword_set(fibers) then fibers=1 else fibers=0

  fcnz = 'gmos_redshift_spec'

; Get fit initialization
  initdat = gmos_initfit_spectra(gal,bin=bin)
  zinit = initdat.zinit
  ncompinit = initdat.ncomp
  infile = initdat.infile
  outdir = initdat.outdir
  fcninitpar = initdat.fcninitpar
  argsinitpar = initdat.argsinitpar
  fcnlinefit = initdat.fcnlinefit
  fcncontfit = initdat.fcncontfit
  argscontfit = initdat.argscontfit
  startempfile = initdat.startempfile
  fitran_rest = initdat.fitran_rest
  vdisp = initdat.vdisp
  checkcomp_sigcut_init=initdat.fluxsigthresh
  siglim = initdat.siglim
  sigguess = initdat.sigguess
  argslinelist = initdat.argslinelist
  
; Get linelist
  if keyword_set(argslinelist) then $
     linelist = call_function('gmos_initlinelist',_extra=argslinelist) $
  else linelist = gmos_initlinelist()

  if fibers then begin
;    Read data
     data = readfits(infile,header,ext=2,/silent)
     var = readfits(infile,ext=3,/silent)
     dq = readfits(infile,ext=4,/silent)
     datasize = size(data)
     nz    = datasize[1]
     ncols = datasize[2]
     nrows = 1
;    Wavelength solution
     wave = dindgen(nz)
     crval = double(sxpar(header,'CRVAL1',/silent))
     crpix = double(sxpar(header,'CRPIX1',/silent))
     cdelt = double(sxpar(header,'CDELT1',/silent))
     wave = crval + cdelt*(wave-crpix+1) 
  endif else begin
     data = readfits(infile,header,ext=1,/silent)
     var = readfits(infile,ext=2,/silent)
     dq = readfits(infile,ext=3,/silent)
     datasize = size(data)
     ncols = datasize[1]
     nrows = datasize[2]
     nz    = datasize[3]
     wave = dindgen(nz)
     crval = double(sxpar(header,'CRVAL3',/silent))
     crpix = double(sxpar(header,'CRPIX3',/silent))
;     cdelt = double(sxpar(header,'CDELT3',/silent))
     cdelt = double(sxpar(header,'CD3_3',/silent))
     wave = crval + cdelt*(wave-crpix+1) 
  endelse

; Add continuum sigma tag to argsinitpar
  if keyword_set(argsinitpar) then $
     argsinitpar = jjadd_tag(argsinitpar,'cont_sig',0d) $
  else $
     argsinitpar = {cont_sig:0d}

  if ~ keyword_set(cols) then cols=[1,ncols] $
  else if n_elements(cols) eq 1 then cols = [cols,cols]
  cols = fix(cols)
  for i=cols[0]-1,cols[1]-1 do begin

     print,'Column ',i+1,' of ',ncols,format='(A,I0,A,I0)'

     if ~ keyword_set(rows) then rows=[1,nrows] $
     else if n_elements(rows) eq 1 then rows = [rows,rows]
     rows = fix(rows)
     for j=rows[0]-1,rows[1]-1 do begin

        if fibers then begin
           flux = data[*,i]
;          absolute value takes care of a few deviant points
           err = sqrt(abs(var[*,i]))
           bad = dq[*,i]
        endif else begin
           print,'  Row ',j+1,' of ',nrows,format='(A,I0,A,I0)'
           flux = data[i,j,*]
           err = sqrt(abs(var[i,j,*]))
           bad = dq[i,j,*]
        endelse

;       Apply DQ plane
        indx_bad = where(bad gt 0,ct)
        if ct gt 0 then begin
           flux[indx_bad] = 0d
           err[indx_bad] = max(err)*100d
        endif

        checkcomp_sigcut = checkcomp_sigcut_init

        nodata = where(flux ne 0d,ct)
        if ct ne 0 then begin

           ncomp = abs(ncompinit[i,j])
           if ncompinit[i,j] lt 0 then forcecomp=1 else forcecomp=0
           disperse = 0
           dzstel = 0d

fit:
           
;       Initialize redshift structure
           if ncomp gt 0 then $
              z = {star:zinit[i,j,0]+dzstel,$
                   gas:dblarr(ncomp)+zinit[i,j,0:ncomp-1]} $
           else $
              z = {star:zinit[i,j,0]+dzstel,$
                   gas:-1}

;       Remove NaI D line for purposes of continuum fit by maximizing
;       error
           if ~keyword_set(keepnad) then begin
              nadran_rest = [5850d,5900d]
              nadran = gmos_redshift_spec(nadran_rest,z)
              indx_nad = where(wave ge nadran[0] AND wave le nadran[1],ct)
              if ct gt 0 then err[indx_nad]=max(err)
           endif

;       Remove strongest part of telluric feature
           ; telran = [6865d,6885d]
           ; indx_tel = where(wave ge telran[0] AND wave le telran[1],ct)
           ; if ct gt 0 then err[indx_tel]=max(err)

;       Redshift fit range
           fitran = gmos_redshift_spec(fitran_rest,z)
;       Add proper redshift structure to continuum fitting routine
           if keyword_set(argscontfit) then begin
              if tag_exist(argscontfit,'z') then begin
                 argscontfit.z.star = z.star
                 argscontfit.z.gas = z.gas[0]
              endif
           endif

;          Reset sigguess
           sigguess = initdat.sigguess

           structinit = fit_spectrum(wave,flux,err,startempfile,z,ncomp,$
                                     /subtract,obj_id=gal,$
                                     linelist=linelist,quiet=quiet,$
                                     time=time,fitran=fitran,$
                                     fcninitpar=fcninitpar,$
                                     argsinitpar=argsinitpar,$
                                     disperse=disperse,$
                                     fcncontfit=fcncontfit,$
                                     argscontfit=argscontfit,$
                                     fcnz=fcnz,vdisp=vdisp,$
                                     fcnlinefit=fcnlinefit,$
                                     sigguess=sigguess)
           
           if ~ quiet then print,'FIT STATUS: ',structinit.fitstatus
           if structinit.fitstatus eq -16 then goto,nofit
           
           if ncomp gt 0 AND ~ keyword_set(onefit) then begin
              
              if ncomp gt 1 then gmos_orderlines,structinit
              gmos_updatez,structinit.param,z
              z.star = z.gas[0] + dzstel
              linepars = sepfitpars(structinit.param,structinit.perror)
              masklines = reform(linepars.wave[*,0:ncomp-1],$
                                 n_elements(linepars.wave[*,0:ncomp-1]))
              masksig=2d
              maskwidths = masksig*reform(linepars.sigma[*,0:ncomp-1],$
                                          n_elements(linepars.sigma[*,0:ncomp-1]))
              fitran = gmos_redshift_spec(fitran_rest,z)
              if keyword_set(argscontfit) then begin
                 if tag_exist(argscontfit,'z') then begin
                    argscontfit.z.star = z.star
                    argscontfit.z.gas = z.gas[0]
                 endif
              endif
              
              struct = fit_spectrum(wave,flux,err,startempfile,z,ncomp,$
                                    /subtract,obj_id=gal,$
                                    linelist=linelist,quiet=quiet,$
                                    time=time,fitran=fitran,$
                                    fcninitpar=fcninitpar,$
                                    argsinitpar=argsinitpar,$
                                    maskwidths=maskwidths,$
                                    masklines=masklines,$
                                    sigguess=linepars.sigma,$
                                    peakguess=linepars.fluxpk,$
                                    disperse=disperse,$
                                    fcncontfit=fcncontfit,$
                                    argscontfit=argscontfit,$
                                    fcnz=fcnz,vdisp=vdisp,$
                                    fcnlinefit=fcnlinefit)
           
              if ~ quiet then print,'FIT STATUS: ',struct.fitstatus
              if struct.fitstatus eq -16 then goto,nofit

              gmos_updatez,struct.param,z

              if ~ forcecomp then begin
                 goodcomp = $
                    gmos_checkcomp(struct,z,$
                                   sigcut=checkcomp_sigcut,$
                                   siglim=siglim)
                 ngood = n_elements(goodcomp)
                 if ngood lt ncomp then begin
                    ncomp--
                    if ncomp eq 1 then checkcomp_sigcut = 1d
                    print,'    Repeating the fit with ',ncomp,$
                          ' component(s).',format='(A,I0,A)'
                    goto,fit
                 endif
              endif

           endif else begin
         
              if ncomp gt 1 then begin
                 gmos_orderlines,structinit
                 gmos_updatez,structinit.param,z
              endif
              struct = structinit
      
           endelse

           struct = jjadd_tag(struct,'z',z)
           struct = jjadd_tag(struct,'fitrange',fitran,/array_tag)

;; ;          Return error to original state
;;            struct.spec_err = err_orig[struct.gd_indx]

           if fibers then $
              save,struct,file=string(outdir,gal,'_',i+1,'.xdr',$
                                      format='(A,A,A,I04,A,I04,A)') $
           else $
              save,struct,file=string(outdir,gal,'_',i+1,'_',j+1,'.xdr',$
                                      format='(A,A,A,I04,A,I04,A)')

nofit:
        
        endif else begin

           print,'  GMOS_FIT_SPECTRA: No data.'

        endelse

     endfor

  endfor

  print,'Total time for calculation: ',systime(1)-starttime,' s.',$
        format='(/,A0,I0,A0,/)'

end
