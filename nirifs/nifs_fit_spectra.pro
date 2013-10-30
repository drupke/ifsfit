pro nifs_fit_spectra,gal,bin,time=time,verbose=verbose,$
                     cols=cols,rows=rows,cpass=cpass
;
; History
;   13mar05  DSNR  created
;
  
  starttime = systime(1)
  if ~ keyword_set(time) then time=0
  if keyword_set(verbose) then quiet=0 else quiet=1
  if ~ keyword_set(cpass) then cpass=0

  fcnz = 'nifs_redshift_spec'
  fcncheckcomp = 'nifs_checkcomp'

; Get fit initialization
  initdat = nifs_initfit_spectra(gal,bin,cpass=cpass)
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
  fitran_rest_lo = initdat.fitran_rest_lo
  fitran_rest_hi = initdat.fitran_rest_hi
  vdisp = initdat.vdisp
  argscheckcomp = initdat.argscheckcomp
  sigguess = initdat.sigguess
  argslinelist = initdat.argslinelist
  doubleline = initdat.doubleline
  fitblr = initdat.fitblr

  if gal eq 'mrk231' then begin
     data = readfits(infile,header,ext=1,/silent)
     var = readfits(infile,ext=2,/silent)
     ;; dq = readfits(infile,ext=3,/silent)
     if bin eq 999 then begin
        datasize = size(data)
        ncols = 1
        nrows = 1
        nz    = datasize[1]
     endif else begin
        datasize = size(data)
        ncols = datasize[1]
        nrows = datasize[2]
        nz    = datasize[3]
     endelse
     wave = dindgen(nz)
     crval = double(sxpar(header,'CRVAL3',/silent))
     crpix = double(sxpar(header,'CRPIX3',/silent))
     cdelt = double(sxpar(header,'CD3_3',/silent))
  endif else if gal eq 'f08572nw' then begin
     data = mrdfits(infile,0,header,/silent)
     err = mrdfits(infile,1,/silent)
     dq = mrdfits(infile,2,/silent)
     datasize = size(data)
     nz = datasize[1]
     ncols = datasize[2]
     nrows = datasize[3]
     wave = dindgen(nz)
;    Convert wavelength solution to A; OSIRIS default is nm
     crval = double(sxpar(header,'CRVAL1',/silent)) * 10d
     crpix = double(sxpar(header,'CRPIX1',/silent))
     cdelt = double(sxpar(header,'CDELT1',/silent)) * 10d
  endif

  wave = crval + cdelt*(wave-crpix+1) 

  if ~ keyword_set(cols) then cols=[1,ncols] $
  else if n_elements(cols) eq 1 then cols = [cols,cols]
  cols = fix(cols)
  for i=cols[0]-1,cols[1]-1 do begin

     print,'Column ',i+1,' of ',ncols,format='(A,I0,A,I0)'

     if ~ keyword_set(rows) then rows=[1,nrows] $
     else if n_elements(rows) eq 1 then rows = [rows,rows]
     rows = fix(rows)
     for j=rows[0]-1,rows[1]-1 do begin

        print,'  Row ',j+1,' of ',nrows,format='(A,I0,A,I0)'

  ; Get linelist
        if keyword_set(argslinelist) then $
           linelist = call_function('nifs_initlinelist',$
                                    _extra=argslinelist[i,j]) $
        else linelist = nifs_initlinelist()


        if gal eq 'mrk231' then begin
           if bin eq 999 then begin
              flux = data
              err_use = sqrt(abs(var))
           endif else begin
              flux = data[i,j,*]
              err_use = sqrt(abs(var[i,j,*]))
           endelse
;          Remove Brd line introduced by telluric standard
           flux_norm = gaussian(wave,[0.14,19451,16])
           flux_norm += 1d
           flux /= flux_norm
;          Remove a couple of areas near H_2 S(3)
           ibd = where((wave gt 20320 AND wave lt 20340) OR $
                       (wave gt 20275 AND wave lt 20295))
           ;; flux[ibd] = 0
           err_use[ibd] = max(err_use)*100d
        endif else if gal eq 'f08572nw' then begin
           flux = data[*,i,j]
           err_use = abs(err[*,i,j])
           ;; bad = dq[*,i,j]
;          Remove a bad region from specific spaxels near H_2 S(2)
           if ((i ge 22 AND i le 25) AND j eq 3) then begin
              flux[364:371] = 0
              err_use[364:371] = max(err_use)*100d
           endif
           if ((i eq 25 AND j eq 9) OR $
               (i eq 26 AND j eq 11)) then begin
              flux[386:397] = 0
              err_use[386:397] = max(err_use)*100d
           endif
           ;; if (i eq 24 and j eq 5) then begin
           ;;    flux[727:733] = 0
           ;;    err_use[727:733] = max(err_use)*100d
           ;; endif
           ;; if (i eq 16 and (j eq 3 OR j eq 4)) then begin
           ;;    flux[1103:1107] = 0
           ;;    err_use[1103:1107] = max(err_use)*100d
           ;; endif
;          Remove region near H_2 S(4)
           ih2s4_lo = value_locate(wave,2.00d)
           ih2s4_hi = value_locate(wave,2.01d)
           err_use[ih2s4_lo] = max(err_use)*100d
           err_use[ih2s4_hi] = max(err_use)*100d
;          Remove REALLY bad regions introduced by mosaicing
           ibad = where(flux gt 5 OR flux lt -1,ctbad)
           if ctbad gt 0 then begin
              for k=0,10 do begin
                 flux[ibad-k] = 0
                 flux[ibad+k] = 0
                 err[ibad-k] = max(err_use)*100d
                 err[ibad+k] = max(err_use)*100d
              endfor
           endif
        endif

;       Apply DQ plane
        ;; indx_bad = where(bad gt 0,ct)
        ;; if ct gt 0 then begin
        ;;    flux[indx_bad] = 0d
        ;;    err_use[indx_bad] = max(err_use)*100d
        ;; endif

        nodata = where(flux ne 0d AND ~ finite(flux,/nan),ct)

        if ct gt 200 then begin

;          Set any NaNed values so that flux is 0 and error is high
           inan = where(finite(flux,/nan),ctnan)
           if ctnan gt 0 then begin
              flux[inan] = 0d
              err_use[inan] = max(err_use)*100d
           endif

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

;       Add 'doubleline' parameter to argsinitpar
           if keyword_set(argsinitpar) then begin
              argsinitpar_use = argsinitpar
              argsinitpar_use = $
                 jjadd_tag(argsinitpar_use,$
                           'doubleline',doubleline[i,j])
           endif else argsinitpar_use = {doubleline:doubleline[i,j]}

;       Redshift fit range
           fitran_rest = [fitran_rest_lo[i,j],fitran_rest_hi[i,j]]
           fitran = nifs_redshift_spec(fitran_rest,z)
;       Add proper redshift structure to continuum fitting routine
           if keyword_set(argscontfit) then begin
              argscontfit_use = argscontfit
              ;; if tag_exist(argscontfit_use,'z') then begin
              ;;    argscontfit_use.z.star = z.star
              ;;    argscontfit_use.z.gas = z.gas ;[0]
              ;; endif
              argscontfit_use = $
                 jjadd_tag(argscontfit_use,$
                           'fitblr',fitblr[i,j])
              argscontfit_use = $
                 jjadd_tag(argscontfit_use,$
                           'z',z)
           endif
           
;          Reset sigguess
           sigguess = initdat.sigguess

           structinit = fit_spectrum(wave,flux,err_use,startempfile,z,ncomp,$
                                     /subtract,obj_id=gal,$
                                     linelist=linelist,quiet=quiet,$
                                     time=time,fitran=fitran,$
                                     fcninitpar=fcninitpar,$
                                     argsinitpar=argsinitpar_use,$
                                     disperse=disperse,$
                                     fcncontfit=fcncontfit,$
                                     argscontfit=argscontfit_use,$
                                     fcnz=fcnz,vdisp=vdisp,$
                                     fcnlinefit=fcnlinefit,$
                                     sigguess=sigguess)
           
           if ~ quiet then print,'FIT STATUS: ',structinit.fitstatus
           if structinit.fitstatus eq -16 then goto,nofit
           
           if ncomp gt 0 AND ~ keyword_set(onefit) then begin
              
              if ncomp gt 2 then $
                 nifs_orderlines,structinit,$
                                 doubleline=doubleline[i,j]
              nifs_updatez,structinit.param,z
              z.star = z.gas[0] + dzstel
              linepars = sepfitpars(structinit.param,structinit.perror)
              masklines = reform(linepars.wave[*,0:ncomp-1],$
                                 n_elements(linepars.wave[*,0:ncomp-1]))
              masksig=2d
              maskwidths = masksig*reform(linepars.sigma[*,0:ncomp-1],$
                                          n_elements(linepars.sigma[*,0:ncomp-1]))
              fitran = nifs_redshift_spec(fitran_rest,z)
              if keyword_set(argscontfit) then begin
                 if tag_exist(argscontfit_use,'z') then begin
                    argscontfit_use.z.star = z.star
                    argscontfit_use.z.gas = z.gas ;[0]
                 endif
              endif
              
              struct = fit_spectrum(wave,flux,err_use,startempfile,z,ncomp,$
                                    /subtract,obj_id=gal,$
                                    linelist=linelist,quiet=quiet,$
                                    time=time,fitran=fitran,$
                                    fcninitpar=fcninitpar,$
                                    argsinitpar=argsinitpar_use,$
                                    maskwidths=maskwidths,$
                                    masklines=masklines,$
                                    sigguess=linepars.sigma,$
                                    peakguess=linepars.fluxpk,$
                                    disperse=disperse,$
                                    fcncontfit=fcncontfit,$
                                    argscontfit=argscontfit_use,$
                                    fcnz=fcnz,vdisp=vdisp,$
                                    fcnlinefit=fcnlinefit)
           
              if ~ quiet then print,'FIT STATUS: ',struct.fitstatus
              if struct.fitstatus eq -16 then goto,nofit

              nifs_updatez,struct.param,z

              if ~ forcecomp then begin
                 goodcomp = $
                    call_function(fcncheckcomp,struct,z,$
                                  doubleline=doubleline[i,j],$
                                  _extra=argscheckcomp)
                 ngood = n_elements(goodcomp)
                 if ngood lt ncomp OR goodcomp[0] eq -1 then begin
                    ncomp--
                    ;; if ncomp eq 1 then checkcomp_sigcut = 3d
                    print,'    Repeating the fit with ',ncomp,$
                          ' component(s).',format='(A,I0,A)'
                    goto,fit
                 endif
              endif

           endif else begin
         
              if ncomp gt 2 then begin
                 nifs_orderlines,structinit,$
                                 doubleline=doubleline[i,j]
                 nifs_updatez,structinit.param,z
              endif
              struct = structinit
      
           endelse

           struct = jjadd_tag(struct,'z',z)
           struct = jjadd_tag(struct,'fitrange',fitran,/array_tag)

;; ;          Return error to original state
;;            struct.spec_err = err_orig[struct.gd_indx]

           save,struct,file=string(outdir,gal,'_',i+1,'_',j+1,'.xdr',$
                                   format='(A,A,A,I04,A,I04,A)')

nofit:
        
        endif else begin

           print,'  NIFS_FIT_SPECTRA: No data.'

        endelse

     endfor

  endfor

  print,'Total time for calculation: ',systime(1)-starttime,' s.',$
        format='(/,A0,I0,A0,/)'

end
