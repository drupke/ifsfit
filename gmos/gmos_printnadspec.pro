;
; History
;  10jul20  DSNR  created
;  11may17  DSNR  changed Weq calculation so emission line ignored
;                 use median redshift to calculate fit range
;  13apr09  DSNR  changed doz = {gas: ...} statements to doz.gas = ...
;

pro gmos_printnadspec,instr,outfile,cfitord,waveout,fluxout,polypars,$
                      weq,qsocntargs=qsocntargs,restcomp=restcomp

  outran_rest    = [5810d,5960d]
; Mrk 231
;  cfit1ran_rest  = [5810d,5865d]
;  cfit2ran_rest  = [5900d,5960d]
  cfit1ran_rest  = [5810d,5865d]
  cfit2ran_rest  = [5905d,5960d]

;
; Redshift wavelength ranges
;
; Get # of components
  if n_elements(instr.param) gt 1 then ncomp = instr.param[1] $
  else ncomp = 0
  gas=1
  icomp=0
  doz=instr.z
; Pick median redshift
  if ncomp eq 0 then gas=0
  if ncomp eq 2 then doz.gas=mean(instr.z.gas)
  if ncomp eq 3 then begin
     zsort = sort(instr.z.gas)
     doz.gas=(instr.z.gas[zsort[1]])
  endif
  if keyword_set(restcomp) then doz.gas=instr.z.gas[restcomp-1]
; Do redshifting
  outran    = gmos_redshift_spec(outran_rest,doz,gas=gas,icomp=icomp)
  cfit1ran  = gmos_redshift_spec(cfit1ran_rest,doz,gas=gas,icomp=icomp)
  cfit2ran  = gmos_redshift_spec(cfit2ran_rest,doz,gas=gas,icomp=icomp)
  dcfit = cfit2ran[0] - cfit1ran[1]

; Compute continuum
  wave = instr.wave
  specstars = instr.spec - instr.specfit
  errstars = instr.spec_err

; For case where we fit a QSO template.
  if keyword_set(qsocntargs) then begin
     qsotmpspec = readspec(qsocntargs.qsotmp)
     qsotmpwave = qsotmpspec[*,0]
     qsotmpflux = qsotmpspec[*,1]
     iqsotmpflux = interpol(qsotmpflux,qsotmpwave,wave)
     gmos_qso_cnt_fcn,wave,instr.ct_coeff,modqso,$
                      fitord=qsocntargs.fitord,qsoflux=iqsotmpflux,$
                      qsoord=qsocntargs.qsoord,expterms=qsocntargs.expterms,$
                      /qsoonly
     specresid = specstars - modqso
     errresid = errstars
  endif else begin
     specresid = specstars
     errresid = errstars
  endelse

; Normalize continuum

  iwavefit = where((wave ge cfit1ran[0] AND wave le cfit1ran[1]) OR $
                   (wave ge cfit2ran[0] AND wave le cfit2ran[1]))
  wavefit = wave[iwavefit]
  fluxfit = specresid[iwavefit]
  errfit = errresid[iwavefit]

  parinfo = replicate({value:0d},cfitord)
  polypars = mpfitfun('poly',wavefit,fluxfit,errfit,$
                      parinfo=parinfo,/quiet)

  iwaveout = where(wave ge outran[0] AND wave le outran[1],ctout)
  waveout = wave[iwaveout]
  fluxout = specresid[iwaveout] / poly(waveout,polypars)
  errout = errresid[iwaveout] / poly(waveout,polypars)

; Compute equivalent width
  iweq2 = value_locate(waveout,cfit1ran[1])
  iweq3 = value_locate(waveout,cfit2ran[0])
  fluxoutnad = 1d -fluxout[iweq2:iweq3]
  erroutnad = errout[iweq2:iweq3]

; Filter out emission lines using 2-sigma cut
  iem = where(fluxoutnad lt -erroutnad*2d,ctem)
  if ctem gt 0 then fluxoutnad[iem] = 0d
  weq = total(fluxoutnad*(waveout[iweq2:iweq3]-waveout[iweq2-1:iweq3-1]))
  weq_e = sqrt(total(errout[iweq2:iweq3]^2*(waveout[iweq2:iweq3]-waveout[iweq2-1:iweq3-1])))
  weq = [weq,weq_e]

; Print spectrum to file

  openw,lun,outfile,/get_lun
  printf,lun,ctout
  for i=0,ctout-1 do $
     printf,lun,waveout[i],fluxout[i],errout[i]
  free_lun,lun

end
