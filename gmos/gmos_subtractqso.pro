;
; History
;  10oct07  DSNR  created
;
function gmos_subtractqso,instr,qsocntargs=qsocntargs,qsomod=qsomod
  
  wave = instr.wave

; polynomial near Ha/[NII]
  ypoly = dblarr(n_elements(wave))
  if n_elements(param) gt 1 then begin
     npoly = instr.param[4]
     wavelo = instr.param[5]
     wavehi = instr.param[6]
     n2hawave = where(wave ge wavelo AND wave le wavehi)
     ypoly[n2hawave] = poly(wave[n2hawave]-mean(wave[n2hawave]),instr.param[7:7+npoly])
  endif

  specstars = instr.spec - instr.specfit + ypoly
  qsotmpspec = readspec(qsocntargs.qsotmp)
  qsotmpwave = qsotmpspec[*,0]
  qsotmpflux = qsotmpspec[*,1]
  iqsotmpflux = interpol(qsotmpflux,qsotmpwave,wave)
  gmos_qso_cnt_fcn,wave,instr.ct_coeff,qsomod,$
                   fitord=qsocntargs.fitord,qsoflux=iqsotmpflux,$
                   qsoord=qsocntargs.qsoord,expterms=qsocntargs.expterms,$
                   /qsoonly
  qsomod += ypoly
  specresid = specstars - qsomod

  return,specresid

end
