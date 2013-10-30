;
; History
;   2009     DSNR  created
;   10jan27  DSNR  major revisions
;
; Update redshift array to reflect most recent emission-line fit.
;
pro lris_updatez,param,z

  ppoff = param[0]
  ncomp = param[1]
  ppoff0 = ppoff - (ncomp-1)

; Collect redshifts (in stellar restframe)
  z.blue[0] = param[2]
  z.gas[0]  = param[3]
  if ncomp gt 1 then begin
     z.gas[1:ncomp-1] = [param[ppoff0:ppoff-1]]
     z.blue[1:ncomp-1] = [param[ppoff0:ppoff-1] + param[2] - param[3]]
  endif
  z.bluetilt = param[4]

end
