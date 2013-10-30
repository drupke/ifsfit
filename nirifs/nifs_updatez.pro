;
; History
;   13mar05  DSNR  created
;
; Update redshift array to reflect most recent emission-line fit.
;
pro nifs_updatez,param,z

  ppoff = param[0]
  ncomp = param[1]
  ppoff0 = ppoff - (ncomp-1)

; Collect redshifts (in stellar restframe)
  z.gas[0]  = param[2]
  if ncomp gt 1 then $
     z.gas[1:ncomp-1] = [param[ppoff0:ppoff-1]]
     
end
