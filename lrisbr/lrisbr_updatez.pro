;
; History
;   10mar18     DSNR  created
;
; Update redshift array to reflect most recent emission-line fit.
;
pro lrisbr_updatez,param,z,bluered

  ppoff = param[0]
  ncomp = param[1]
  ppoff0 = ppoff - (ncomp-1)

  if bluered eq 'blue' then begin
     z.blue[0] = param[2]
     z.bluetilt = param[3]
     if ncomp gt 1 then $
        z.blue[1:ncomp-1] = param[ppoff0:ppoff-1]
  endif else begin
     z.gas[0] = param[2]
     z.redtilt = param[3]
     if ncomp gt 1 then $
        z.gas[1:ncomp-1] = param[ppoff0:ppoff-1]
  endelse

end
