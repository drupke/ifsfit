;
; History
;  10jan28  DSNR  created
;
; Redshift a rest-frame spectrum using the specified redshift.

function lrisbr_redshift_spec,lambda,z,bluered=bluered,$
                              gas=gas,icomp=icomp

; This keyword selects other velocity components for the "base"
; redshift if desired.
  if keyword_set(icomp) then icomp = icomp else icomp=0

  if keyword_set(bluered) then begin

;    BLUE
     if bluered eq 'blue' then begin
        if keyword_set(gas) then zbase = z.blue[icomp] else zbase = z.star
        newlambda = lambda * (1d + zbase + z.bluetilt*(lambda - z.blueint))
;        newlambda = lambda * (1d + zbase)

;    RED
     endif else begin
        if keyword_set(gas) then zbase = z.gas[icomp] else zbase = z.star
        newlambda = lambda * (1d + zbase + z.redtilt*(lambda - z.redint))
;        newlambda = lambda * (1d + zbase)
     endelse

  endif else begin
     print,'LRISBR_REDSHIFT_SPEC: ERROR: BLUERED keyword not set.'
  endelse

  return,newlambda

end
