;
; History
;  10apr01  DSNR  created
;
; Redshift a rest-frame spectrum using the specified redshift.

function sp1_redshift_spec,lambda,z,gas=gas,icomp=icomp

; This keyword selects other velocity components for the "base"
; redshift if desired.
  if ~ keyword_set(icomp) then icomp=0
; Set base wavelength based on gas or stellar redshifts.  Default is
; stellar.
  if keyword_set(gas) then zbase = z.gas[icomp] else zbase = z.star
; Redshift everything using 1 redshift ...
  newlambda = lambda * (1d + zbase)

  return,newlambda

end
