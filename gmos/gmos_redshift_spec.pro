;
; History
;  10may27  DSNR  created
;
; Redshift data using the specified redshift.

function gmos_redshift_spec,lambda,z,gas=gas,icomp=icomp

; This keyword selects other velocity components for the "base"
; redshift if desired.
  if keyword_set(icomp) then icomp = icomp else icomp=0
; Set base wavelength based on gas or stellar redshifts.  Default is
; stellar.
  if keyword_set(gas) then zbase = z.gas[icomp] else zbase = z.star

  newlambda = lambda * (1d + zbase)

  return,newlambda

end
