;
; History
;  10jan28  DSNR  created
;
; Redshift a rest-frame spectrum using the specified redshift.

function lris_redshift_spec,lambda,z,stitchwave=stitchwave,$
                            gas=gas,icomp=icomp

; This keyword selects other velocity components for the "base"
; redshift if desired.
  if keyword_set(icomp) then icomp = icomp else icomp=0
; Set base wavelength based on gas or stellar redshifts.  Default is
; stellar.
  if keyword_set(gas) then zbase = z.gas[icomp] else zbase = z.star
; Redshift everything using 1 redshift ...
  if ~ keyword_set(stitchwave) then newlambda = lambda * (1d + zbase) $
; ... or include red/blue differences if desired.
  else begin

     iblue = where(lambda lt stitchwave * (1d + zbase))
     ired = where(lambda ge stitchwave * (1d + zbase))
     
     newlambda = lambda
     newlambda[iblue] = lambda[iblue] * (1d + zbase + (z.blue[0] - z.gas[0]) + $
                                         z.bluetilt*(lambda[iblue] - z.blueint))
     newlambda[ired] = lambda[ired] * (1d + zbase)

  endelse

  return,newlambda

end
