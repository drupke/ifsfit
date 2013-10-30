;
; History
;  13apr30  DSNR  copies from gmos_fit_flat_continuum
;
function fit_no_continuum,lambda,flux,weight,template_flux,index,$
                          ct_coeff,no_dust=no_dust,quiet=quiet,$
                          addnorm=addnorm

  continuum = dblarr(n_elements(lambda))
  if keyword_set(addnorm) then continuum += addnorm
  ct_coeff=0d
  return,continuum

end
