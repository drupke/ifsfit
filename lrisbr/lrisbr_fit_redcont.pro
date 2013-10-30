function lrisbr_fit_redcont,lambda,flux,weight,template_flux,index,$
                            ct_coeff,no_dust=no_dust,quiet=quiet,$
                            order=order,ctBcoeff=ctBcoeff

;+
; History
;  10mar19  DSNR  created
;-
  
  err = 1d/sqrt(weight)

  continuum  = ctBcoeff##template_flux
; Fit line to red part of continuum and data for normalization
  if keyword_set(order) then fitord=order else fitord=5
  ctfitord=5
  fluxfit = poly_fit(lambda[index],flux[index],fitord,measure=err[index])
  fluxfit = mpfitfun('polycomp',lambda[index],flux[index],err[index],$
                     fluxfit,/quiet)
  cntfit = poly_fit(lambda,continuum,ctfitord,measure=err)
  cntfit = mpfitfun('polycomp',lambda,continuum,err,cntfit,/quiet)
; Normalize the fitted blue continuum to the red flux.
  continuum /= polycomp(lambda,cntfit)
  continuum *= polycomp(lambda,fluxfit)

  ct_coeff = ctBcoeff
  return,continuum
   
end
