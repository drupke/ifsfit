function gmos_fit_continuum,lambda,flux,weight,template_flux,index,$
                            ct_coeff,no_dust=no_dust,quiet=quiet
; REVISION HISTORY:
;     09aug14  DSNR  created
;     09dec11  DSNR  fixed bug in call to ibackfit: invar->invvar

  backfit = ibackfit(flux[index],lambda[index],template_flux[index, *], $
                     invvar=weight[index], nodust=no_dust, quiet=quiet)
  continuum = backfit.starcoeff##template_flux
;  Convert weights to errors
  err = 1d/sqrt(weight)
;  Fit line to red part of continuum and data for normalization
  fitord=3
  fluxfit = poly_fit(lambda[index],flux[index],$
                     fitord,measure=err[index])
  cntfit = poly_fit(lambda[index],continuum[index],$
                    fitord,measure=err[index])
;  Normalize the fitted continuum to the flux.
  continuum /= polycomp(lambda,cntfit)
  continuum *= polycomp(lambda,fluxfit)

  ct_coeff=backfit.starcoeff
  
  return,continuum

end
