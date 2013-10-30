function gmos_fit_continuum,lambda,flux,weight,template_flux,index,$
                            no_dust=no_dust,quiet=quiet
;+
; NAME:
;     FIT_CONTINUUM
;
; PURPOSE:
;     Fit stellar continuum to spectrum.
;
; EXPLANATION:
;
; CALLING SEQUENCE
;     lris_fit_continuum,lambda,flux,weight,template_flux,index,\no_dust,
;                        \quiet
;
; INPUTS:
;     lambda - wavelength array
;     flux - flux array
;     weight - inverse variance array
;     template_flux - flux array for stellar templates
;     index - array of indices containing continuum regions to fit
;     no_dust - select this to turn off fitting dust extinction to stellar cont.
;     quiet - select this to suppress verbose fitting output
;
; OUTPUT:
;     The stellar continuum model
;
; METHOD:
;
; REVISION HISTORY:
;     09aug14  DSNR  created
;     09dec11  DSNR  fixed bug in call to ibackfit: invar->invvar
;-

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

   return,continuum

end
