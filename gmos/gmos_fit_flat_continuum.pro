function gmos_fit_flat_continuum,lambda,flux,weight,template_flux,index,$
                                 no_dust=no_dust,quiet=quiet,cfitran=cfitran
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
;-

;  GMOS
;  if ~keyword_set(cfitran) then cfitran=[6375,6475]
;  ifit = where(lambda ge cfitran[0] AND lambda le cfitran[1])
;  KPNO
  if ~keyword_set(cfitran) then cfitran=[6375,6475,6650,6700,6750,6850]
  ifit = where((lambda ge cfitran[0] AND lambda le cfitran[1]) OR $
               (lambda ge cfitran[2] AND lambda le cfitran[3]) OR $
               (lambda ge cfitran[4] AND lambda le cfitran[5]))

  err = 1d/sqrt(weight)
  fitord=3
  fluxfit = poly_fit(lambda[ifit],flux[ifit],$
                     fitord,measure=err[ifit])
  continuum = polycomp(lambda,fluxfit)
  
  return,continuum

end
