function fit_two_continua,lambda,flux,weight,template_flux,index,$
                          no_dust=no_dust,quiet=quiet,$
                          stitchwave=stitchwave,rfitoff=rfitoff
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

  if ~ keyword_set(stitchwave) then stitchwave = 5600d
  if ~ keyword_set(rfitoff) then rfitoff = 300d

;  Indices for continuum fitting
  index_blue = cmset_op(index,'AND',where(lambda lt stitchwave-100d))
  index_red = cmset_op(index,'AND',where(lambda gt stitchwave+rfitoff))
;  Indices locating all blue / red data
  index_allblue = where(lambda lt stitchwave)
  index_allred = where(lambda ge stitchwave)

;  Continuum fit blue part of spectrum only
   backfit_blue = ibackfit(flux[index_blue], lambda[index_blue], $
                           template_flux[index_blue, *],$
                           invvar=weight[index_blue], $
                           nodust=no_dust, quiet=quiet)
   backfit_red = ibackfit(flux[index_red], lambda[index_red], $
                          template_flux[index_red, *],$
                          invvar=weight[index_red], $
                          nodust=no_dust, quiet=quiet)

   continuum_blue = backfit_blue.starcoeff##template_flux[index_allblue,*]
   continuum_red = backfit_red.starcoeff##template_flux[index_allred,*]
   continuum = [continuum_blue,continuum_red]

   return,continuum

end
