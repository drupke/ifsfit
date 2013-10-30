function lris_fit_continuum,lambda,flux,weight,template_flux,index,$
                            ct_coeff,no_dust=no_dust,quiet=quiet,$
                            stitchwave=stitchwave,redord=redord
;+
; NAME:
;     LRIS_FIT_CONTINUUM
;
; PURPOSE:
;     Fit stellar continuum to LRIS data.
;
; EXPLANATION:
;
; CALLING SEQUENCE
;     lris_fit_continuum,lambda,flux,err,weight,template_flux,index,\no_dust,
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
;     stitchwave - wavelength at which LRIS blue / red modules were stitched
;
; OUTPUT:
;     The stellar continuum model
;
; METHOD:
;
; REVISION HISTORY:
;     09aug14  DSNR  created
;     09sep03  DSNR  changed red lower limit from 5800 to stitchwave+200
;     09dec11  DSNR  fixed bug in call to ibackfit: invar->invvar
;     10jan12  DSNR  changed red lower limit to stitchwave+100
;                    changed index_allred to index_red in fitting of
;                      red data
;-

   if ~ keyword_set(no_dust) then no_dust=0
   if ~ keyword_set(quiet) then quiet=0
   if ~ keyword_set(stitchwave) then stitchwave = 5600d

;  Convert weights to errors
   err = 1d/sqrt(weight)

;  Indices for continuum fitting
   index_blue = cmset_op(index,'AND',where(lambda lt stitchwave-100d))
   index_red = cmset_op(index,'AND',where(lambda gt stitchwave))
;  remove 5-sigma outliers for continuum fitting
;   sigflux = flux/err
;   index_red = cmset_op(index_red,'AND',where(sigflux gt 5d))
;  Indices locating all blue / red data
   index_allblue = where(lambda lt stitchwave)
   index_allred = where(lambda ge stitchwave)

;  Continuum fit blue part of spectrum only
   backfit = ibackfit(flux[index_blue], lambda[index_blue], $
                      template_flux[index_blue, *],$
                      invvar=weight[index_blue], $
                      nodust=no_dust, quiet=quiet)

;  Compute blue and red continua.
   continuum_blue = backfit.starcoeff##template_flux[index_allblue,*]
   continuum_red  = backfit.starcoeff##template_flux[index_allred,*]

;  Fit line to red part of continuum and data for normalization
   if keyword_set(redord) then fitord=redord else fitord=5
   ctfitord=5
   fluxfit = poly_fit(lambda[index_red],flux[index_red],$
                      fitord,measure=err[index_red])
   fluxfit = mpfitfun('polycomp',lambda[index_red],flux[index_red],$
                      err[index_red],fluxfit,/quiet)
   cntfit = poly_fit(lambda[index_allred],$
                     backfit.starcoeff##template_flux[index_allred,*],$
                     ctfitord,measure=err[index_allred])
   cntfit = mpfitfun('polycomp',lambda[index_allred],$
                     backfit.starcoeff##template_flux[index_allred,*],$
                     err[index_allred],cntfit,/quiet)

;  Test fitting
;   set_plot,'x'
;   loadct,0,/silent
;   plot,lambda[index_red],flux[index_red]
;   loadct,13,/silent
;   oplot,lambda[index_allred],polycomp(lambda[index_allred],fluxfit),color=255

;  For red continuum, normalize the fitted blue continuum to the red
;  flux.
   continuum_red /= polycomp(lambda[index_allred],cntfit)
   continuum_red *= polycomp(lambda[index_allred],fluxfit)

;; ;  Fit line to red part of continuum and data for normalization
;;    fluxfit = mpcurvefit(lambda[index_red],flux[index_red],$
;;                         measure=err[index_red])
;;    cntfit = linfit(lambda[index_red],$
;;                    backfit.starcoeff##template_flux[index_red,*],$
;;                    measure=err[index_red])
;; ;  For red continuum, normalize (with linear fit) the fitted blue
;; ;  continuum to the red flux.
;;    continuum_red /= cntfit[0] + cntfit[1]*lambda[index_allred]
;;    continuum_red *= fluxfit[0] + fluxfit[1]*lambda[index_allred]
      
   continuum = [continuum_blue,continuum_red]

   ct_coeff=backfit.starcoeff
   return,continuum
   
end
