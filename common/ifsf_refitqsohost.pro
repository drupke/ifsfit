; docformat = 'rst'
;
;+
;
; Re-fit host part of quasar+host fit with stellar templates.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    The best fit continuum spectrum (over all wavelengths, not just
;    those fit).
;
; :Params:
;    lambda: in, required, type=dblarr(N)
;    flux: in, required, type=dblarr(N)
;    weight: in, required, type=dblarr(N)
;    template_flux: in, required, type=any
;      Ignored. When stellar templates are fit, contatins templates.
;    index: in, required, type=intarr
;      Contains indices of continuum regions to fit
;    ct_coeff: out, required, type=integer, default=0
;      Best-fit coefficients.
;    ct_coeff_in: in, required, type=integer
;      Best-fit coefficients for IFSF_FITQSOHOST fit.
;
; :Keywords:
;    blrpar: in, optional, type=dblarr(3N)
;      Initial guesses for broad Gaussian components. N is number of 
;      components. Parameters are flux, wavelength, and sigma in wavelength.
;    expterms: in, optional, type=integer, default=0
;      Number of exponential terms by which to normalize, up to 4
;    fitord: in, optional, type=integer, default=3
;      Specifies order of additive renormalization
;    quiet: in, optional, type=byte
;    qsoord: in, optional, type=integer, default=3
;      Specifies order of multiplicative renormalization
;    qsoxdr: in, required, type=string
;      XDR file containing IFSF fit to QSO continuum (i.e., the QSO template).
; 
; :Author:
;    David S. N. Rupke::
;      Rhodes College
;      Department of Physics
;      2000 N. Parkway
;      Memphis, TN 38104
;      drupke@gmail.com
;
; :History:
;    Change History::
;      2016sep22, DSNR, copied from IFSF_FITQSOHOST
;    
; :Copyright:
;    Copyright (C) 2016 David S. N. Rupke
;
;    This program is free software: you can redistribute it and/or
;    modify it under the terms of the GNU General Public License as
;    published by the Free Software Foundation, either version 3 of
;    the License or any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;    General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program.  If not, see
;    http://www.gnu.org/licenses/.
;
;-
function ifsf_refitqsohost,lambda,flux,weight,template_wave,template_flux,index,$
                           ct_coeff,ct_coeff_in,$
                           linelistz,maskwidths,zstar,$
                           quiet=quiet,expterms=expterms,fitord=fitord,$
                           qsoxdr=qsoxdr,qsoord=qsoord,blrpar=blrpar,$
                           siginit_stars=siginit_stars

   c = 299792.458d        ; speed of light, km/s
                           
   if not keyword_set(quiet) then quiet=0b
   if ~ keyword_set(fitord) then fitord=3
   if ~ keyword_set(qsoord) then qsoord=3
   if ~ keyword_set(qsoonly) then qsoonly=0b
   if ~ keyword_set(expterms) then expterms=0
   if ~ keyword_set(blrpar) then begin
      blrpar=0b
      blrterms=0
   endif else blrterms=n_elements(blrpar)

   ilambda=lambda[index]
   iflux=flux[index]
   iweight=weight[index]
   err=1d/sqrt(weight)
   ierr = err[index]

   restore,file=qsoxdr
   qsowave = struct.wave
   qsoflux = struct.cont_fit

;  re-normalize template
   qsoflux /= median(qsoflux)

   ifsf_qsohostfcn,lambda,ct_coeff_in,qsomod,fitord=fitord,$
                   qsoord=qsoord,expterms=expterms,qsoflux=qsoflux,$
                   blrterms=blrterms,/qsoonly
   resid=flux-qsomod

;  Log rebin galaxy spectrum
   log_rebin,[lambda[0],lambda[n_elements(lambda)-1]],resid,$
             resid_log,lambda_log,velscale=velscale
   log_rebin,[lambda[0],lambda[n_elements(lambda)-1]],err^2d,$
             errsq_log
   err_log = sqrt(errsq_log)
        
;  Interpolate template to same grid as data
   temp = ifsf_interptemp(lambda,templatelambdaz,template.flux)
   temp_log = ifsf_interptemp(lambda_log,alog(templatelambdaz),$
                              template.flux)

;  Mask emission lines in log space
   ct_indx_log = $
      ifsf_masklin(exp(lambda_log), linelistz, maskwidths,$
                   nomaskran=nomaskran)

   polyterms = 30

   ppxf,temp_log,resid_log,err_log,velscale,$
        [0,siginit_stars],sol,$
        goodpixels=ct_indx_log,bestfit=continuum_log,moments=2,$
        degree=polyterms,polyweights=polyweights,quiet=quiet,$
        weights=ct_coeff

;  Resample the best fit into linear space
   cont_resid=interpol(continuum_log,lambda_log,ALOG(lambda))
   continuum = qsomod+cont_resid

;  Adjust stellar redshift based on fit
   zstar += sol[0]/c

   return,continuum

end
