; docformat = 'rst'
;
;+
;
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
;    template_wave: in, required, type=dblarr
;      Redshifted stellar template wavelength array.
;    template_flux: in, required, type=dblarr
;      Stellar templates fluxes.
;    index: in, required, type=intarr
;      Contains indices of continuum regions to fit
;    ct_coeff: out, required, type=integer, default=0
;      When stellar templates are fit, contains best-fit coefficients.
;    zstar: in, required, type=double
;      Stellar redshift.
;
; :Keywords:
;    blrpar: in, optional, type=dblarr(3N)
;      Initial guesses for broad Gaussian components. N is number of 
;      components. Parameters are flux, wavelength, and sigma in wavelength.
;    fitord: in, optional, type=integer, default=3
;      Specifies order of additive Legendre renormalization
;    index_log: in, optional, type=intarr
;      Same as index, but in log(lambda) space.
;    refit: in, optional, type=byte
;      Refit using stellar templates.
;    add_poly_degree: in, optional, type=int, def=30
;      Max. order of additive Legendre polynomial in PPXF.
;    quiet: in, optional, type=byte
;    qsoord: in, optional, type=integer, default=3
;      Specifies order of multiplicative Legendre renormalization
;    qsoonly: in, optional, type=byte
;      Do not add stellar continuum.
;    qsoxdr: in, required, type=string
;      XDR file containing IFSF fit to QSO continuum (i.e., the QSO template).
;    siginit_stars: in, optional, type=double, def=50
;      Initial sigma for stellar template fitting
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
;      2015jan21, DSNR, copied from IFSF_MULTICONT
;      2015sep04, DSNR, added option to include broad components as a way
;                       to either fit a BLR or to deal with scattered
;                       light issues with BLR emission (PG1411+442 et al.)
;      2016sep22, DSNR, now allow refit with stellar continuum; fixed bug 
;                       where QSOONLY input to fit
;      2016nov11, DSNR, removed EXPTERMS, added normalization
;    
; :Copyright:
;    Copyright (C) 2015--2016 David S. N. Rupke
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
function ifsf_fitqsohost,lambda,flux,weight,template_wave,template_flux,index,$
                         ct_coeff,zstar,quiet=quiet,$
                         blrpar=blrpar,qsoxdr=qsoxdr,qsoonly=qsoonly,$
                         index_log=index_log,refit=refit,$
                         add_poly_degree=add_poly_degree,$
                         siginit_stars=siginit_stars,$
                         polyspec_refit=polyspec_refit,fitran=fitran,$
                         fittol=fittol

   c = 299792.458d        ; speed of light, km/s
                         
   if ~ keyword_set(quiet) then quiet=0b
   if ~ keyword_set(qsoonly) then qsoonly=0b
   if ~ keyword_set(blrpar) then begin
      blrpar=0b
      blrterms=0
   endif else blrterms=n_elements(blrpar)
   if ~ keyword_set(add_poly_degree) then add_poly_degree = 30
   if ~ keyword_set(siginit_stars) then siginit_stars=50d

   restore,file=qsoxdr
   qsowave = qsotemplate.wave
   qsoflux_full = qsotemplate.flux

   iqsoflux = where(qsowave ge fitran[0]*0.99999d AND $
                    qsowave le fitran[1]*1.00001d)
   qsoflux = qsoflux_full[iqsoflux]

;  re-normalize template
   qsoflux /= median(qsoflux)

   err=1d/sqrt(weight)
   ilambda=lambda[index]
   iflux=flux[index]
   iweight=weight[index]
   ierr=err[index]

   itot = 16
   npar = itot+blrterms

   fcn = 'ifsf_qsohostfcn'
   fcnargs = {qsoflux: qsoflux[index],$
              qsoonly: qsoonly,$
              blrterms: blrterms,$
              nxfull: n_elements(lambda),$
              ixfull: index}
   parinfo = REPLICATE({limited:[0b,0b],value:0d,limits:[0d,0d],$
                        fixed:0b,mpprint:0b},npar)
   parinfo[0:15].limited=[1b,0b]
   param = dblarr(npar)
   if keyword_set(blrpar) then begin
      ngauss = fix(blrterms/3)
      param[itot:npar-1] = blrpar
      for i=0,ngauss-1 do begin
         parinfo[itot+3*i].limited = [1b,0b]
         parinfo[itot+3*i+1].fixed = 1b
;         parinfo[itot+3*i+1].limited = [1b,1b]
;         parinfo[itot+3*i+1].limits = [blrpar[3*i+1]-50d,blrpar[3*i+1]+50d]
         parinfo[itot+3*i+2].limited = [1b,1b]
         parinfo[itot+3*i+2].limits = [2000d,6000d]/299792d*blrpar[3*i+1]
      endfor
   endif
    

;;  Normalize data so it's near 1
;   fluxmed = median(iflux)
;   iflux /= fluxmed
;   iweight *= fluxmed^2d
   
   if ~ keyword_set(fittol) then fittol=1d-10
   fluxfit = mpcurvefit(ilambda,iflux,iweight,param,function_name=fcn,$
                        functargs=fcnargs,/noderiv,parinfo=parinfo,$
                        quiet=quiet,status=status,errmsg=errmsg,$
                        itmax=10000,tol=fittol)
   
   if status eq 0 OR status eq -16 then message,'MPFIT: '+errmsg
   if status eq 5 then message,'MPFIT: Max. iterations reached.',/cont

;;  Un-normalize
;   param[0:itot-1]*=fluxmed
;   iflux *= fluxmed
;   iweight /= fluxmed^2d
;   fluxfit *= fluxmed
;   if blrterms gt 0 then begin
;      iblrnorm = indgen(blrterms)*3
;      param[itot+iblrnorm]*=fluxmed
;   endif

   ct_coeff = param

   if keyword_set(refit) then begin                       

      ifsf_qsohostfcn,lambda,ct_coeff,qsomod,qsoflux=qsoflux,$
                      blrterms=blrterms,/qsoonly
      resid=flux-qsomod

;     Log rebin residual
      log_rebin,fitran,resid,resid_log,lambda_log,velscale=velscale
      log_rebin,fitran,err^2d,errsq_log
      err_log = sqrt(errsq_log)
        
;     Interpolate template to same grid as data
      temp_log = ifsf_interptemp(lambda_log,alog(template_wave),$
                                 template_flux)

;;     Normalize data so it's near 1
;      fluxmed = median(resid_log)
;      resid_log /= fluxmed
;      err_log /= fluxmed

      ppxf,temp_log,resid_log,err_log,velscale,$
           [0,siginit_stars],sol,$
           goodpixels=index_log,bestfit=continuum_log,moments=2,$
           degree=add_poly_degree,polyweights=polyweights,quiet=quiet,$
           weights=ct_coeff_refit

;;     un-normalize
;      continuum_log *= fluxmed
;      polyweights *= fluxmed
;      ct_coeff_refit *= fluxmed

      ct_coeff = {qso_host: param,$
                  stel: ct_coeff_refit,$
                  poly: polyweights,$
                  ppxf_sigma: sol[1]}

;     host can't be negative
      ineg = where(continuum_log lt 0,ctneg)
      if ctneg gt 0 then continuum_log[ineg]=0d

;     Resample the best fit into linear space
      cont_resid=interpol(continuum_log,lambda_log,ALOG(lambda))
      continuum = qsomod+cont_resid

;     Adjust stellar redshift based on fit
      zstar += sol[0]/c

;     Refit QSO residual???
;      newresid=flux-cont_resid
;      fcnargs.qsoonly=1b
;      newparam = dblarr(npar)
;      if keyword_set(blrpar) then newparam[npar-blrterms:npar-1] = blrpar
;      fluxfit = mpcurvefit(ilambda,newresid[index],iweight,newparam,$
;                           function_name=fcn,/quiet,$
;                           functargs=fcnargs,/noderiv,parinfo=parinfo)
;      ifsf_qsohostfcn,lambda,newparam,newqsomod,expterms=expterms,$
;                      fitord=fitord,qsoflux=qsoflux,qsoord=qsoord,$
;                      /qsoonly,blrterms=blrterms
;      continuum = newqsomod+cont_resid
;      ct_coeff = newparam
      
   endif else begin
      
      ifsf_qsohostfcn,lambda,param,continuum,qsoflux=qsoflux,$
                      qsoonly=qsoonly,blrterms=blrterms
   endelse

   return,continuum

end
