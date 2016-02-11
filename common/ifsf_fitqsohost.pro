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
;    template_flux: in, required, type=any
;      Ignored. When stellar templates are fit, contatins templates.
;    index: in, required, type=intarr
;      Contains indices of continuum regions to fit
;    ct_coeff: out, required, type=integer, default=0
;      When stellar templates are fit, contains best-fit coefficients.
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
;    qsoonly: in, optional, type=byte
;      Do not add stellar continuum.
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
;      2015jan21, DSNR, copied from IFSF_MULTICONT
;      2015sep04, DSNR, added option to include broad components as a way
;                       to either fit a BLR or to deal with scattered
;                       light issues with BLR emission (PG1411+442 et al.)
;    
; :Copyright:
;    Copyright (C) 2015 David S. N. Rupke
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
function ifsf_fitqsohost,lambda,flux,weight,template_flux,index,$
                         ct_coeff,quiet=quiet,expterms=expterms,fitord=fitord,$
                         qsoxdr=qsoxdr,qsoord=qsoord,qsoonly=qsoonly,$
                         blrpar=blrpar

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
   ierr = 1d/sqrt(iweight)

   restore,file=qsoxdr
   qsowave = struct.wave
   qsoflux = struct.cont_fit

;  re-normalize template
   qsoflux /= median(qsoflux)

   npar = (fitord+qsoord+2*expterms)*2+blrterms
   param = dblarr(npar)
   if keyword_set(blrpar) then param[npar-blrterms:npar-1] = blrpar

   fcn = 'ifsf_qsohostfcn'
   fcnargs = {fitord: fitord,$
              qsoflux: qsoflux[index],$
              qsoord: qsoord,$
              expterms: expterms,$
              blrterms: blrterms}
   parinfo = REPLICATE({limited:[1b,0b],limits:[0d,0d],$
                        fixed:0b,mpprint:0b},npar)

   fluxfit = mpcurvefit(ilambda,iflux,iweight,param,function_name=fcn,$
                        functargs=fcnargs,/noderiv,parinfo=parinfo,$
                        /quiet)
                       
   ifsf_qsohostfcn,lambda,param,continuum,expterms=expterms,$
                   fitord=fitord,qsoflux=qsoflux,qsoord=qsoord,$
                   qsoonly=qsoonly,blrterms=blrterms

   ct_coeff = param
   return,continuum

end
