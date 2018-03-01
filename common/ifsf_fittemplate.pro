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
;      2018feb16, DSNR, created
;    
; :Copyright:
;    Copyright (C) 2018 David S. N. Rupke
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
function ifsf_fittemplate,lambda,flux,weight,template_wave,template_flux,index,$
                          ct_coeff,z,quiet=quiet,$
                          templatefcn=templatefcn,fcnargs=fcnargs,$
                          parinfo=parinfo,npar=npar
                         
   if ~ keyword_set(quiet) then quiet=0b
   if ~ keyword_set(parinfo) then parinfo=0b

   err=1d/sqrt(weight)
   ilambda=lambda[index]
   iflux=flux[index]
   iweight=weight[index]
   ierr=err[index]

   fcnargsnew = {twave:template_wave[index],tflux:template_flux[index]}
   if keyword_set(fcnargs) then fcnargscat = [fcnargs,fcnargsnew] $
   else fcnargscat = fcnargsnew

   if keyword_set(parinfo) then begin
      if tag_exist(parinfo,'value') then param=parinfo.value $
      else param = dblarr(npar)
   endif else param = dblarr(npar)
   yfit = mpcurvefit(ilambda,iflux,iweight,param,function_name=templatefcn,$
                     functargs=fcnargscat,parinfo=parinfo,$
                     quiet=quiet,status=status,errmsg=errmsg,$
                     /noderiv)
   
   if status eq 0 OR status eq -16 then message,'MPFIT: '+errmsg
   if status eq 5 then message,'MPFIT: Max. iterations reached.',/cont

   fcnargsnew = {twave:template_wave,tflux:template_flux}
   if keyword_set(fcnargs) then fcnargscat = [fcnargs,fcnargsnew] $
   else fcnargscat = fcnargsnew
   call_procedure,templatefcn,wave,param,continuum,_extra=fcnargscat
               
   ct_coeff = param
   return,continuum

end
