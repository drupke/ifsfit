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
;    template_lambdaz: in, required, type=any
;      Ignored. When stellar templates are fit, contains redshifted
;      template wavelength.
;    template_flux: in, required, type=any
;      Ignored. When stellar templates are fit, contatins templates.
;    index: in, required, type=intarr
;      Contains indices of continuum regions to fit
;    ct_coeff: out, required, type=integer, default=0
;      When stellar templates are fit, contains best-fit coefficients.
;    zstar: in, required, type=any
;      Ignored. When stellar templates are fit, ...
;
; :Keywords:
;    quiet: in, optional, type=byte
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
;      2015jan21, DSNR, copied from IFSF_FITPOLY
;      2015may11, DSNR, fixed bug when only one region is used
;      2016oct21, DSNR, added new (but ignored) parameters to match 
;                       IFSF_FITSPEC calling rubric
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
function ifsf_fitmulticont,lambda,flux,weight,template_lambdaz,$
                           template_flux,index,ct_coeff,zstar,$
                           quiet=quiet,fitreg=fitreg,fitfcn=fitfcn,$
                           fitargs=fitargs

   if not keyword_set(quiet) then quiet=0b
   
   continuum = dblarr(n_elements(lambda))
   sizetmp = size(fitreg)
   if sizetmp[0] eq 1 then begin
      nfitreg=1
      fitreg_use = rebin(fitreg,2,1)
   endif else begin
      nfitreg=sizetmp[2]
      fitreg_use = fitreg
   endelse
   for i=0,nfitreg-1 do begin
      jlo = value_locate(lambda,fitreg_use[0,i])
      jhi = value_locate(lambda,fitreg_use[1,i])
      if jhi eq -1 OR jlo eq n_elements(lambda)-1 then begin
         print,'IFSF_FITMULTICONT: Error -- no data in fit region.'
         stop
      endif
      if jlo eq -1 then jlo = 0
      tmplambda = lambda[jlo:jhi]
      tmpflux = flux[jlo:jhi]
      tmpweight = weight[jlo:jhi]
      if template_lambdaz ne !NULL then $
         tmp_template_lambdaz = template_lambdaz[jlo:jhi]
      if template_flux ne !NULL then $
         tmp_template_flux = template_flux[jlo:jhi]
      tmpindex = index[where(index ge jlo AND index le jhi)]
      tmpindex -= jlo
      tmpcont = call_function(fitfcn[i],tmplambda,tmpflux,tmpweight,$
                              tmp_template_lambdaz,tmp_template_flux,tmpindex,$
                              ct_coeff,quiet=quiet,$
                              _extra=fitargs[string('reg',i+1,format='(A0,I0)')])
      continuum[jlo:jhi] = tmpcont
   endfor

   ct_coeff = 0
   return,continuum

end
