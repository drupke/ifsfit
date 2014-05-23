; docformat = 'rst'
;
;+
;
; Find best-fit Gaussian sigma with which to convolve stellar templates.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    bestsig: out, type=double
;      Best fit sigma.
;    bestcont: out, type=dblarr
;      Best fit continuum.
;    besttemp: out, typte=dblarr(N\,M)
;      Best fit convolved templates.
;
; :Params:
;    lambda: in, required, type=dblarr(N)
;      Wavelengths.
;    flux: in, required, type=dblarr(N)
;      Flux of spectrum to fit.
;    weight: in, required, type=dblarr(N)
;      Weights of spectrum to fit.
;    ct_indx: in, required, type=dblarr
;      Indices of input spectrum over which to perform fit.
;    ct_coeff: in, required, type=dblarr(M)
;      Coefficients of templates for best fit continuum.
;    template: in, required, type=dblarr(N\,M)
;      Undispersed templates, with N pixels and M templates
;    sigfitvals: in, required, type=dblarr
;      Sigma values to try for dispersion.
;    initstr: in, required, type=structure
;      Structure input to IFSF_FITSPEC.
;
; :Keywords:
;    quiet: in, optional, type=byte, default=0
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
;    ChangeHistory::
;      2013oct09, DSNR, created
;      2013nov21, DSNR, renamed, added license and copyright 
;    
; :Copyright:
;    Copyright (C) 2013 David S. N. Rupke
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
pro ifsf_optstelsig,lambda,flux,weight,ct_indx,ct_coeff,$
                    template,sigfitvals,initstr,$
                    bestsig,bestcont,besttemp,$
                    quiet=quiet

  if ~ keyword_set(quiet) then quiet=0b

  bestcc = 0                    ; highest correlation coefficient
  bestsig = 0                   ; sigma at minimum chi2
  bestcont = 0
  besttemp = 0

  for i=0,n_elements(sigfitvals)-1 do begin

;    Convolve template with Gaussian
     if initstr.sigfitvals[i] gt 0 then temp_tmp = $
        ifsf_disptemp(template,lambda,sigfitvals[i],$
                      loglam=loglam) $
     else temp_tmp = template
;    Add polynomials if requested
     if tag_exist(initstr,'argsaddpoly2temp') then temp_tmp = $
        call_function('ifsf_addpoly2temp',temp_tmp,$
                      _extra=initstr.argsaddpoly2temp) $
     else temp_tmp = $
        call_function('ifsf_addpoly2temp',temp_tmp)
;    Compute cont_tmp without re-fitting
     if tag_exist(initstr,'argscontfit') then cont_tmp = $
        call_function(initstr.fcncontfit,lambda,flux,weight,$
                      temp_tmp,ct_indx,ct_coeff,$
                      quiet=quiet,/ctinput,_extra=initstr.argscontfit) $
     else cont_tmp = $
        call_function(initstr.fcncontfit,lambda,flux,weight,$
                      temp_tmp,ct_indx,ct_coeff,$
                      quiet=quiet,/ctinpu)
     
;    Compute correlation coefficient
     cc = correlate(flux[ct_indx],cont_tmp[ct_indx],/double)
     if cc gt bestcc then begin
        bestcc = cc
        bestsig = sigfitvals[i]
        bestcont = cont_tmp
        besttemp = temp_tmp
     endif
           
  endfor

  if not quiet then $
     print,'IFSF_OPTSTELSIG: Using sigma = ',string(bestsig,'(I0)'),$
           ' km/s to disperse stellar template.'
  
end
