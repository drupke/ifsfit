; docformat = 'rst'
;
;+
;
; Compute the continuum spectrum fit by PPXF.
; 
; Unfortunately, this doesn't yet work, for reasons unclear. Though at least 
; one reason is the fact that I don't do the convolution with the best-fit 
; sigma.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Best-fit continuum spectrum, in linear space.
;
; :Params:
;    lambda: in, required, type=dblarr(nlam)
;      Wavelengths in linear space.
;    lambda_log: in, required, type=dblarr(nlam_log)
;      Wavelengths in ln(lambda) space.
;    temp: in, required, type=dblarr(nlam\,ntemp)
;      Continuum templates.
;    tempweights: in, required, type=dblarr(ntemp)
;      Best-fit weights for continuum templates.
;    polydeg: in, required, type=double
;      Degree of additive polynomials added to fit.
;    polyweights: in, required, type=dblarr(polydeg)
;      Weights of additive polynomials added to fit.
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
;    ChangeHistory::
;      2013dec09, DSNR, created
;      2013dec11, DSNR, testing and bug fix
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
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
;    General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program. If not, see
;    http://www.gnu.org/licenses/.
;
;-
function ifsf_cmpcontppxf,lambda,lambda_log,temp,tempweights,$
                          polydeg,polyweights
  
; Compute polynomial in log space
  x = cap_range(-1d,1d,n_elements(lambda_log))
  apoly_log = 0d
  for j=1,polydeg do apoly_log += legendre(x,j)*polyweights[j]
; Interpolate polynomial to linear space
  apoly = interpol(apoly_log,lambda_log,ALOG(lambda))
; Compute template combination in linear space
  spec = temp # tempweights

  spec += apoly

  return,spec

end
