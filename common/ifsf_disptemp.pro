; docformat = 'rst'
;
;+
;
; Convolve stellar templates with a Gaussian in wavelength space. To
; work properly, 
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Convolved templates.
;
; :Params:
;    template: in, required, type=dblarr(N\,M)
;      Stellar continuum templates for continuum fitting, with N
;      pixels and M templates. The template must have constant
;      dispersion in wavelength space (default) or velocity space (if
;      LOGLAM set).
;    lambda: in, required, type=dblarr(N)
;      Wavelengths
;    sigma: in, required, type=double
;      Sigma, in velocity velocity, of Gaussian to use in convolution.
;
; :Keywords:
;    loglam: in, optional, type=byte, default=0
;      Set to do convolution in velocity, or log(lambda), space. 
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
;      2009nov24  DSNR  added log(lam) treatment
;      2013oct10, DSNR, added documentation
;      2013nov13, DSNR, renamed, added license and copyright 
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
function ifsf_disptemp, template, lambda, sigma, loglam=loglam

  npoints = size(template)
  new_temp = dblarr(npoints[1], npoints[2])

  if keyword_set(loglam) then begin

     uselam = alog(lambda)
     disp = uselam[1] - uselam[0]
     sigma_pix = (sigma / 299792d) / disp

  endif else begin
     
     uselam = lambda
     disp = uselam[1] - uselam[0]
     sigma_pix = (sigma / 299792d) * mean(lambda) / disp

  endelse

  kernelsize = 10d*sigma_pix
  
; Properly center the kernel!  Otherwise one gets a wavelength shift 
; in the resulting convolved array.
;  kernel = dindgen(kernelsize) - kernelsize/2d
  kernel = [dindgen(kernelsize)-fix(kernelsize),$
            dindgen(kernelsize +1)]
  kernel = exp( -kernel^2d/(2d*sigma_pix^2d) )
  kernel /= total(kernel)

  for i = 0, npoints[2] - 1 do $
     new_temp[0,i] = convol(template[*,i], kernel, /edge_wrap)

  return, new_temp

end
