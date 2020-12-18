; docformat = 'rst'
;
;+
;
; Copied from Cappellari's python gaussian_filter1d, in file log_rebin.py:
; 
; "Convolve a spectrum by a Gaussian with different sigma for every
;  pixel, given by the vector "sigma" with the same size as "spec".
;  If all sigma are the same this routine produces the same output as
;  scipy.ndimage.gaussian_filter1d, except for the border treatment.
;  Here the first/last p pixels are filled with zeros.
;  When creating  template library for SDSS data, this implementation
;  is 60x faster than the naive loop over pixels."
;
; In my case the first/last p pixels are filled with the original spectrum.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;
; :Params:
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
;      2020dec07, DSNR, created
;
; :Copyright:
;    Copyright (C) 2020 David S. N. Rupke
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
function ifsf_filter_gauss1d,spec,sig

   n = n_elements(sig)
   for i=0,n-1 do sig[i] = (sig[i] eq 0d) ? 0.01d : sig[i] ; forces sig=0 to be 0.01
   p = fix(ceil(max(3d*sig)))
   m = 2*p + 1 ; kernel size
   x2 = rebin((indgen(m) - p)^2,m,n)
   
   a = dblarr(m,n)
   for i=0,m-1 do a[i,p:n-p-1] = spec[i:n-m+i] ; loop over the small size of the kernel
   
   gau = exp(-x2/(2d*rebin(reform(sig,1,n),m,n)^2))
   gau /= rebin(reform(total(gau,1),1,n),m,n) ; normalize kernel
  
   conv_spectrum = total(a*gau,1)
   
   conv_spectrum[0:p-1] = spec[0:p-1]
   conv_spectrum[n-p:n-1] = spec[n-p:n-1]
   
   return,conv_spectrum

end
