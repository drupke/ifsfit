; docformat = 'rst'
;
;+
;
; Convolve stellar hizea template, which currently match HIRES data, to match 
; ESI data. Output sigma from PPXF should then be approximately correct.
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
;      2021dec14, DSNR, created
;
; :Copyright:
;    Copyright (C) 2021 David S. N. Rupke
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
pro ifsf_conv_resamp_hizea_esi,infile,cols,outxdr,fitran


   ;  Read input
   fxbopen,tablun,infile,1
   fxbreadm,tablun,cols,name,wave,cont

   ;  normalize continuum
   cont /= max(cont)
   ;  range to keep
   iran = where(wave ge fitran[0] AND wave le fitran[1])
   waveuse = wave[iran]
   contuse = cont[iran]

   sigtofwhm = 2d*sqrt(2d*alog(2d))
   ; multiply here b/c R is bigger for smaller dlambda or dvel
   sig_hires_R = 37500d * sigtofwhm ; 1.148 x 7 arcsec slit, the C5 decker.
   sig_esi_R = 4000d * sigtofwhm ; https://www2.keck.hawaii.edu/inst/esi/echmode.html
   sig_hires_kms = 299792d / sig_hires_R
   sig_esi_kms = 299792d / sig_esi_R
   sig_diff_kms = sqrt(sig_esi_kms^2d - sig_hires_kms^2d)
   sig_diff_pix = sig_diff_kms / 10d ; dispersion is dv = c * d(ln lam) = 10 km/s
   sig_diff_pix_arr = dblarr(n_elements(waveuse))+sig_diff_pix

   ; convolve with Gaussian that's constant in d(ln lam), or dv
   conv_contuse = ifsf_filter_gauss1d(contuse,sig_diff_pix_arr)

   template = {lambda: waveuse,$
               flux: [[conv_contuse]],$
               ages: [0d]}
   save,template,file=outxdr

end
