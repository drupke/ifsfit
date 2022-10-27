; docformat = 'rst'
;
;+
;
; Convolve CKC14 templates from Christy (which she convolved to sigma = 30 km/s)
; to resolution of KCWI data with BM grating for input to PPXF.
; Output sigma from PPXF should then be approximately correct.
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
;      2021jun25, DSNR, created
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
function ifsf_convolve_gmos_ckc,templam,temp,datlam,dat,newdat,z

   sigtofwhm = 2d*sqrt(2d*alog(2d))

   ; assuming here that Christy convolved things for constant R in log-lambda space
   ; multiply here b/c R is bigger for smaller dlambda or dvel
   sig_mod_R = 299792d/30d ; 30 km/s is sigma already! http://www.lco.cl/wp-content/uploads/2021/02/MAGEhandout2021.pdf
   sig_dat_A = 2.5d / sigtofwhm
   ; model sigma will be higher by (1 + z) because templam is redshifted
   sig_mod_A = templam/sig_mod_R

   ;  Convolve templates with difference in sigmas, since
   ;  templates have higher resolution over all wavelengths
   size_temp = size(temp)
   ntemplam = n_elements(templam)
   if size_temp[0] eq 1 then $
      ntemp = 1 $
   else $
      ntemp = size_temp[2] ; assume more than one template

   ; copies of old templates, data
   newtemp = temp
   newdat = dat

   sig_datminmod_A = sqrt(sig_dat_A^2d - sig_mod_A^2d)
   sig_datminmod_pix = sig_datminmod_A / (templam[1]-templam[0]) ; assume template dispersion is constant
   if ntemp eq 1 then $
      newtemp = ifsf_filter_gauss1d(temp,sig_datminmod_pix) $
   else $
      for i=0,ntemp-1 do $
         newtemp[*,i] = ifsf_filter_gauss1d(temp[*,i],sig_datminmod_pix)

   return,newtemp

end
