; docformat = 'rst'
;
;+
;
; Convolve v2.3 BPASS templates to resolution of KCWI data with BM grating for 
; input to PPXF.
; 
; Resolution appears to be 1 A. Spacing of data points certainly is (see
; Byrne & Stanway 2023). Chisholm et al. 2019 discuss "resolution" of BPASS
; v2.1 models.
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
;      2024mar20, DSNR, created
;
; :Copyright:
;    Copyright (C) 2024 David S. N. Rupke
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
function ifsf_convolve_kcwibm_bpass,templam,temp,datlam,dat,newdat,z

   sigtofwhm = 2d*sqrt(2d*alog(2d))

   ; assuming here that Christy convolved things for constant R in log-lambda space
   ; multiply here b/c R is bigger for smaller dlambda or dvel
   ; This is equal to 10000
   sig_mod_A = 1.0d*(1d + z) / sigtofwhm
   sig_dat_A = 2.5d / sigtofwhm

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
      newtemp = gauss_smooth(temp,sig_datminmod_pix,/edge_truncate) $
   else $
      for i=0,ntemp-1 do $
         newtemp[*,i] = gauss_smooth(temp[*,i],sig_datminmod_pix,/edge_truncate)

   return,newtemp

end
