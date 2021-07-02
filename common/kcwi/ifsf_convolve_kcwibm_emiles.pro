; docformat = 'rst'
;
;+
;
; Convolve KCWI data with BM grating for fitting with E-MILES models,
; for input to PPXF.
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
function ifsf_convolve_kcwibm_emiles,templam,temp,datlam,dat,newdat,z

   sigtofwhm = 2d*sqrt(2d*alog(2d))

   ; model resolution higher in A because templates are redshifted
   sig_mod1_A = 3.00d*(1d + z) / sigtofwhm
   sig_mod2_A = 5.00d*(1d + z) / sigtofwhm
   sig_dat_A = 2.50d / sigtofwhm
   ; locations in observed frame where model resolution changes
   mergelam1_A = 3061d*(1d + z)
   mergelam2_A = 3541d*(1d + z)

   size_temp = size(temp)
   ntemp = size_temp[2] ; assume more than one template

   ; copies of old templates, data
   newtemp = temp
   newdat = dat

;  Convolve data with difference in sigmas in two blue regions since data has 
;  higher resolution here. In third region, resolutions of data and models
;  are identical.
   mergeinddat1 = value_locate(datlam,mergelam1_A)
   sig_modmindat1_A = sqrt(sig_mod1_A^2d - sig_dat_A^2d)
   sig_modmindat1_pix = sig_modmindat1_A / (datlam[1]-datlam[0]) ; assume data dispersion is constant

   mergeinddat2 = value_locate(datlam,mergelam2_A)
   sig_modmindat2_A = sqrt(sig_mod2_A^2d - sig_dat_A^2d)
   sig_modmindat2_pix = sig_modmindat2_A / (datlam[1]-datlam[0]) ; assume data dispersion is constant

   ; Have to edge_truncate (replicates last point at edges) or result is set to 0
   ; at edges
   if mergeinddat1 ne -1 then $
      newdat[0:mergeinddat1] = $
         gauss_smooth(dat[0:mergeinddat1],sig_modmindat1_pix,/edge_truncate)
   if mergeinddat2 ne -1 then $
      newdat[mergeinddat1+1:mergeinddat2] = $
         gauss_smooth(dat[mergeinddat1+1:mergeinddat2],sig_modmindat2_pix,$
         /edge_truncate)

   return,newtemp

end
