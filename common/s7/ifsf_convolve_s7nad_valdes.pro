; docformat = 'rst'
;
;+
;
; Convolve Valdes templates and S7 data to common resolution for input to PPXF.
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
function ifsf_convolve_s7nad_valdes,templam,temp,datlam,dat,newdat,z

   sigtofwhm = 2d*sqrt(2d*alog(2d))

   sig_valdes_A = 1.35d*(1d + z) / sigtofwhm
   ; multiply here b/c R is bigger for smaller dlambda or dvel
   sig_s7b_R = 3000d * sigtofwhm
   sig_s7r_R = 7000d * sigtofwhm
   mergelam_A = 5600d ; where blue and red spectra are stitched

;  Convolve templates with difference in sigmas in blue, since
;  templates have higher resolution over most wavelengths
   mergeindtemp = value_locate(templam,mergelam_A)
   sig_s7b_A = templam[0:mergeindtemp]/sig_s7b_R
   
   sig_s7bminv_A = dblarr(mergeindtemp+1)
   for i=0,mergeindtemp do  $
      sig_s7bminv_A[i] = $
         sqrt(((sig_s7b_A[i]^2d - sig_valdes_A^2d) lt 0d) ? 0d : (sig_s7b_A[i]^2d - sig_valdes_A^2d))
   sig_s7bminv_pix = sig_s7bminv_A / (templam[1]-templam[0]) ; assume template dispersion is constant

   size_temp = size(temp)
   ntemp = size_temp[2] ; assume more than one template
   newtemp = temp
   for i=0,ntemp-1 do $
      newtemp[0:mergeindtemp,i] = ifsf_filter_gauss1d(temp[0:mergeindtemp,i],sig_s7bminv_pix)

;  Convolve data with difference in sigmas in red, since data has higher resolution
;  here
   mergeinddat = value_locate(datlam,mergelam_A)
   sig_s7r_A = datlam[mergeinddat+1:n_elements(datlam)-1]/sig_s7r_R
   sig_vmins7r_A = sqrt(sig_valdes_A^2d - sig_s7r_A^2d)
   sig_vmins7r_pix = sig_vmins7r_A / (datlam[1]-datlam[0]) ; assume data dispersion is constant

   newdat = dat
   newdat[mergeinddat+1:n_elements(dat)-1] = $
      ifsf_filter_gauss1d(dat[mergeinddat+1:n_elements(dat)-1],sig_vmins7r_pix)

;  Could convolve data further in blue where templates degrade further, but
;  this is only important for a few hundred A there; let's ignore

   return,newtemp

end
