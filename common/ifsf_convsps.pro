; docformat = 'rst'
;
;+
;
; Convolve stellar templates with a velocity dispersion. This assumes that the 
; convolution happens in log(lambda) space like in PPXF. LOSVD assumed to be
; Gaussian.
; 
; PPXF has a faster FFT algorithm, could implement for non-Gaussian LOSVD.
; 
;
; :Categories:
;    IFSFIT
;
; :Returns:
;     Convolved stellar model.
;     
; :Params:
;     lambda: in, required, type=dblarr(Nwave)
;        linear wavelength array, constant dispersion
;     temps: in, required, type=dblarr(Ntemp,Nwave)
;        template array
;     coeffs: in, required, type=dblarr(Ntemp)
;        best-fit linear coefficients (from, e.g., PPXF)
;     sigma: in, required, type=double
;        LOSVD velocity dispersion with which to convolve templates.
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
;      2020dec04, DSNR, created
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
function ifsf_convsps,lambda,temps,coeffs,sigma

;   ;  Data comes in as a linear wavelength grid, sigma as a velocity dispersion
;   smtemps = temps
;   ;  Because ppxf loops over templates, convolving each one separately, that's
;   ;  what I'm doing. Not clear that it's necessary, though it is slower ...
;   for i=0,n_elements(coeffs)-1 do begin
;      ; rebin to log space
;      log_rebin,[lambda[0],lambda[n_elements(lambda)-1]],$
;         temps[*,i],temp_log,lambda_log,velscale=velscale
;      ; compute sigma in pixels
;      sigma_pix = double(sigma)/double(velscale)
;      ; gaussian smooth; i.e., assuming LOSVD is Gaussian
;      smtemp_log = gauss_smooth(temp_log,sigma_pix)
;      ; interpolate back to linear
;      smtemps[*,i] = interpol(smtemp_log,lambda_log,ALOG(lambda))
;   endfor
;   ; multiply by coefficients and sum
;   smsps = smtemps # coeffs
       
   sps = temps # coeffs
   ; rebin to log space
   log_rebin,[lambda[0],lambda[n_elements(lambda)-1]],$
      sps,sps_log,lambda_log,velscale=velscale
   ; compute sigma in pixels
   sigma_pix = double(sigma)/double(velscale)
   ; gaussian smooth; i.e., assuming LOSVD is Gaussian
   smsps_log = gauss_smooth(sps_log,sigma_pix)
   ; interpolate back to linear
   smsps = interpol(smsps_log,lambda_log,ALOG(lambda))

   return,smsps

end
