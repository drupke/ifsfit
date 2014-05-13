; docformat = 'rst'
;
;+
; This function is input to MPFIT, which uses it to compute the
; NaD + HeI spectrum.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    An NaD + HeI spectrum.
;
; :Params:
;    wave: in, required, type=dblarr(N)
;      Wavelengths
;    param: in, required, type=dblarr
;      Best-fit parameter array output by MPFIT.
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
;      2014may09, DSNR, created
;    
; :Copyright:
;    Copyright (C) 2014 David S. N. Rupke
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
function ifsf_nadfcn, wave, param

   c = 299792.458d
;  NaD wavelength ratio (red to blue)
   lratio = 1.001014158d

;  Numbers of components
   nhei = param[0]
   nnadabs = param[1]
   nnadem = param[2]

   nwave = n_elements(wave)

;  Emission
   modflux = dblarr(nwave)+1d
   for i=0,nhei-1 do begin
      arg = ((wave-param[3+i*3])/$
            (2d*param[3+i*3]*param[4+i*3]/c))^2d
      modflux += param[5+i*3]*exp(-arg)
   endfor
   ilo = 3+nhei*3+nnadabs*4
   for i=0,nnadem-1 do begin
      arg1 = ((wave-param[ilo+i*4])/$
             (2d*param[ilo+i*4]*param[1+ilo+i*4]/c))^2d
      arg2 = ((wave-lratio*param[ilo+i*4])/$
             (2d*lratio*param[ilo+i*4]*param[1+ilo+i*4]/c))^2d
      modflux += param[2+ilo+i*4]*(exp(-arg1) + param[3+ilo+i*4]*exp(-arg2))
   endfor
;  Absorption
   modabs = dblarr(nwave)+1d
   ilo = 3+nhei*3
   for i=0,nnadabs-1 do modabs *= ifsf_cmpnad(wave,param[ilo+i*4:ilo+(i+1)*4-1])
   modflux *= modabs
   
   return,modflux

end
