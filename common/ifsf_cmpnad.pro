; docformat = 'rst'
;
;+
;
; Compute Na D profile. (Presently, absorption only.)
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Flat spectrum normalized to unity, Na D absorption included.
;
; :Params:
;    wave: in, required, type=dblarr(N)
;      Wavelength array over which to compute profile.
;    pars: in, required, type=dblarr(4)
;      Parameters of profile: covering factor, optical depth, central
;      wavelength, and sigma (in km/s). Optical depth and wavelength
;      apply to D2 line (blue).
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
;      2010jul22, DSNR, created
;      2013nov21, DSNR, documented, renamed, added license and copyright
;      2014may13, DSNR, now uses sigma parameter instead of b (Doppler param).
;      2014may14, DSNR, fixed floating underflow
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
function ifsf_cmpnad,wave,pars
  
   c = 299792.458d
;  wavelength ratio (red to blue)
   lratio = 1.001014158d
;  optical depth ratio (blue to red)
   tratio = 2.0093d
  
   denom = pars[2]*pars[3]/c
   arg1 = (wave-pars[2])/denom
   arg2 = (lratio * wave - pars[2])/denom
   arg1 = arg1^2d
   arg2 = arg2^2d
;  The "mask" parameter eliminates floating underflow by removing very large
;  negative exponents. See http://www.idlcoyote.com/math_tips/underflow.html
;  for more details.
   mask1 = (arg1 lt 80)
   mask2 = (arg2 lt 80)
   arg12 = mask1*exp(-arg1/2d*mask1) + mask2*tratio*exp(-arg2/2d*mask2)
   exparg12 = exp(-pars[1]*arg12)
  
   yabs = 1d - pars[0]*(1d -exparg12)

   return,yabs

end
