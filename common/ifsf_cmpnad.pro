; docformat = 'rst'
;
;+
;
; Compute Na D profile.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Flat spectrum Na D absorption included.
;
; :Params:
;    wave: in, required, type=dblarr(N)
;      Wavelength array over which to compute profile.
;    pars: in, required, type=dblarr(4)
;      Parameters of profile: covering factor, optical depth, central
;      wavelength, and Doppler parameter (in km/s)
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
;      10jul22  DSNR  created
;      2013nov21, DSNR, documented, renamed, added license and copyright 
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
; wavelength ratio (red to blue)
  lratio = 1.001014158d
; optical depth ratio (blue to red)
  tratio = 2.0093d
  
  denom = pars[2]*pars[3]/c
  arg1 = (wave-pars[2])/denom
  arg2 = (lratio * wave - pars[2])/denom
  arg1 = arg1^2d
  arg2 = arg2^2d
  arg12 = exp(-arg1) + tratio*exp(-arg2)
  exparg12 = exp(-pars[1]*arg12)
  
  yabs = 1d - pars[0]*(1d -exparg12)

  return,yabs

end
