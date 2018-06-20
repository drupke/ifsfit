; docformat = 'rst'
;
;+
;
; Compute singlet absorption profile for the red, low-tau line in a doublet.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Flat spectrum normalized to unity.
;
; :Params:
;    wave: in, required, type=dblarr(N)
;      Wavelength array over which to compute profile.
;    pars: in, required, type=dblarr(4)
;      Parameters of profile: covering factor, optical depth, central
;      wavelength, and sigma (in km/s). Optical depth and wavelength
;      apply to red line.
;    tratio: in, required, type=double
;      Optical depth ratio of blue line to that of red line.
;    lratio: in, required, type=double
;      Wavelength ratio of red line to that of blue line.
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
;      2018jun05, DSNR, copied from IFSF_CMPDOUBLET
;    
; :Copyright:
;    Copyright (C) 2018 David S. N. Rupke, Anthony To
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
function ifsf_cmpsinglet,wave,pars,tratio,lratio
  
   c = 299792.458d
   
   denom = pars[2]*pars[3]/c
   arg1 = (wave-pars[2])/denom
;   arg2 = (lratio * wave - pars[2])/denom
   arg1 = arg1^2d
;   arg2 = arg2^2d
;  The "mask" parameter eliminates floating underflow by removing very large
;  negative exponents. See http://www.idlcoyote.com/math_tips/underflow.html
;  for more details.
   mask1 = (arg1 lt 80)
;   mask2 = (arg2 lt 80)
;   arg12 = mask1*exp(-arg1/2d*mask1) + mask2*tratio*exp(-arg2/2d*mask2)
   arg12 = mask1*exp(-arg1/2d*mask1)
   exparg12 = exp(-pars[1]*arg12)
  
   yabs = 1d - pars[0]*(1d -exparg12)

   return,yabs

end
