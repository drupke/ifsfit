; docformat = 'rst'
;
;+
;
; Compute emission line profile of a given line and velocity
; component.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Array of fluxes representing emission line profile.
;
; :Params:
;    instr: in, required, type=structure
;      Contains parameters of line profile. Required tags are PARAM,
;      which is an array of best fit line parameters output from MPFIT;
;      LINELABEL, which is an array of line labels; and WAVE, which is
;      the wavelength array of the full spectrum.
;    line: in, required, type=string
;      Name of line for which to compute profile.
;    comp: in, required, type=integer
;      Number of velocity component for which to compute line profile.
;
; :Keywords:
;    velsig: in, optional, type=byte, default=0
;      Set if line sigma in PARAM array is in velocity space (km/s).
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
;      2013sep12  DSNR  made into stand-alone routine
;      2013oct09, DSNR, added documentation
;      2013nov13, DSNR, renamed, added license and copyright 
;      2014apr10, DSNR, check if sigma = 0 to prevent IEEE exception in GAUSSIAN
;      2014apr14, DSNR, fixed bug due to change in ordering of lines in PARAM
;                       array; switched from ordered hashes to ordinary hashes
;                       a while ago, which unsynced ordering of lines from place
;                       to place
;      2014apr15, DSNR, turned off /DOUBLE in call to Gaussian b/c of floating
;                       point underflow; culprit seems to be MACHAR() ...
;    
; :Copyright:
;    Copyright (C) 2013-2014 David S. N. Rupke
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
function ifsf_cmplin,instr,line,comp,velsig=velsig

  c = 299792.458d
  
  iline = where(instr.linelabel eq line,ct)
  ppoff = instr.param[0]
  ncomp = instr.param[1]

  indices = where(instr.parinfo.line eq line AND instr.parinfo.comp eq comp)
  gausspar = instr.param[indices]
  if keyword_set(velsig) then gausspar[2] *= gausspar[1]/c
  if gausspar[2] eq 0d then flux = 0d else $
     flux = double(gaussian(instr.wave,gausspar))

  return,flux

end
