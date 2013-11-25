; docformat = 'rst'
;
;+
;
; Compute total equivalent widths (over all components) for the
; specified emission lines. Uses models of emission lines and
; continuum, and integrates over both using the "rectangle rule."
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Array of equivalent widths.
;
; :Params:
;    instr: in, required, type=structure
;      Contains output of IFSF_FITSPEC.
;    lines: in, required, type=strarr
;      Names of line for which to compute equivalent widths.
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
;      2013nov25, DSNR, copied from PRINTWEQ
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
function ifsf_cmpweq,instr,lines

  ncomp = instr.param[1]
  nlam = n_elements(instr.wave)
  nlines = n_elements(lines)

  weq = dblarr(nlines)

  for i=0,nlines-1 do begin
     modlines = dblarr(nlam)
     for j=1,ncomp do modlines += uhsf_cmplin(instr,lines[i],j)
     dwave = instr.wave[1:nlam-1] - instr.wave[0:nlam-2]
     weq = total(-modlines[1:nlam-1]/instr.cont_fit[1:nlam-1]*dwave)
  endfor

  return,weq

end
