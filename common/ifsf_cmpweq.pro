; docformat = 'rst'
;
;+
;
; Compute equivalent widths for the specified emission lines. 
; Uses models of emission lines and continuum, and integrates over both using 
; the "rectangle rule."
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
;      2019mar20, DSNR, rewritten to match treatment of flux in IFSF_SEPFITPARS
;    
; :Copyright:
;    Copyright (C) 2013--2019 David S. N. Rupke
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
function ifsf_cmpweq,instr,linelist,doublets=doublets

  ncomp = instr.param[1]
  nlam = n_elements(instr.wave)
  lines = linelist->keys()
  nlines = n_elements(lines)

  tot = hash()
  comp = hash()
  dwave = instr.wave[1:nlam-1] - instr.wave[0:nlam-2]

  foreach line,lines do begin
     tot[line] = 0d
     comp[line] = dblarr(ncomp)
     for j=1,ncomp do begin
        modlines = ifsf_cmplin(instr,line,j,/velsig)
        if n_elements(modlines) ne 1 then $
           comp[line,j-1] = total(-modlines[1:nlam-1]/instr.cont_fit[1:nlam-1]*dwave) $
        else comp[line,j-1] = 0d
        tot[line] += comp[line,j-1]
     endfor     
  endforeach

; Special doublet cases: combine fluxes from each line
  if keyword_set(doublets) then begin
     sdoub = size(doublets)
     if sdoub[0] eq 1 then ndoublets = 1 else ndoublets = sdoub[2]
     for i=0,ndoublets-1 do begin
        if linelist.haskey(doublets[0,i]) AND linelist.haskey(doublets[1,i]) then begin
;          new line label
           dkey = doublets[0,i]+'+'+doublets[1,i]
;          add fluxes
           tot[dkey] = tot[doublets[0,i]]+tot[doublets[1,i]]
           comp[dkey] = comp[doublets[0,i]]+comp[doublets[1,i]]
        endif
     endfor
  endif

  return,{tot: tot,comp: comp}

end
