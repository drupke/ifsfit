; docformat = 'rst'
;
;+
;
; Parse the output parameter file from the Na D fit.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    A structure containing the fit parameters.
;
; :Params:
;    parfile: in, required, type=string
;      Path and file name for parameter file.
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
;      2010jul21, DSNR, created
;      2013nov22, DSNR, documented, renamed, added copyright and
;                       license, moved output to single structure
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
function ifsf_readnadpar,parfile

  openr,lun,parfile,/get_lun
  line = ''
  readf,lun,line
  readf,lun,line
  linespl = strsplit(line,' ',/extract)
  nabs = fix(linespl[0])
  nem = fix(linespl[1])
  llo = double(linespl[2])
  lhi = double(linespl[3])
  readf,lun,line
  readf,lun,line
  readf,lun,line
  if nabs gt 0 then begin
     abspars = dblarr(4,nabs)
     for i=0,nabs-1 do begin
        for j=0,3 do begin
           readf,lun,line
           linespl = strsplit(line,' ',/extract)
           abspars[j,i] = double(linespl[4])
        endfor
     endfor
  endif
  if nem gt 0 then begin
     empars = dblarr(3,nem)
     for i=0,nem-1 do begin
        for j=0,2 do begin
           readf,lun,line
           linespl = strsplit(line,' ',/extract)
           empars[j,i] = double(linespl[4])
        endfor
     endfor
  endif
  free_lun,lun

  pars = {abs: abspars,$
          em: empars,$
          nabs: nabs,$
          nem: nem,$
          llo: llo,$
          lhi: lhi $
         }

  return, pars

end
