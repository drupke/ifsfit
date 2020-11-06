;
;+
;
; This function gets the line syntax right for emission lines with brackets or
; other characters that need special treatment in file names.
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
;      2019apr03, DSNR, moved from IFSF_MAKEMAPS to standalone routine
;      2020jul02, DSNR, added treatment of asterisks
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
function ifsf_linesyntax,line
   ;  Change first, second occurrence of [ and ] to \[ and \]
   ;  Otherwise call to DEVICE can choke if the string includes [ rather
   ;  than \[
;   COMPILE_OPT IDL2, HIDDEN
   linelab=line
   ilb = strpos(linelab,'[')
   ilb2 = strpos(linelab,'[',/reverse_search) ; second occurrence
   if ilb ne -1 then $
      linelab = strmid(linelab,0,ilb)+'\'+strmid(linelab,ilb)
   if ilb ne -1 AND ilb ne ilb2 then $
      linelab = strmid(linelab,0,ilb2+1)+'\'+strmid(linelab,ilb2+1)
   irb = strpos(linelab,']')
   irb2 = strpos(linelab,']',/reverse_search) ; second occurrence
   if irb ne -1 then $
      linelab = strmid(linelab,0,irb)+'\'+strmid(linelab,irb)
   if irb ne -1 AND irb ne irb2 then $
      linelab = strmid(linelab,0,irb2+1)+'\'+strmid(linelab,irb2+1)
    iast = strpos(linelab,'*')
    if iast ne -1 then $
      linelab = strmid(linelab,0,iast)+'\'+strmid(linelab,iast)
   return,linelab
end
