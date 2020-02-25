;
; Helper routine for IFSF_MAKEMAPS.
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
;      2019dec13, DSNR, moved from IFSF_MAKEMAPS to standalone routine
;
; :Copyright:
;    Copyright (C) 2019 David S. N. Rupke
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
pro ifsf_plotcompass,xarr,yarr,carr=carr,nolab=nolab,hsize=hsize,hthick=hthick,$
                     thick=thick
  COMPILE_OPT IDL2, HIDDEN
  if ~ keyword_set(carr) then carr='Black'
  if ~ keyword_set(hsize) then hsize=!D.X_SIZE / 64
  if ~ keyword_set(hthick) then hthick=1d
  if ~ keyword_set(thick) then thick=!P.thick
  cgarrow,xarr[0],yarr[0],xarr[1],yarr[1],/data,/solid,color=carr,hsize=hsize,$
    hthick=hthick,thick=thick
  cgarrow,xarr[0],yarr[0],xarr[2],yarr[2],/data,/solid,color=carr,hsize=hsize,$
    hthick=hthick,thick=thick
  if ~ keyword_set(nolab) then begin
    cgtext,xarr[3],yarr[3],'N',color=carr,align=0.5
    cgtext,xarr[4],yarr[4],'E',color=carr,align=0.5
  endif
end

