;
; Helper routine for IFSF_MAKEMAPS.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Postscript plots.
;
; :Params:
;    initproc: in, required, type=string
;      Name of procedure to initialize the fit.
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
function ifsf_plotrange,auto=auto,rline=rline,matline=matline,$
                        rcomp=rcomp,matcomp=matcomp,$
                        rquant=rquant,matquant=matquant,$
                        rncbdiv=rncbdiv,rlo=rlo,rhi=rhi,$
                        mapgd=mapgd,divinit=divinit,ncbdivmax=ncbdivmax
  ;
  if keyword_set(auto) then doauto=1 else doauto=0
  if ~ doauto then begin
    if keyword_set(rline) AND keyword_set(matline) AND $
      keyword_set(rquant) AND keyword_set(matquant) AND $
      keyword_set(rncbdiv) AND (keyword_set(rlo) OR keyword_set(rhi)) then begin
      if keyword_set(rcomp) AND keyword_set(matcomp) then $
        ithisline = where(rline eq matline AND $
        rcomp eq matcomp AND $
        rquant eq matquant,ctthisline) $
      else $
        ithisline = where(rline eq matline AND $
        rquant eq matquant,ctthisline)
      if ctthisline eq 1 then begin
        zran = [rlo[ithisline],rhi[ithisline]]
        ncbdiv = rncbdiv[ithisline]
        ncbdiv = ncbdiv[0]
      endif else doauto=1
    endif else doauto=1
  endif
  if doauto then begin
    if keyword_set(mapgd) AND $
      keyword_set(divinit) AND $
      keyword_set(ncbdivmax) then begin
      zran = [min(mapgd),max(mapgd)]
      divarr = ifsf_cbdiv(zran,divinit,ncbdivmax)
      ncbdiv = divarr[0]
    endif else begin
      message,'Proper keywords not specified.'
    endelse
  endif
  return,[zran,zran[1]-zran[0],ncbdiv]
end