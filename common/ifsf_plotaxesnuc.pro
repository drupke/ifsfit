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
pro ifsf_plotaxesnuc,xran_kpc,yran_kpc,xnuc,ynuc,nolab=nolab,toplab=toplab,$
                     noxlab=noxlab,charsize=charsize,rightlab=rightlab,$
                     nonuc=nonuc,colornuc=colornuc,noylab=noylab,$
                     colorax=colorax

  COMPILE_OPT IDL2, HIDDEN
  if not keyword_set(colorax) then colornuc = !P.color
  if not keyword_set(charsize) then charsize=!P.charsize
  if not keyword_set(nolab) then begin
    if keyword_set(toplab) OR keyword_set(noxlab) then $
      cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save,xtickn=replicate(' ',60),$
      color=colorax
    if not keyword_set(noxlab) AND not keyword_set(toplab) then $
      cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save,charsize=charsize,$
      color=colorax
    if keyword_set(rightlab) then $
      cgaxis,yaxis=1,yran=yran_kpc,/ysty,/save,charsize=charsize,$
      color=colorax $
    else $
      cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save,charsize=charsize,$
      color=colorax
  endif else begin
    cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save,xtickn=replicate(' ',60),$
      color=colorax
    cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save,ytickn=replicate(' ',60),$
      color=colorax
  endelse
  if not keyword_set(toplab) OR keyword_set(noxlab) then $
    cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty,$
    color=colorax
  if keyword_set(toplab) AND $
    not keyword_set(noxlab) AND $
    not keyword_set(nolab) then $
    cgaxis,xaxis=1,xran=xran_kpc,/xsty,charsize=charsize,$
    color=colorax
  if not keyword_set(rightlab) then $
    cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty,$
    color=colorax $
  else $
    cgaxis,yaxis=0,yran=yran_kpc,ytickn=replicate(' ',60),/ysty,$
    color=colorax
  if not keyword_set(colornuc) then colornuc = !P.color
  if not keyword_set(nonuc) then cgoplot,xnuc,ynuc,psym=1,color=colornuc
end
;
