; docformat = 'rst'
;
;+
;
; This procedure ...
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    IDL save file (.xdr)
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
;      2014jan24, DSNR, created
;    
; :Copyright:
;    Copyright (C) 2014 David S. N. Rupke
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
pro ifsf_makemaps,initproc

  fwhm2sig = 2d*sqrt(2d*alog(2d))

; Get fit initialization
  initdat=call_function(initproc)

; Get linelist
  linelist = ifsf_linelist(initdat.lines)

; Restore line maps
  if not tag_exist(initdat,'outlines') then outlines = linelist->keys() $
  else outlines = initdat.outlines
  restore,file=initdat.outdir+initdat.label+'.lin.xdr'

  size_tmp = size(linmaps[outlines[0])
  dx = size_tmp[1]
  dy = size_tmp[2]

  foreach line,outlines do begin

    npx = 3
    npy = initdat.maxncomp
    !P.multi=[0,npx,npy]
    ps_start,file=initdat.mapdir+initdat.label+line+'.eps',charsize=2,/encap,$
       /inches,xs=sizeunit*npx,ys=sizeunit*npy,/qui
    cbform = '(I0)'
    mmar=[0.25,0,2.25,3.2]
    
    for i=0,initdat.maxncomp-1 do begincbform = '(D0.1)'

;      Flux
       map = linemaps[line,*,*,i,0]
       ibd = where(map eq bad,ctbd)
       igd = where(map ne bad,ctgd)
       zran = [min(map),max(map[igd])]
       mapscl = bytscl(rebin(linmaps[line,*,*,i,0],dx*20,dy*20,/sample),$
                       min=zran[0],max=zran[1])
       dzran = zran[1]-h_zran[0]
       ncbdiv = 5
       cgimage,mapscl,/keep,multimar=mmar,oposition=p
       cgplot,[0],xsty=5,ysty=5,xran=xran,yran=yran,position=p,$
          /nodata,/noerase

    endfor

    !P.multi=0
    ps_end

  endforeach

end
