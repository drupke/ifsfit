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
  plotquantum = 2.5 ; in inches
  bad = 1d99

; Get fit initialization
  initdat=call_function(initproc)

; Get linelist
  linelist = ifsf_linelist(initdat.lines)

; Restore line maps
  if not tag_exist(initdat,'outlines') then outlines = linelist->keys() $
  else outlines = initdat.outlines
  restore,file=initdat.outdir+initdat.label+'.lin.xdr'

  size_tmp = size(linmaps[outlines[0]])
  dx = size_tmp[1]
  dy = size_tmp[2]
  center = [double(dx)/2d,double(dy)/2d]

; Luminosity and angular size distances
  ldist = lumdist(initdat.zsys,H0=73,Omega_m=0.27,Lambda0=0.73,/silent)
  kpc_per_as = ldist/(1+initdat.zsys^2)*1000d/206265d

; scales in kpc
  xran_kpc = double([-(center[0]-0.5),dx-(center[0]-0.5)]) $
             * platescale * kpc_per_as
  yran_kpc = double([-(center[1]-0.5),dy-(center[1]-0.5)]) $
             * platescale * kpc_per_as

; Size of plot grid
  npx = initdat.maxncomp
  npy = 3

; Linemap indices to plot below: flux, wavelength (converted to velocity), sigma
  ilinmap = [0,2,3]

  foreach line,outlines do begin

;   Get syntax of linelabel right; otherwise call to DEVICE chokes
    linelab=line
    ilb = strpos(linelab,'[')
    if ilb ne -1 then $
       linelab = strmid(linelab,0,ilb)+'\'+strmid(linelab,ilb)
    irb = strpos(linelab,']')
    if irb ne -1 then $
       linelab = strmid(linelab,0,irb)+'\'+strmid(linelab,irb)

    cgps_open,initdat.mapdir+initdat.label+linelab+'.eps',charsize=2,/encap,$
       /inches,xs=plotquantum*npx,ys=plotquantum*npy,/qui
    
    for i=0,initdat.maxncomp-1 do begin

       for j=0,2 do begin

          map = linmaps[line,*,*,i,ilinmap[j]]
          ibd = where(map eq bad,ctbd)
          igd = where(map ne bad,ctgd)
          if j eq 1 then map[igd] = map[igd]/...
          zran = [min(map),max(map[igd])]
          mapscl = bytscl(rebin(linmaps[line,*,*,i,ilinmap[j]],dx*20,dy*20,$
                                /sample),min=zran[0],max=zran[1])
          dzran = zran[1]-zran[0]
          cgimage,mapscl,/keep,oposition=p,layout=[npx,npy,j*npx+i+1]
          cgplot,[0],xsty=5,ysty=5,xran=xran,yran=yran,position=p,$
                 /nodata,/noerase
          cgaxis,xaxis=0,xran=xran_kpc,/xsty
          cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
          cgaxis,yaxis=0,yran=yran_kpc,/ysty
          cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty

;       cbform = '(D0.1)'
;       ncbdiv = 5

       endfor

    endfor

    cgps_close

  endforeach

end
