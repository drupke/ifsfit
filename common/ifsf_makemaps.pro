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
;    comprange: in, optional, type=byte
;      Adjust ranges separately by component, rather than simultaneously for all
;      components.
;    rangefile: in, optional, type=string
;      File containing unique range specifications.
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
;      2014apr15, DSNR, moved colorbar division code to subroutine IFSF_CBDIV;
;                       added ability to automatically fix ranges separately
;                       for different components; added ability to change ranges
;                       using an input file; cleared up some floating point
;                       errors
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
pro ifsf_makemaps,initproc,comprange=comprange,rangefile=rangefile

  fwhm2sig = 2d*sqrt(2d*alog(2d))
  plotquantum = 2.5 ; in inches
  bad = 1d99
  c_kms = 299792.458d
  ncbdivmax = 7

; Get fit initialization
  initdat=call_function(initproc)

; Get linelist
  linelist = ifsf_linelist(initdat.lines)

; Get range file
;
; plot types, in order; used for correlating with input ranges (array 
; rangequant)
  plottypes = ['flux','velocity','sigma']
  hasrangefile=0
  if keyword_set(rangefile) then begin
     if file_test(rangefile) then begin
        readcol,rangefile,rangeline,rangecomp,rangequant,rangelo,rangehi,$
        rangencbdiv,format='(A,I,A,D,D,I)',/silent
        hasrangefile=1
     endif else print,'IFSF_MAKEMAPS: Range file not found.'
  endif

; Restore line maps
  if not tag_exist(initdat,'outlines') then outlines = linelist->keys() $
  else outlines = initdat.outlines
  restore,file=initdat.outdir+initdat.label+'.lin.xdr'

  size_tmp = size(linmaps[outlines[0]])
  dx = size_tmp[1]
  dy = size_tmp[2]
  center = [double(dx)/2d,double(dy)/2d]

; Luminosity and angular size distances
  ldist = lumdist(initdat.zsys_gas,H0=73,Omega_m=0.27,Lambda0=0.73,/silent)
  kpc_per_as = ldist/(1+initdat.zsys_gas^2)*1000d/206265d

; scales in kpc
  xran_kpc = double([-(center[0]-0.5),dx-(center[0]-0.5)]) $
             * initdat.platescale * kpc_per_as
  yran_kpc = double([-(center[1]-0.5),dy-(center[1]-0.5)]) $
             * initdat.platescale * kpc_per_as

; Size of plot grid
  npx = initdat.maxncomp
  npy = 3

; Linemap indices to plot below: flux, wavelength (converted to velocity), sigma
  ilinmap = [0,2,3]

; Loop through emission lines
  foreach line,outlines do begin

;   Get syntax of linelabel right; otherwise call to DEVICE chokes
    linelab=line
    ilb = strpos(linelab,'[')
    if ilb ne -1 then $
       linelab = strmid(linelab,0,ilb)+'\'+strmid(linelab,ilb)
    irb = strpos(linelab,']')
    if irb ne -1 then $
       linelab = strmid(linelab,0,irb)+'\'+strmid(linelab,irb)

    cgps_open,initdat.mapdir+initdat.label+linelab+'.eps',charsize=1,/encap,$
              /inches,xs=plotquantum*npx,ys=plotquantum*npy,/qui
    pos = cglayout([npx,npy],ixmar=[3,3],iymar=[3,3],oxmar=[0,0],oymar=[0,0],$
                   xgap=0,ygap=0)

;   loop through plot types
    for j=0,2 do begin

;      Set up colorbar labeling
       if j eq 0 then cbform = '(D0.1)' else cbform = '(I0)'


;      Set up ranges for all components at once
       if ~ keyword_set(comprange) then begin
        
          mapallcomp = linmaps[line,*,*,0:initdat.maxncomp-1,ilinmap[j]]
          ibd = where(mapallcomp eq bad,ctbd)
          igd = where(mapallcomp ne bad AND mapallcomp ne 0,ctgd)
          if j eq 1 then begin
;            redshift with respect to galaxy systemic
             zdiff = mapallcomp[igd]/linelist[line]-1d - initdat.zsys_gas
;            relativistic velocity shift;
;            see http://hyperphysics.phy-astr.gsu.edu/hbase/relativ/reldop2.html
             mapallcomp[igd] = c_kms * ((zdiff+1d)^2d - 1d) / ((zdiff+1d)^2d + 1d)
          endif
          zran = [min(mapallcomp[igd]),max(mapallcomp[igd])]
          if j eq 0 then begin
             zmax_flux = zran[1]
             zran=[0,1]
             dzran = 1
             ncbdiv = 5
          endif else begin
             divarr = ifsf_cbdiv(zran,100d,ncbdivmax)
             ncbdiv = divarr[0]
             dzran = zran[1]-zran[0]
          endelse

       endif

;      Loop through velocity components
       for i=0,initdat.maxncomp-1 do begin

;         Plot index
          iplot = j*npx+i

;         Get map and scale
          map = linmaps[line,*,*,i,ilinmap[j]]
          ibd = where(map eq bad AND ~ finite(map),ctbd)
          inan = where(~finite(map),ctnan)
          igd = where(map ne bad AND map ne 0 AND finite(map),ctgd)

          if ctgd gt 0 then begin
            
             if j eq 1 then begin
                zdiff = map[igd]/linelist[line]-1d - initdat.zsys_gas
                map[igd] = c_kms * ((zdiff+1d)^2d - 1d) / ((zdiff+1d)^2d + 1d)
             endif

;            Set up ranges for each component separately

;            Check for manual range first ...
             hasrange = 0
             if hasrangefile then begin
                ithisline = where(line eq rangeline AND $
                                  i+1 eq rangecomp AND $
                                  plottypes[j] eq rangequant,ctthisline)
                if ctthisline eq 1 then begin
                   zran = [rangelo[ithisline],rangehi[ithisline]]
                   dzran = zran[1]-zran[0]
                   ncbdiv = rangencbdiv[ithisline]
                   ncbdiv = ncbdiv[0]
                   hasrange = 1
                endif
             endif
;            otherwise set it automagically.
             if keyword_set(comprange) AND ~hasrange then begin
                zran = [min(map[igd]),max(map[igd])]
                if j eq 0 then begin
                   zmax_flux = zran[1]
                   zran=[0,1]
                   dzran = 1
                   ncbdiv = 5
                endif else begin
                   divarr = ifsf_cbdiv(zran,100d,ncbdivmax)
                   ncbdiv = divarr[0]
                   dzran = zran[1]-zran[0]
                endelse
             endif

             if j eq 0 then map[igd] = map[igd]/zmax_flux
             if ctnan gt 0 then map[inan] = bad
             mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
                             min=zran[0],max=zran[1])

;            Plot image
             if j eq 0 then begin
                cgloadct,65,/reverse
                title='c'+string(i+1,format='(I0)')+' flux'
                title += ' ('+string(zmax_flux,format='(E0.2)')+')'
             endif
             if j eq 1 then begin
                cgloadct,74,/reverse
                title='c'+string(i+1,format='(I0)')+' velocity'
             endif
             if j eq 2 then begin
                cgloadct,65,/reverse
                title='c'+string(i+1,format='(I0)')+' sigma'
             endif
             cgimage,mapscl,/keep,pos=pos[*,iplot],opos=truepos,$
                     noerase=iplot ne 0,missing_value=bad,missing_index=255,$
                     missing_color='white'
;            Plot axes in kpc
             cgplot,[0],xsty=5,ysty=5,xran=xran,yran=yran,position=truepos,$
                    /nodata,/noerase,title=title
             cgaxis,xaxis=0,xran=xran_kpc,/xsty
             cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
             cgaxis,yaxis=0,yran=yran_kpc,/ysty
             cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
;            Colorbar
             cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
             ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
                                (dzran - zran[1]),format=cbform)
             cgcolorbar,position=cbpos,divisions=ncbdiv,$
                        ticknames=ticknames,/ver,/right

          endif

       endfor

    endfor

    cgps_close

  endforeach

end
