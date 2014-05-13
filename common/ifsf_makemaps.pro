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
;                       for different components; cleared up some floating point
;                       errors
;      2014apr21, DSNR, added line ratio maps and VO plots
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
pro ifsf_makemaps,initproc,comprange=comprange

   fwhm2sig = 2d*sqrt(2d*alog(2d))
   plotquantum = 2.5 ; in inches
   bad = 1d99
   c_kms = 299792.458d
   ncbdivmax = 7

;  Get fit initialization
   initdat=call_function(initproc)

;  Get galaxy-specific parameters from initialization file
   center_axes = -1
   center_nuclei = -1
   hasrangepar=0
   hasratpar=0
   if tag_exist(initdat,'argsmakemap') then begin
      if tag_exist(initdat.argsmakemap,'center_axes') then $
         center_axes = initdat.argsmakemap.center_axes
      if tag_exist(initdat.argsmakemap,'center_nuclei') then $
         center_nuclei = initdat.argsmakemap.center_nuclei
      if tag_exist(initdat.argsmakemap,'argslinratmaps') then begin
         argslinratmaps = initdat.argsmakemap.argslinratmaps
         hasratpar=1
      endif
      if tag_exist(initdat.argsmakemap,'rangefile') then begin
         rangefile = initdat.argsmakemap.rangefile
         hasrangepar=1
      endif
   endif

;  Get linelist
   linelist = ifsf_linelist(initdat.lines)

;  Get range file
;
;  plot types, in order; used for correlating with input ranges (array 
;  rangequant)
   plottypes = ['flux','velocity','sigma']
   hasrangefile=0
   if hasrangepar then begin
      if file_test(rangefile) then begin
         readcol,rangefile,rangeline,rangecomp,rangequant,rangelo,rangehi,$
         rangencbdiv,format='(A,I,A,D,D,I)',/silent
         hasrangefile=1
      endif else print,'IFSF_MAKEMAPS: Range file not found.'
   endif
   if keyword_set(rangefile) then begin
      if file_test(rangefile) then begin
         readcol,rangefile,rangeline,rangecomp,rangequant,rangelo,rangehi,$
            rangencbdiv,format='(A,I,A,D,D,I)',/silent
         hasrangefile=1
      endif else print,'IFSF_MAKEMAPS: Range file not found.'
   endif

;  Restore line maps
   if not tag_exist(initdat,'outlines') then outlines = linelist->keys() $
   else outlines = initdat.outlines
   restore,file=initdat.outdir+initdat.label+'.lin.xdr'

   size_tmp = size(linmaps[outlines[0]])
   dx = size_tmp[1]
   dy = size_tmp[2]
   if center_axes[0] eq -1 then center_axes = [double(dx)/2d,double(dy)/2d]
   if center_nuclei[0] eq -1 then center_nuclei = center_axes

;  Luminosity and angular size distances
   ldist = lumdist(initdat.zsys_gas,H0=73,Omega_m=0.27,Lambda0=0.73,/silent)
   kpc_per_as = ldist/(1+initdat.zsys_gas^2)*1000d/206265d

;  coordinates in kpc
   xran_kpc = double([-(center_axes[0]-0.5),dx-(center_axes[0]-0.5)]) $
              * initdat.platescale * kpc_per_as
   yran_kpc = double([-(center_axes[1]-0.5),dy-(center_axes[1]-0.5)]) $
              * initdat.platescale * kpc_per_as
   center_nuclei_kpc_x = (center_nuclei[0,*] - center_axes[0]) $
                         * initdat.platescale * kpc_per_as  
   center_nuclei_kpc_y = (center_nuclei[1,*] - center_axes[1]) $
                         * initdat.platescale * kpc_per_as

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plots of individual emission lines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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
             cgplot,[0],xsty=5,ysty=5,position=truepos,$
                    /nodata,/noerase,title=title
;            Plot axes in kpc
             cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
             cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
             cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
             cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
;            Plot nuclei
             cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
;            Colorbar
             cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
             ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
                                (dzran - zran[1]),format=cbform)
             cgcolorbar,position=cbpos,divisions=ncbdiv,$
                        ticknames=ticknames,/ver,/right,charsize=0.6

           endif

        endfor

     endfor
 
     cgps_close

   endforeach

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plots of line ratios
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if hasratpar then begin

;     Compute line ratios
      linrats = ifsf_lineratios(linmaps,linelist)

      keys = argslinratmaps->keys()

;     Loop through plot files
      foreach plt,keys do begin

;        Size of plot grid
         arrsiz = size(argslinratmaps[plt])
         npx = arrsiz[1]
         if arrsiz[0] gt 1 then begin
            npy = arrsiz[2]
            nplots = npx*npy
         endif else begin
            npy = 1
            nplots = npx
         endelse

         cgps_open,initdat.mapdir+initdat.label+plt+'.eps',charsize=1,/encap,$
            /inches,xs=plotquantum*npx,ys=plotquantum*npy,/qui
         pos = cglayout([npx,npy],ixmar=[3,3],iymar=[3,3],oxmar=[0,0],oymar=[0,0],$
            xgap=0,ygap=0)
         cbform = '(D0.2)'

;        Loop through plot panels
         for i=0,nplots-1 do begin

            tmpstr = strsplit(argslinratmaps[plt,i],'_',/extract)
            vcomp = fix(tmpstr[0])
            ptype = tmpstr[1]

            map = linrats[ptype[0],*,*,vcomp-1]
            ibd = where(map eq bad AND ~ finite(map),ctbd)
            inan = where(~finite(map),ctnan)
            igd = where(map ne bad AND finite(map),ctgd)

;           Line ratio maps
            if n_elements(tmpstr) eq 2 AND ctgd gt 0 then begin
               
               hasrange = 0
               if hasrangefile then begin
                  ithisline = where(ptype eq rangeline AND $
                     vcomp eq rangecomp,ctthisline)
                  if ctthisline eq 1 then begin
                     zran = [rangelo[ithisline],rangehi[ithisline]]
                     dzran = zran[1]-zran[0]
                     ncbdiv = rangencbdiv[ithisline]
                     ncbdiv = ncbdiv[0]
                     hasrange = 1
                  endif
               endif
               if keyword_set(comprange) AND ~hasrange then begin
                  zran = [min(map[igd]),max(map[igd])]
                  dzran = zran[1] - zran[0]
                  ncbdiv = ifsf_cbdiv(zran,0.5,7)
               endif

               title='c'+string(vcomp,format='(I0)')+' '
               if ptype eq 'n2ha' then title+=textoidl('[NII]/H\alpha')
               if ptype eq 'o3hb' then title+=textoidl('[OIII]/H\beta')
               if ptype eq 'ebv' then title+=textoidl('E(B-V)')

               if ctnan gt 0 then map[inan] = bad
               mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
                               min=zran[0],max=zran[1])

               cgloadct,65,/reverse
               cgimage,mapscl,/keep,pos=pos[*,i],opos=truepos,$
                       noerase=i ne 0,missing_value=bad,missing_index=255,$
                       missing_color='white'
               cgplot,[0],xsty=5,ysty=5,position=truepos,$
                      /nodata,/noerase,title=title
               cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
               cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
               cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
               cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
               cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
               cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
               ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
                                  (dzran - zran[1]),format=cbform)
               cgcolorbar,position=cbpos,divisions=ncbdiv,$
                  ticknames=ticknames,/ver,/right,charsize=0.6

;           VO plots
            endif else if n_elements(tmpstr) eq 4 AND ctgd gt 0 then begin
               
               map2 = linrats[tmpstr[3],*,*,vcomp-1]
               ibd2 = where(map2 eq bad AND ~ finite(map2),ctbd2)
               inan2 = where(~finite(map2),ctnan2)
               igd2 = where(map2 ne bad AND finite(map2),ctgd2)


               if ctgd2 gt 0 then begin

                  title='c'+string(vcomp,format='(I0)')
                  ptype = ptype+'_vs_'+tmpstr[3]
                  if ptype eq 'n2ha_vs_o3hb' then begin
                     xran = [-1.99d,0.99d]
                     yran = [-1.19d,1.49d]
                     xkew1 = 0.05d*indgen(110)-5d
                     ykew1 = 0.61d / (xkew1-0.47d)+1.19d
                     xkew2 = xkew1
                     ykew2 = xkew1-xkew1-99d
                     xtit = textoidl('[NII]/H\alpha')
                     ytit = textoidl('[OIII]/H\beta')
                  endif

                  if ctnan gt 0 then map[inan] = bad
                  if ctnan2 gt 0 then map2[inan2] = bad

                  cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,pos=pos[*,i],$
                     /nodata,noerase=i ne 0,title=title,xtit=xtit,ytit=ytit,$
                     /iso
                  cgoplot,map,map2,psym=1
                  cgoplot,xkew1,ykew1
                  cgoplot,xkew2,ykew2
                  
               endif
               
            endif

         endfor

         cgps_close

      endforeach

   endif

end
