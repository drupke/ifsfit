; docformat = 'rst'
;
;+
;
; Spatially smooth initial guesses in # components, linewidths, and redshifts.
; Start by median-smoothing maps with FILTER_IMAGE (default box size of 3), and 
; then dilating the result by the following kernel:
; 
;   0 0 1 0 0
;   0 1 1 1 0
;   1 1 1 1 1
;   0 1 1 1 0
;   0 0 1 0 0
;
; :Categories:
;    IFSFIT
;
; :Returns:
;
; :Params:
;   infile: in, required, type=string
;     Name of IDL save file produced by IFSFA (*.lininit.xdr), containing 
;     number of emission-line components in each spaxel and linewidth sigma
;     and redshift in each velocity component of each spaxel.
;   outfile: in, required, type=string
;     Name of output IDL save file to store spatially-smoothed arrays of 
;     # of components, linewidths, and redshifts.
;     
; :Keywords:
;   plot: in, optional, type=boolean
;     Produce plots to screen of original and spatially-smoothed maps
;   ywinsize: in, optional, type=integer, default=600
;     Y size of output window in dpi (?).
;   smbox: in, optional, type=integer, default=3
;     Size of box for moving median filter of maps.
; 
; :Author:
;    David Rupke
;
; :History:
;    Change History::
;      2021jun10, DSNR, created
;      2022mar23, DSNR, commented
;
; :Copyright:
;    Copyright (C) 2021-2022 David S. N. Rupke
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
pro ifsf_smoothinit,infile,outfile,plot=plot,ywinsize=ywinsize,smbox=smbox

   bad = 1d99

   if not keyword_set(ywinsize) then ywinsize=double(600) $
      else ywinsize=double(ywinsize)
   if not keyword_set(smbox) then smbox=3

   restore,infile
      
   if keyword_set(plot) then begin

      cgdelete,/all
      set_plot,'x'
      cleanplot,/silent

      ; hashes to hold window ids and titles
      wids = hash()
      wtis = hash()

      ; cycle through lines
      foreach line, emlncomp.keys() do begin

         ; size of input cube
         cubesize = size(emlncomp[line])
         ; max. number of velocity components
         maxncomp = max(emlncomp[line])

         waspect_ncomp = double(cubesize[2]/cubesize[1])/2d
         waspect_zsig = double(cubesize[2]/cubesize[1])/3d*double(maxncomp)*1.2

         ; three windows for each line: ncomp, z, sig
         cgwindow,waspect=waspect_ncomp,$
            wysize=ywinsize/2d,wtitle=line+': Ncomp',wypos=0d,wxpos=0d
         cgwindow,waspect=waspect_zsig,$
            wysize=ywinsize,wtitle=line+': z',wypos=ywinsize/2d + 100d
         cgwindow,waspect=waspect_zsig,$
            wysize=ywinsize,wtitle=line+': sig',wxpos=ywinsize/waspect_zsig

         wids[line] = cgQuery(TITLE=titles)
         wtis[line] = titles
 
      endforeach

   endif
   
   ; hashes to hold output z values and numbers of components
   newemlncomp = hash()
   newemlz = hash()
   newemlsiginit = hash()

   ; kernel for dilation
   str_el = [[0b,0b,1b,0b,0b],$
      [0b,1b,1b,1b,0b],$
      [1b,1b,1b,1b,1b],$
      [0b,1b,1b,1b,0b],$
      [0b,0b,1b,0b,0b]]

   ; cycle through lines
   foreach line, emlncomp.keys() do begin

      ; size of input cube
      cubesize = size(emlncomp[line])
      ; max. number of velocity components
      maxncomp = max(emlncomp[line])
      ; array to hold output z values
      newemlz[line] = dblarr(cubesize[1],cubesize[2],maxncomp)+bad
      newemlsiginit[line] = dblarr(cubesize[1],cubesize[2],maxncomp)+bad
      newemlncomp[line] = intarr(cubesize[1],cubesize[2])

      ; plot original number of copmonents
      if keyword_set(plot) then begin
         mapscl = cgimgscl(emlncomp[line],minval=0,maxval=maxncomp,$
            ncolors=maxncomp+1)
         cgset,wids[line,0]
         cgloadct,8,/brewer,ncolors=maxncomp+1,/window
         cgimage,mapscl,/axes,axkey = {ticklen:1.},/keep,layout=[2,1,1],$
            multimargin=[10,1,1,1],chars=2,/window
      endif
      
      ; loop through number of components for this line
      for i=0,max(emlncomp[line])-1 do begin

         ; do smoothing and dilating

         ; z map for this line and component
         mapz = emlz[line,*,*,i]
         mapsig = emlsiginit[line,*,*,i]
         ; where there are measured z values
         ibd = where(mapz eq bad)
         igd = where(mapz ne bad)
         ; range of measured values; for z, round min/max to nearest thousandth and
         ; add some padding; for sigma, round after dividing by 10
         zran = [min(mapz[igd]),max(mapz[igd])]
         zran = double(zran)*1000d
         zran = [floor(zran[0]),ceil(zran[1])]/1000d
         dzran = zran[1]-zran[0]
         sigran = [min(mapsig[igd]),max(mapsig[igd])]
         sigran = double(sigran)/10d
         sigran = [floor(sigran[0]),ceil(sigran[1])]*10d
         dsigran = sigran[1]-sigran[0]
         ; byte scaled raw map
         mapzscl = cgimgscl(mapz,minval=zran[0], maxval=zran[1], $
            missing_index=0b, missing_value=bad)
         mapsscl = cgimgscl(mapsig,minval=sigran[0], maxval=sigran[1], $
            missing_index=0b, missing_value=bad)
;         mapscl[ibd] = 0b
         ; smooth map with moving median filter
         mapzsm = filter_image(mapz,med=smbox,/all)
         mapssm = filter_image(mapsig,med=smbox,/all)
         ; then byte scale smoothed map
         mapzsclsm = cgimgscl(mapzsm,minval=zran[0],maxval=zran[1], $
            missing_index=0b, missing_value=bad)
         mapssclsm = cgimgscl(mapssm,minval=sigran[0],maxval=sigran[1], $
            missing_index=0b, missing_value=bad)

         ; dilate byte-scaled map using kernel str_el
         ; input must be byte-scaled or boolean
         ; /gray for byte-scaled case
         ; /constrained to only dilate into unfit points
         mapzscldil = dilate(mapzsclsm,str_el,/gray,/constrained)
         mapsscldil = dilate(mapssclsm,str_el,/gray,/constrained)
         ; scale back to decimal; use the inverse of the function in BYTSCL
         newemlz[line,*,*,i] = zran[0] + (double(mapzscldil)/255.9999d)*dzran
         newemlsiginit[line,*,*,i] = sigran[0] + (double(mapsscldil)/255.9999d)*dsigran
         

         ; calculate number of new components in each spaxel
         igddil = where(mapzscldil gt 0b)
         comptmp = newemlncomp[line]
         comptmp[igddil] += 1
         newemlncomp[line] = comptmp
;         
         ; plot everything

         if keyword_set(plot) then begin

            ; z

            ; unsmoothed image
            cgset,wids[line,1]
            cgloadct,74,/reverse,/window
            cgimage,mapzscl,/axes,axkey = {ticklen:1.},/keep,/add,$
               layout=[3,maxncomp,1+(3*i)],multimargin=[10,1,1,1],$
               chars=2,/window
            ; for some reason oposition keyword in cgimage not working with
            ; call to window; plan B:
            common fsc_$cgimage
            truepos = _cgimage_position
            ncbdiv = 3
            cbform = '(D0.3)'
            cbpos = [truepos[0],truepos[1]-0.1*(truepos[3]-truepos[1]),$
               truepos[2],truepos[1]-0.05*(truepos[3]-truepos[1])]
            ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
               (dzran - zran[1]),format=cbform)
            cgcolorbar,/addcmd,pos=cbpos,ticknames=ticknames,div=ncbdiv
            ; smoothed map
            cgimage,mapzsclsm,/axes,axkey = {ticklen:1.},/keep,$
               layout=[3,maxncomp,2+(3*i)],multimargin=[10,1,1,1],/addcmd,$
               chars=2
            ; dilated map
            cgimage,mapzscldil,/axes,axkey = {ticklen:1.},/keep,$
               layout=[3,maxncomp,3+(3*i)],multimargin=[10,1,1,1],/addcmd,$
               chars=2

            ; sigma

            ; unsmoothed image
            cgset,wids[line,2]
            cgloadct,74,/reverse,/window
            cgimage,mapsscl,/axes,axkey = {ticklen:1.},/keep,$
               layout=[3,maxncomp,1+(3*i)],multimargin=[10,1,1,1],/addcmd,$
               chars=2,/window
            ; for some reason oposition keyword in cgimage not working with
            ; call to window; plan B:
            common fsc_$cgimage
            truepos = _cgimage_position
            ncbdiv = 3
            cbform = '(I0)'
            cbpos = [truepos[0],truepos[1]-0.1*(truepos[3]-truepos[1]),$
               truepos[2],truepos[1]-0.05*(truepos[3]-truepos[1])]
            ticknames = string(dindgen(ncbdiv+1)*dsigran/double(ncbdiv) - $
               (dsigran - sigran[1]),format=cbform)
            cgcolorbar,/addcmd,pos=cbpos,ticknames=ticknames,div=ncbdiv
            ; smoothed map
            cgimage,mapssclsm,/axes,axkey = {ticklen:1.},/keep,$
               layout=[3,maxncomp,2+(3*i)],multimargin=[10,1,1,1],/addcmd,$
               chars=2
            ; dilated map
            cgimage,mapsscldil,/axes,axkey = {ticklen:1.},/keep,$
               layout=[3,maxncomp,3+(3*i)],multimargin=[10,1,1,1],/addcmd,$
               chars=2

         endif

      endfor

      ; plot # of new components
      if keyword_set(plot) then begin
         mapscl = cgimgscl(newemlncomp[line],minval=0,maxval=maxncomp,$
            ncolors=maxncomp+1)
         cgset,wids[line,0]
         cgloadct,8,/brewer,ncolors=maxncomp+1,/window
         cgimage,mapscl,/axes,axkey = {ticklen:1.},/keep,layout=[2,1,2],$
            multimargin=[10,1,1,1],chars=2,/window,/addcmd
      endif

   endforeach
   
   save,newemlz,newemlncomp,newemlsiginit,file=outfile
   
end
