; docformat = 'rst'
;
;+
;
; Plot continuum fit and output to JPG.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    None.
;
; :Params:
;    instr: in, required, type=structure
;      Contains results of fit.
;    outfile: in, required, type=string
;      Full path and name of output plot.
;
; :Keywords:
;    ps: in, optional, type=byte, default=0
;      Set to get postscript, instead of default jpg, output.
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
;      2009, DSNR, created
;      2013oct09, DSNR, documented
;      2013nov21, DSNR, renamed, added license and copyright
;      2014jan16, DSNR, changed to Coyote Graphics routines
;      2016aug31, DSNR, added overplotting of continuum ranges masked during
;                       continuum fit with thick cyan line
;      2019jan24, DSNR, added option to have yrange default to min/max values
;                       rather than 0 to max
;    
; :Copyright:
;    Copyright (C) 2013--2019 David S. N. Rupke
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
pro ifsf_pltcont,instr,outfile,compspec=compspec,comptitles=comptitles,$
                 ps=ps,title=title,fitran=fitran,yranminmax=yranminmax
                 

  if keyword_set(compspec) then begin
     sizecomp = size(compspec)
     if sizecomp[0] gt 1 then ncomp = sizecomp[2] else ncomp = 1
     compcolors = ['Cyan','Orange','Green']
  endif else ncomp = 0
  if keyword_set(ps) then dops=1 else dops=0

  if dops then begin
     set_plot,'ps',/copy,/interpolate
     device,filename=outfile+'.eps',/encapsulated,xsize=10,ysize=7.5,$
            bits_per_pixel=8,/color,/inches
     !P.charsize=1
     !P.charthick=2
     !P.thick=2
  endif else begin
     set_plot,'Z'
     device,decomposed=0,set_resolution=[1280,960],set_pixel_depth=24
     !P.charsize=1.5
     !P.charthick=1
     !P.thick=4
     cgerase,'Black'
  endelse

  defaultXtickint=!X.tickinterval
  defaultXminor=!X.minor
  !X.tickinterval=200
  !X.minor=50
  
  wave = instr.wave
  specstars = instr.cont_dat
  speclines = instr.emlin_dat
  modstars = instr.cont_fit
  
  if keyword_set(fitran) then xran = fitran else xran = instr.fitran
  dxran = xran[1] - xran[0]
  xran1 = [xran[0],xran[0]+dxran/3d]
  xran2 = [xran[0]+dxran/3d,xran[0]+2d*dxran/3d]
  xran3 = [xran[0]+2d*dxran/3d,xran[1]]
  i1 = where(wave gt xran1[0] AND wave lt xran1[1],ct1)
  i2 = where(wave gt xran2[0] AND wave lt xran2[1],ct2)
  i3 = where(wave gt xran3[0] AND wave lt xran3[1],ct3)

;  Find masked regions during continuum fit
   nmasked = 0 ; # of masked regions
;  Find consecutive unmasked regions
   consec,instr.ct_indx,lct,hct,nct
;  Set interior masked regions
   if nct gt 1 then begin
      nmasked = nct-1
      masklam = dblarr(2,nmasked)
      for i=0,nmasked-1 do $
         masklam[*,i] = [wave[instr.ct_indx[hct[i]]],$
                         wave[instr.ct_indx[lct[i+1]]]]
   endif
;  Set masked region if it occurs at beginning of lambda array
   if instr.ct_indx[0] ne 0 then begin
      nmasked++
      masklam = [[wave[0],wave[instr.ct_indx[lct[0]-1]]],[masklam]]
   endif
;  Set masked region if it occurs at end of lambda array
   if instr.ct_indx[n_elements(instr.ct_indx)-1] ne $
      n_elements(instr.wave)-1 then begin
      nmasked++
      masklam = [[masklam],$
                 [wave[instr.ct_indx[hct[nct-1]]],$
                  wave[n_elements(instr.wave)-1]]]
   endif

  maxthresh=0.2
  ntop = 20
  nbottom = 20
  if n_elements(wave) lt 100 then begin
    ntop = 10
    nbottom = 10
  endif
  ntop++
  nbottom--

  xtit = 'Observed Wavelength (!3' + STRING(197B) + '!X)'
;xtit = textoidl('Observed Wavelength (!6!sA!r!u!9 %!6!n )')
;  ytit = textoidl('F_\lambda')
  ytit=''
  multiplot,[1,3],/doyaxis,/doxaxis,ygap=0.02,mxtitle=xtit,mytitle=ytit,$
            mxtitsize=2,mytitsize=2,mxtitoff=1
  if ct1 gt 0 then begin
     ydat = specstars
     ymod = modstars
     if keyword_set(yranminmax) then $
        yran = [min([ydat[i1],ymod[i1]]),max([ydat[i1],ymod[i1]])] $
     else $
        yran = [0,max([ydat[i1],ymod[i1]])]
     ydi = ydat[i1]
     ymodi = ymod[i1]
     y = [ydi-ymodi]
     ny = n_elements(y)
     iysort = sort(y)
     ysort = y[iysort]
     ymodisort = ymodi[iysort]
     if ysort[ny-ntop] lt ysort[ny-1]*maxthresh then $
        yran[1] = max(ysort[0:ny-ntop]+ymodisort[0:ny-ntop])
;  if ysort[nbottom] gt ysort[0]*maxthresh then $
;     yran[0] = min(ysort[nbottom:ny-1]+ymodisort[nbottom:ny-1])
;  if (yran[0] lt 0) then yran[0]=0
     cgplot,wave,ydat,xran=xran1,yran=yran,/xsty,/ysty,$
            color='White',axiscol='White',thick=1,backg='Black'
     if ncomp gt 0 then $
        for i=0,ncomp-1 do $
           cgoplot,wave,compspec[*,i],color=compcolors[i],linesty=0,thick=2
     cgoplot,wave,ymod,color='Red'
     if nmasked gt 0 then $
        for j=0,nmasked-1 do $
           cgoplot,[masklam[0,j],masklam[1,j]],[yran[0],yran[0]],thick=8,$
                   color='Cyan'
  endif
  multiplot,/doyaxis,/doxaxis
  if ct2 gt 0 then begin
     ydat = specstars
     ymod = modstars
     if keyword_set(yranminmax) then $
        yran = [min([ydat[i2],ymod[i2]]),max([ydat[i2],ymod[i2]])] $
     else $
        yran = [0,max([ydat[i2],ymod[i2]])]
     ydi = ydat[i2]
     ymodi = ymod[i2]
     y = [ydi-ymodi]
     ny = n_elements(y)
     iysort = sort(y)
     ysort = y[iysort]
     ymodisort = ymodi[iysort]
     if ysort[ny-ntop] lt ysort[ny-1]*maxthresh then $
        yran[1] = max(ysort[0:ny-ntop]+ymodisort[0:ny-ntop])
;  if ysort[nbottom] gt ysort[0]*maxthresh then $
;     yran[0] = min(ysort[nbottom:ny-1]+ymodisort[nbottom:ny-1])
;  if (yran[0] lt 0) then yran[0]=0
     cgplot,wave,ydat,xran=xran2,yran=yran,/xsty,/ysty,$
            color='White',axiscol='White',thick=1
     if ncomp gt 0 then $
        for i=0,ncomp-1 do $
           cgoplot,wave,compspec[*,i],color=compcolors[i],linesty=0,thick=2
     cgoplot,wave,ymod,color='Red'
     if nmasked gt 0 then $
        for j=0,nmasked-1 do $
           cgoplot,[masklam[0,j],masklam[1,j]],[yran[0],yran[0]],thick=8,$
                   color='Cyan'
  endif
  multiplot,/doyaxis,/doxaxis
  if ct3 gt 0 then begin
     ydat = specstars
     ymod = modstars
     if keyword_set(yranminmax) then $
        yran = [min([ydat[i3],ymod[i3]]),max([ydat[i3],ymod[i3]])] $
     else $
        yran = [0,max([ydat[i3],ymod[i3]])]
     ydi = ydat[i3]
     ymodi = ymod[i3]
     y = [ydi-ymodi]
     ny = n_elements(y)
     iysort = sort(y)
     ysort = y[iysort]
     ymodisort = ymodi[iysort]
     if ysort[ny-ntop] lt ysort[ny-1]*maxthresh then $
        yran[1] = max(ysort[0:ny-ntop]+ymodisort[0:ny-ntop])
;     if ysort[nbottom] gt ysort[0]*maxthresh then $
;        yran[0] = min(ysort[nbottom:ny-1]+ymodisort[nbottom:ny-1])
;     if (yran[0] lt 0) then yran[0]=0
     cgplot,wave,ydat,xran=xran3,yran=yran,/xsty,/ysty,$
          color='White',axiscol='White',thick=1
     if ncomp gt 0 then $
        for i=0,ncomp-1 do $
           cgoplot,wave,compspec[*,i],color=compcolors[i],linesty=0,thick=2
     cgoplot,wave,ymod,color='Red'
     if nmasked gt 0 then $
        for j=0,nmasked-1 do $
           cgoplot,[masklam[0,j],masklam[1,j]],[yran[0],yran[0]],thick=8,$
                   color='Cyan'
  endif
  multiplot,/reset

  if keyword_set(title) then $
     cgtext,0.75,0.96,title,/norm,align=0.5,charsize=2,charthick=2
  if ncomp gt 0 AND keyword_set(comptitles) then begin
     comptitles = [comptitles,'total']
     compcolors = [compcolors[0:ncomp-1],'Red']
     al_legend,comptitles,linesty=dblarr(ncomp+1),colors=compcolors,$
               /horiz,pos=[0.05,0.99],/norm,linsiz=0.2
  endif

  tmpfile = outfile
  if dops then device,/close_file $
  else img = cgsnapshot(filename=tmpfile,/jpeg,/nodialog,quality=100)

  !X.tickinterval=defaultXtickint
  !X.minor=defaultXminor
  
end
