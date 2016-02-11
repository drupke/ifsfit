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
;    
; :Copyright:
;    Copyright (C) 2013 David S. N. Rupke
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
pro ifsf_pltcont,instr,outfile,ps=ps

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
     !P.charsize=1
     !P.charthick=1
     erase
  endelse
  
  wave = instr.wave
  spectot = instr.spec
  specstars = instr.cont_dat
  speclines = instr.emlin_dat
  modstars = instr.cont_fit
  modlines = instr.emlin_fit
  modtot = modstars + modlines

  if keyword_set(outstelfit) then outstelfit=modstars

  norm = max(modstars)
  spectot /= norm
  specstars /= norm
  speclines /= norm
  modtot /= norm
  modstars /= norm
  modlines /= norm
  
  xran = instr.fitran
  dxran = xran[1] - xran[0]
  xran1 = [xran[0],xran[0]+dxran/3d]
  xran2 = [xran[0]+dxran/3d,xran[0]+2d*dxran/3d]
  xran3 = [xran[0]+2d*dxran/3d,xran[1]]
  i1 = where(wave gt xran1[0] AND wave lt xran1[1],ct1)
  i2 = where(wave gt xran2[0] AND wave lt xran2[1],ct2)
  i3 = where(wave gt xran3[0] AND wave lt xran3[1],ct3)
  
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
  ytit = textoidl('Normalized F_\lambda')
  multiplot,[1,3],/doyaxis,/doxaxis,ygap=0.02,mxtitle=xtit,mytitle=ytit,$
            mxtitsize=2,mytitsize=2,mxtitoff=1
  ydat = specstars
  ymod = modstars
  yran = [min([ydat[i1],ymod[i1]]),max([ydat[i1],ymod[i1]])]
  ydi = ydat[i1]
  ymodi = ymod[i1]
  y = [ydi-ymodi]
  ny = n_elements(y)
  iysort = sort(y)
  ysort = y[iysort]
  ymodisort = ymodi[iysort]
  if ysort[ny-ntop] lt ysort[ny-1]*maxthresh then $
     yran[1] = max(ysort[0:ny-ntop]+ymodisort[0:ny-ntop])
  if ysort[nbottom] gt ysort[0]*maxthresh then $
     yran[0] = min(ysort[nbottom:ny-1]+ymodisort[nbottom:ny-1])
  if (yran[0] lt 0) then yran[0]=0
  cgplot,wave,ydat,xran=xran1,yran=yran,/xsty,/ysty,$
         color='White',axiscol='White'
  cgoplot,wave,ymod,color='Red'
  multiplot,/doyaxis,/doxaxis
  ydat = specstars
  ymod = modstars
  yran = [min([ydat[i2],ymod[i2]]),max([ydat[i2],ymod[i2]])]
  ydi = ydat[i2]
  ymodi = ymod[i2]
  y = [ydi-ymodi]
  ny = n_elements(y)
  iysort = sort(y)
  ysort = y[iysort]
  ymodisort = ymodi[iysort]
  if ysort[ny-ntop] lt ysort[ny-1]*maxthresh then $
     yran[1] = max(ysort[0:ny-ntop]+ymodisort[0:ny-ntop])
  if ysort[nbottom] gt ysort[0]*maxthresh then $
     yran[0] = min(ysort[nbottom:ny-1]+ymodisort[nbottom:ny-1])
  if (yran[0] lt 0) then yran[0]=0
  cgplot,wave,ydat,xran=xran2,yran=yran,/xsty,/ysty,$
         color='White',axiscol='White'
  cgoplot,wave,ymod,color='Red'
  if ct3 gt 0 then begin
     multiplot,/doyaxis,/doxaxis
     ydat = specstars
     ymod = modstars
     yran = [min([ydat[i3],ymod[i3]]),max([ydat[i3],ymod[i3]])]
     ydi = ydat[i3]
     ymodi = ymod[i3]
     y = [ydi-ymodi]
     ny = n_elements(y)
     iysort = sort(y)
     ysort = y[iysort]
     ymodisort = ymodi[iysort]
     if ysort[ny-ntop] lt ysort[ny-1]*maxthresh then $
        yran[1] = max(ysort[0:ny-ntop]+ymodisort[0:ny-ntop])
     if ysort[nbottom] gt ysort[0]*maxthresh then $
        yran[0] = min(ysort[nbottom:ny-1]+ymodisort[nbottom:ny-1])
     if (yran[0] lt 0) then yran[0]=0
     cgplot,wave,ydat,xran=xran3,yran=yran,/xsty,/ysty,$
          color='White',axiscol='White'
     cgoplot,wave,ymod,color='Red'
  endif else multiplot
  multiplot,/reset

  tit = 'STELLAR CONTINUUM FIT'
  cgtext,0.5,0.96,tit,/norm,align=0.5,charsize=2,charthick=2

  tmpfile = outfile
  if dops then device,/close_file $
  else img = cgsnapshot(filename=tmpfile,/jpeg,/nodialog,quality=100)
  
end
