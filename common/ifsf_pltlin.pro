; docformat = 'rst'
;
;+
;
; Plot emission line fit and output to JPG.
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
;    pltpar: in, required, type=structure
;      Contains parameters to control plot. Tags:
;      label: type=strarr(Nlines), line labels for plot
;      wave: type=dblarr(Nlines), rest wavelengths of lines
;      linoth: type=dblarr(Notherlines,Ncomp), wavelengths of other
;              lines to plot
;      nx: type=int, # of plot columns
;      ny: type=int, # of plot rows
;      
;    outfile: in, required, type=string
;      Full path and name of output plot.
;
; :Keywords:
;    micron: in, optional, type=byte
;      Label output plots in um rather than A. Input wavelengths still assumed
;      to be in A.
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
;      13sep12, DSNR, re-written
;      2013oct, DSNR, documented
;      2013nov21, DSNR, renamed, added license and copyright 
;      2015may13, DSNR, switched from using LAYOUT keyword to using CGLAYOUT
;                       procedure to fix layout issues
;      2016aug31, DSNR, added overplotting of continuum ranges masked during
;                       continuum fit with thick cyan line
;      2016sep13, DSNR, added MICRON keyword
;    
; :Copyright:
;    Copyright (C) 2013--2016 David S. N. Rupke
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
pro ifsf_pltlin,instr,pltpar,outfile

   set_plot,'Z'
   device,decomposed=0,set_resolution=[1500,900],set_pixel_depth=24
   !P.charsize=1.5
   !P.charthick=1
   !P.thick=4
   erase

   defaultXtickint=!X.tickinterval
   defaultXminor=!X.minor
   !X.tickinterval=50
   !X.minor=10
   if tag_exist(pltpar,'micron') then begin
      !X.tickinterval /= 0.5d4
      !X.minor /= 0.5d4
   endif
   pos = cglayout([pltpar.nx,pltpar.ny],$ ; ixmar=[5d,0d],iymar=[-5d,0d],$
                  oxmar=[22,0],oymar=[10,0],xgap=20,ygap=6)

  ncomp = instr.param[1]
  colors = ['Magenta','Green','Orange','Teal']

  wave = instr.wave
  spectot = instr.spec
  specstars = instr.cont_dat
  speclines = instr.emlin_dat
  modstars = instr.cont_fit
  modlines = instr.emlin_fit
  modtot = modstars + modlines

;  norm = max(modstars)
;  spectot /= norm
;  specstars /= norm
;  speclines /= norm
;  modtot /= norm
;  modstars /= norm
;  modlines /= norm

  zbase = instr.zstar

   if tag_exist(pltpar,'micron') then wave /= 10000d

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

  nlin = n_elements(pltpar.label)
  linlab = pltpar.label
  linwav = pltpar.wave
  off = pltpar.off
  if tag_exist(pltpar,'linoth') then linoth = pltpar.linoth $
  else linoth = strarr(1,nlin)
  for i=0,nlin-1 do begin

     linwavtmp = linwav[i]
     xran = (linwavtmp[0] + off[*,i]) * (1d + zbase)
     ind = where(wave gt xran[0] AND wave lt xran[1],ct)

     cgplot,[0],/nodata,xsty=4,ysty=4,pos=pos[*,i],noerase=i ne 0,backg='Black'
     xwin = [pos[0,i],pos[2,i]]
     ywin = [pos[1,i],pos[3,i]]
     dxwin = xwin[1]-xwin[0]
     dywin = ywin[1]-ywin[0]

     if ct gt 0 then begin
        pos_fit = [xwin[0],ywin[0]+0.3*dywin,$
                   xwin[1],ywin[1]]
        ydat = spectot
        ymod = modtot
        yran = [min([ydat[ind],ymod[ind]]),max([ydat[ind],ymod[ind]])]
        icol = double(i)/double(pltpar.nx)
        if icol eq fix(icol) then ytit = 'Fit' else ytit = ''
        cgplot,wave,ydat,xran=xran,yran=yran,pos=pos_fit,$
               xtickn=replicate(' ',60),ytit=ytit,/noerase,$
               axiscol='White',col='White',/norm,/xsty,/ysty,thick=1
        cgoplot,wave,ymod,color='Red',thick=4
        for j=1,ncomp do begin
           flux = ifsf_cmplin(instr,linlab[i],j,/velsig)
           cgoplot,wave,yran[0]+flux,color=colors[j-1],linesty=2,thick=2
           if linoth[0,i] ne '' then begin
              for k=0,n_elements(linoth[*,i])-1 do begin
                 if linoth[k,i] ne '' then begin
                    flux = ifsf_cmplin(instr,linoth[k,i],j,/velsig)
                    cgoplot,wave,yran[0]+flux,color=colors[j-1],linesty=2,$
                            thick=2
                 endif
              endfor
           endif
        endfor
        cgtext,xran[0]+(xran[1]-xran[0])*0.05d,$
               yran[0]+(yran[1]-yran[0])*0.85d,$
               linlab[i],charsize=1.5,charthick=2,/dat
        if nmasked gt 0 then $
           for j=0,nmasked-1 do $
              cgoplot,[masklam[0,j],masklam[1,j]],[yran[0],yran[0]],thick=8,$
                      color='Cyan'
        pos_res = [xwin[0],ywin[0],$
                   xwin[1],ywin[0]+0.3*dywin]
        ydat = specstars
        ymod = modstars
        yran = [min([ydat[ind],ymod[ind]]),max([ydat[ind],ymod[ind]])]
        if icol eq fix(icol) then ytit = 'Residual' else ytit = ''
        cgplot,wave,ydat,xran=xran,yran=yran,/noerase,ytit=ytit,$
               axiscol='White',col='White',/norm,pos=pos_res,/xsty,/ysty,thick=1
        cgoplot,wave,ymod,color='Red',thick=4
     endif

  endfor

  if tag_exist(pltpar,'micron') then $
    xtit = 'Observed Wavelength ($\mu$m)' $
  else $
     xtit = 'Observed Wavelength (!3' + STRING(197B) + '!X)'
  cgtext,0.5,0.02,xtit,/norm,align=0.5,charsize=2,charthick=2
  
  tmpfile = outfile
  img = cgsnapshot(filename=tmpfile,/jpeg,/nodialog,quality=100)

  !X.tickinterval=defaultXtickint
  !X.minor=defaultXminor
  
end
