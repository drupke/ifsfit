; docformat = 'rst'
;
;+
;
; Plot multiplet data and fit, and output plot information into three files.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    None.
;
; :Params:
;    gal: in, required, type=string
;      Label specifying output location and filename.
;    continuum: in, required, type=dblarr(N)
;      Model continuum flux.
;    directoryname: in, required, type=string
;      Output directory.
;    flux: in, required, type=dblarr(N)
;      Fluxes.
;    outfile: in, required, type=string
;      Full path and name of output plot.
;    param: in, required, type=dblarr
;      Parameters output by Na D fitting routine.
;    relativeflux: in, required, type=dblarr(N)
;      Normalized fluxes.
;    wave: in, required, type=dblarr(N)
;      Wavelengths.
;    zsys: in, required, type=double
;      Systemic redshift.
;
; :Keywords:
;    xran: in, optional, type=dblarr(2), default=all
;      Wavlength range of output plot.
;    yran: in, optional, type=dblarr(2), default=[0\,1.5]
;      Normalized flux range of output plot.
;    init: in, optional, type=byte
;      Plot initial guess.
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
;      2010jul21, DSNR, created
;      2013nov21, DSNR, documented, renamed, added license and
;                       copyright
;      2014may15, DSNR, re-written for new fitting routines
;      2016xxxYY, AT, adopted for UV data
;      2016aug08, DSNR, added option to plot initial guess
;      2019nov20, DSNR, copied from IFSF_PLTDOUBLET
;
; :Copyright:
;    Copyright (C) 2013--2019 David S. N. Rupke, Anthony To
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
pro ifsf_pltmultiplet,gal,wave,relativeflux,continuum,flux,param,linelist,$
                      initstr,multiplet,directoryname,outfile,zsys,$
                      xran=xran,yran=yran,init=init,zres=zres,smooth=smooth

; Default wavelengths
  c = 299792.458d
  
  wavelines = linelist[initstr.reflinename]*(1d +zsys)
  wavenames = initstr.reflinename
  for i=0,n_elements(initstr.linenames)-1 do begin
     wavelines = [wavelines,linelist[initstr.linenames[i]]*(1d +zsys)]
     wavenames = [wavenames,initstr.linenames[i]]
  endfor
  moddoubletabs=1
  modflux = ifsf_multipletfcn(wave,param,refloglf=initstr.refloglf,$
                              refmultwave=initstr.refmultwave,$
                              loglf=initstr.loglf,multwave=initstr.multwave,$
                              modabs=modabs)
  size_abs = size(modabs)
  if size_abs[0] eq 1 then nabs = 1 $
  else if size_abs[0] eq 2 then nabs = fix(size_abs[2]) $
  else nabs=0l

; PLOT
  
  if keyword_set(init) then begin
     set_plot,'x'
     device,decomposed=0
  endif else begin
     set_plot,'Z'
     device,decomposed=0,set_resolution=[1280,960],set_pixel_depth=24
  endelse
  !P.charsize=1.5
  !P.charthick=1
  erase

  nplt = 1
  if ~ keyword_set(xran) then xran = [min(wave),max(wave)] $
  else begin
     xsize = size(xran)
     if xsize[0] gt 1 then nplt = xsize[2]
  endelse
  if ~ keyword_set(yran) then yran = [0d,1.5d]
  
  if nplt gt 1 then begin
    if nplt eq 2 OR nplt eq 3 then begin
       nx = nplt
       ny = 1
    endif else if nplt eq 4 then begin
       nx = 2
       ny = 2
    endif else if nplt eq 5 OR nplt eq 6 then begin
       nx = 3
       ny = 2
    endif
  endif else begin
     nx = 1
     ny = 1
  endelse
  pos = cglayout([nx,ny],oxmar=[12,2],oymar=[10,2],xgap=10,ygap=5)

  if keyword_set(smooth) then relfluxuse = smooth(relativeflux,smooth) $
  else relfluxuse = relativeflux

  for i=0,nplt-1 do begin
     xrantmp = xran[*,i]
     xtit = 'Observed Wavelength (!3' + STRING(197B) + '!X)'
     if i eq 0 then ytit='Normalized F$\down$$\lambda$' else ytit=''
     cgplot,wave,relfluxuse,xran=xrantmp,yran=yran,xstyle=1,ystyle=1,$
            backg='Black',axiscolor='White',color='White',position=pos[*,i],$
            xtit=xtit,ytit=ytit,noerase=i ne 0,/norm
     for j=0,n_elements(wavelines)-1 do begin
        if wavelines[j] gt xrantmp[0] AND wavelines[j] lt xrantmp[1] then begin
           cgoplot,[wavelines[j],wavelines[j]],yran,color='Green',linesty=2
           cgtext,wavelines[j]-1,0.1,wavenames[j],orient=90,align=0,color='Green'
        endif
     endfor
     for j=0,nabs-1 do cgoplot,wave,modabs[*,j],color='Cyan',thick=2
     cgoplot,wave,modflux,color='Red',thick=4
  endfor

  if ~keyword_set(init) then $
      img = cgsnapshot(filename=outfile,/jpeg,/nodialog,quality=100)
  
; Printing plot data
  if ~keyword_set(init) then begin

     output=directoryname+'/'+gal+'/'+gal+multiplet+'_fit_data'
     openw, lun, output+'.txt',/GET_LUN
     FOR M =0, N_ELEMENTS(wave)-1 DO $
        printf, lun, wave[M], modflux[M], continuum[M], flux[M], $
                relfluxuse[M], $
                FORMAT='(F10.4,2X,F7.4,2X,E13.6,2X,E13.6,2X,F7.4)'
     FREE_LUN,lun
  
;     openw,lun,output+'param.txt',/GET_LUN
;     printf,lun,xran[0],xran[1],yran[0],yran[1],w_doublet1,w_doublet2,$
;            FORMAT='(6(F10.2))'
;     FREE_LUN,lun
  
     openw,lun, output+'modabs'+'.txt',/GET_LUN
     abssize = size(modabs)
     abssize=abssize(1)
     printf, lun, nabs, abssize, FORMAT='(I4,2X,I4)'
     for i=0,nabs-1 do printf, lun, modabs[*,i], FORMAT='(F10.7)'
     free_lun, lun

  endif

end
