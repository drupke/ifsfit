; docformat = 'rst'
;
;+
;
; Plot doublet data and fit, and output plot information into three files.
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
;    linelist: in, required, type=hash
;      Line list and labels output by IFSF_LINELIST.
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
;      2020aug03, DSNR, added linelist parameter, rather than calling again
;                       from this routine
;
; :Copyright:
;    Copyright (C) 2013--2020 David S. N. Rupke, Anthony To
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
pro ifsf_pltdoublet,gal,wave,relativeflux,continuum,flux,param,doublet,$
                    directoryname,outfile,zsys,linelist,$
                    xran=xran,yran=yran,init=init,specres=specres,$
                    upsample=upsample

  if ~  keyword_set(specres) then specres=0b
  if ~  keyword_set(upsample) then upsample=0b

; Default wavelengths
  c = 299792.458d
  IF (doublet eq 'NV') THEN BEGIN
; Observed-frame wavelengths
    w_doublet1 = linelist['NV1238']*(1d +zsys)
    w_doublet2 = linelist['NV1242']*(1d +zsys)
  END
  IF (doublet eq 'OVI') THEN BEGIN
; Observed-frame wavelengths
    w_doublet1 = linelist['OVI1031']*(1d +zsys)
    w_doublet2 = linelist['OVI1037']*(1d +zsys)
  END
  IF (doublet eq 'PV') THEN BEGIN
     ; Observed-frame wavelengths
     w_doublet1 = linelist['PV1117']*(1d +zsys)
     w_doublet2 = linelist['PV1128']*(1d +zsys)
  END
  IF (doublet eq 'MgII') THEN BEGIN
     ; Observed-frame wavelengths
     w_doublet1 = linelist['MgII2796']*(1d +zsys)
     w_doublet2 = linelist['MgII2803']*(1d +zsys)
  END
  IF (doublet eq 'CaII') THEN BEGIN
    ; Observed-frame wavelengths
    w_doublet1 = linelist['CaIIK']*(1d +zsys)
    w_doublet2 = linelist['CaIIH']*(1d +zsys)
  END
  IF (doublet eq 'FeIIUV1') THEN BEGIN
     ; Observed-frame wavelengths
     w_doublet1 = linelist['FeII2585']*(1d +zsys)
     w_doublet2 = linelist['FeII2599']*(1d +zsys)
  END
  IF (doublet eq 'FeIIUV2') THEN BEGIN
     ; Observed-frame wavelengths
     w_doublet1 = linelist['FeII2373']*(1d +zsys)
     w_doublet2 = linelist['FeII2382']*(1d +zsys)
  END
  modabs=1
  modem=1
  modflux = ifsf_doubletfcn(wave,param,doubletname=doublet,$
                            modabs=modabs,$
                            modem=modem,specres=specres,upsample=upsample)
  size_doubletabs = size(modabs)
  if size_doubletabs[0] eq 1 then nabs = 1 $
  else if size_doubletabs[0] eq 2 then nabs = fix(size_doubletabs[2]) $
  else nabs=0l
  size_doubletem = size(modem)
  if size_doubletem[0] eq 1 then nem = 1 $
  else if size_doubletem[0] eq 2 then nem = fix(size_doubletem[2]) $
  else nem=0l

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

  if ~ keyword_set(xran) then xran = [min(wave),max(wave)]
  if ~ keyword_set(yran) then yran = [0d,1.5d]
  
  cgplot,wave,relativeflux,xran=xran,yran=yran,xstyle=1,ystyle=1,$
         backg='Black',axiscolor='White',color='White',$
         xtit='Wavelength ($\Angstrom$)',ytit='Normalized F$\down$$\lambda$'
  cgoplot,[w_doublet1,w_doublet1],yran,color='Green',linesty=2
  cgoplot,[w_doublet2,w_doublet2],yran,color='Green',linesty=2
;  for i=0,nhei-1 do cgoplot,wave,modhei[*,i],color='Cyan',thick=2
  for i=0,nabs-1 do cgoplot,wave,modabs[*,i],color='Cyan',thick=2
  for i=0,nem-1 do cgoplot,wave,modem[*,i],color='Cyan',thick=2
  cgoplot,wave,modflux,color='Red',thick=4

  if ~keyword_set(init) then $
      img = cgsnapshot(filename=outfile,/jpeg,/nodialog,quality=100)
  
; Printing plot data
  if ~keyword_set(init) then begin

     output=directoryname+gal+doublet+'_fit_data'
     openw, lun, output+'.txt',/GET_LUN
     FOR M =0, N_ELEMENTS(wave)-1 DO $
        printf, lun, wave[M], modflux[M], continuum[M], flux[M], $
                relativeflux[M], $
                FORMAT='(F10.4,2X,F7.4,2X,E13.6,2X,E13.6,2X,F10.4)'
     FREE_LUN,lun
  
     openw,lun,output+'param.txt',/GET_LUN
     printf,lun,xran[0],xran[1],yran[0],yran[1],w_doublet1,w_doublet2,$
            FORMAT='(6(F10.2))'
     FREE_LUN,lun
  
     openw,lun, output+'modabs'+'.txt',/GET_LUN
     abssize = size(modabs)
     abssize=abssize(1)
     printf, lun, nabs, abssize, FORMAT='(I4,2X,I4)'
     for i=0,nabs-1 do printf, lun, modabs[*,i], FORMAT='(F10.7)'
     free_lun, lun

  endif

end
