; docformat = 'rst'
;
;+
;
; Plot Na D data and fit.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    None.
;
; :Params:
;    wave: in, required, type=dblarr(N)
;      Wavelengths.
;    flux: in, required, type=dblarr(N)
;      Fluxes.
;    param: in, required, type=dblarr
;      Parameters output by Na D fitting routine.
;    outfile: in, required, type=string
;      Full path and name of output plot.
;    zsys: in, required, type=double
;      Systemic redshift.
;
; :Keywords:
;    xran: in, optional, type=dblarr(2), default=all
;      Wavlength range of output plot.
;    yran: in, optional, type=dblarr(2), default=[0\,1.5]
;      Normalized flux range of output plot.
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
;
; :Copyright:
;    Copyright (C) 2013-2014 David S. N. Rupke
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
pro ifsf_pltuvabsfitovi,wave,flux,param,outfile,zsys,xran=xran,yran=yran

; Default wavelengths
  c = 299792.458d
  linelist = ifsf_linelist(['[OVI1]1032','[OVI2]1038','[CII]1347'])

; Observed-frame wavelengths
  w_nad1 = linelist['[OVI1]1032']*(1d +zsys)
  w_nad2 = linelist['[OVI2]1038']*(1d +zsys)
  w_hei = linelist['[CII]1347']*(1d +zsys)

  modhei=1
  modnadabs=1
  modnadem=1
  modflux = ifsf_uvabsfcn(wave,param,modhei=modhei,modnadabs=modnadabs,$
                        modnadem=modnadem)
;  size_hei = size(modhei)
;  if size_hei[0] eq 1 then nhei = 1 $
;  else if size_hei[0] eq 2 then nhei = fix(size_hei[2]) $
;  else nhei=0l
  size_nadabs = size(modnadabs)
  if size_nadabs[0] eq 1 then nnadabs = 1 $
  else if size_nadabs[0] eq 2 then nnadabs = fix(size_nadabs[2]) $
  else nnadabs=0l
  size_nadem = size(modnadem)
  if size_nadem[0] eq 1 then nnadem = 1 $
  else if size_nadem[0] eq 2 then nnadem = fix(size_nadem[2]) $
  else nnadem=0l

; PLOT
  
  set_plot,'Z'
  device,decomposed=0,set_resolution=[1280,960],set_pixel_depth=24
  !P.charsize=1.5
  !P.charthick=1
  erase

;  if ~ keyword_set(xran) then xran = [min(wave),max(wave)]
;  if ~ keyword_set(yran) then yran = [0d,1.5d]
  yran = [0d,1.5d]
  xran = [min(wave),max(wave)]
  
  cgplot,wave,flux,xran=xran,yran=yran,xstyle=1,ystyle=1,$
         backg='Black',axiscolor='White',color='White',$
         xtit='Wavelength ($\Angstrom$)',ytit='Normalized F$\down$$\lambda$'
  cgoplot,[w_nad1,w_nad1],yran,color='Green',linesty=2
  cgoplot,[w_nad2,w_nad2],yran,color='Green',linesty=2
  cgoplot,[w_hei,w_hei],yran,color='Green',linesty=2
;  for i=0,nhei-1 do cgoplot,wave,modhei[*,i],color='Cyan',thick=2
  for i=0,nnadabs-1 do cgoplot,wave,modnadabs[*,i],color='Cyan',thick=2
  for i=0,nnadem-1 do cgoplot,wave,modnadem[*,i],color='Cyan',thick=2
  cgoplot,wave,modflux,color='Red',thick=4

  img = cgsnapshot(filename=outfile,/jpeg,/nodialog,quality=100)

end
