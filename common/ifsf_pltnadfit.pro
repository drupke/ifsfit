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
;    datfile: in, required, type=
;      Full path and name of file containing data (with emission lines
;      subtracted), wavelength, and error in three columns.
;    parfile: in, required, type=
;      Parameter file output by Na D fitting routine.
;    outfile: in, required, type=string
;      Full path and name of output plot.
;    z: in, required, type=double
;      Redshift.
;
; :Keywords:
;    xran: in, optional, type=dblarr(2), default=[lower fit limit-20A\,upper fit limit-20A]
;      Wavelength range of output plot
;    yran: in, optional, type=dblarr(2), default=[0\,1.5]
;      Normalized flux range of output plot 
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
;      2013nov21, DSNR, documented, renamed, added license and copyright 
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
pro ifsf_pltnadfit,datfile,parfile,outfile,z,xran=xran,$
                   yran=yran

; Default wavelengths
  c = 299792.458d
  nad1_rest = 5895.92d
  nad2_rest = 5889.95d
  he_rest = 5875.661d

; Observed-frame wavelengths
  nad1 = nad1_rest*(1d +z)
  nad2 = nad2_rest*(1d +z)
  he = he_rest*(1d +z)

  spec = ifsf_cmpnadfull(datfile,parfile,z)

; PLOT
  
  set_plot,'Z'
  device,decomposed=0,set_resolution=[1280,960],set_pixel_depth=24
  !P.charsize=1
  !P.charthick=1

  if ~ keyword_set(xran) then xran = [opars.llo-20d,opars.lhi+20d]
  if ~ keyword_set(yran) then yran = [0d,1.5d]

  cgplot,spec[*,0],spec[*,1],xran=xran,yran=yran,/xsty,/ysty,$
         backg='Black',color='White'
  cgoplot,spec[*,0],spec[*,4],color='Blue',thick=4
  cgoplot,spec[*,0],spec[*,3],color='Red',thick=4
  cgoplot,spec[*,0],spec[*,2],color='Blue'
  cgoplot,[nad1,nad1],yran,color='Green',linesty=2
  cgoplot,[nad2,nad2],yran,color='Green',linesty=2
  cgoplot,[he,he],yran,color='Green',linesty=2

  img = cgsnapshot(filename=outfile,/jpeg,/nodialog,quality=100)

end
