; docformat = 'rst'
;
;+
;
; Plot emission-line subtracted continuum and normalized continuum
; around Na D.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    None.
;
; :Params:
;    normdat: in, required, type=dblarr(Npts,4)
;      Results of normalization (IFSF_NORMNAD). Includes wavelengths, 
;      unnormalized continuum (minus emission lines), normalized continuum,
;      and normalized error spectra.
;    fitpars: in, required, type=dblarr(M)
;      Polynomial coefficients of continuum normalization.
;    z: in, required, type=double
;      Redshift to use in computing default fit ranges.
;    outfile: in, required, type=string
;      Full path and name of output plot.
;
; :Keywords:
;    pltran: in, optional, type=double
;      Plot range.
;    fitranlo: in, optional, type=dblarr(2), default=[5810\,5865]*(1+z)
;      Wavelength limits of region below Na D used in fit.
;    fitranhi: in, optional, type=dblarr(2), default=[5905\,5960]*(1+z)
;      Wavelength limits of region above Na D used in fit.
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
;      2010jul19, DSNR, created
;      2013nov21, DSNR, documented, renamed, added license and copyright
;      2014may08, DSNR, tweaked for clarity
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
pro ifsf_pltnaddat,normdat,fitpars,z,outfile,pltran=pltran,$
                   fitranlo=fitranlo,fitranhi=fitranhi

   if ~ keyword_set(pltran) then pltran = (1d +z)*[5810d,5960d]
   if ~ keyword_set(fitranlo) then fitranlo = (1d +z)*[5810d,5865d]
   if ~ keyword_set(fitranhi) then fitranhi = (1d +z)*[5905d,5960d]

   linelist = ifsf_linelist(['NaD1','NaD2','HeI5876'])
   nad1 = (1d + z)*linelist['NaD1']
   nad2 = (1d + z)*linelist['NaD2']
   hei = (1d + z)*linelist['HeI5876']

   wave = normdat[*,0]
   flux = normdat[*,1]
   nflux = normdat[*,2]
   iplt = where(wave ge pltran[0] AND wave le pltran[1],ctplt)

   set_plot,'Z'
   device,decomposed=0,set_resolution=[1280,960],set_pixel_depth=24
   !P.charsize=1
   !P.charthick=1
   erase
   
   yran = [min(flux[iplt]),max(flux[iplt])]
   cgplot,wave[iplt],flux[iplt],xran=pltran,yran=yran,/xsty,/ysty,$
          color='White',axisc='White',layout=[1,2,1]
   cgoplot,wave[iplt],poly(wave[iplt],fitpars),color='Red'
   cgoplot,[nad1,nad1],yran,color='Green',linesty=2
   cgoplot,[nad2,nad2],yran,color='Green',linesty=2
   cgoplot,[hei,hei],yran,color='Green',linesty=2
   ytmp = yran[0]+(yran[1]-yran[0])*0.05
   cgoplot,fitranlo,[ytmp,ytmp],color='Cyan',linesty=2
   cgoplot,fitranhi,[ytmp,ytmp],color='Cyan',linesty=2

   yran = [0,1.3]
   cgplot,wave[iplt],nflux[iplt],xran=pltran,yran=yran,/xsty,/ysty,$
          color='White',axisc='White',layout=[1,2,2]
   cgoplot,pltran,[1,1],color='Red'
   cgoplot,[nad1,nad1],yran,color='Green',linesty=2
   cgoplot,[nad2,nad2],yran,color='Green',linesty=2
   cgoplot,[hei,hei],yran,color='Green',linesty=2

   img = cgsnapshot(filename=outfile,/jpeg,/nodialog,quality=100)
 
end
