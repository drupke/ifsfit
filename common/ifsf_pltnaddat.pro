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
;    autoindices: in, optional, type=lonarr(4)
;      Locations of boundaries for empirically finding lines.
;    emwid: in, optional, type=double, default=15d
;      Wavelength range over which flux was integrated to estimate upper limit for 
;      emission line equivalent width and flux.
;    iabsoff: in, optional, type=long, default=4l
;      Index offset from absorption line for calculating emission line upper limit.
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
;      2014jul29, DSNR, added AUTOINDICES, IABSOFF, EMWID keywords, and legend
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
                   fitranlo=fitranlo,fitranhi=fitranhi,$
                   autoindices=autoindices,emwid=emwid,iabsoff=iabsoff,$
                   yran=yran,vacuum=vacuum

   if ~ keyword_set(iabsoff) then iabsoff = 4l
   if ~ keyword_set(emwid) then emwid = 20d 
   if ~ keyword_set(vacuum) then vacuum = 0b

   if ~ keyword_set(pltran) then pltran = (1d +z)*[5810d,5960d]
   if ~ keyword_set(fitranlo) then fitranlo = (1d +z)*[5810d,5865d]
   if ~ keyword_set(fitranhi) then fitranhi = (1d +z)*[5905d,5960d]

   linelist = ifsf_linelist(['NaD1','NaD2','HeI5876'],vacuum=vacuum,/quiet)
   nad1 = (1d + z)*linelist['NaD1']
   nad2 = (1d + z)*linelist['NaD2']
   hei = (1d + z)*linelist['HeI5876']

   wave = normdat.wave
   nflux = normdat.nflux
   flux = normdat.flux
   iplt = where(wave ge pltran[0] AND wave le pltran[1],ctplt)

   set_plot,'Z'
   cgdisplay,1280,960
   cgerase,'Black'
   !P.charsize=1
   !P.charthick=1
   !P.thick=2

   autocol='Red'
   cntcol='Yellow'
   lincol='Cyan'
   
   yran = [min(flux[iplt]),max(flux[iplt])]
   if yran[0] lt 0 then yran[0] = 0
   cgplot,wave[iplt],flux[iplt],xran=pltran,yran=yran,/xsty,/ysty,$
          color='White',axisc='White',layout=[1,2,1],backg='Black',$
          ytit='Actual Flux',xtit='Observed Wavelength ($\Angstrom$)'
   cgoplot,pltran,[0,0],color='White'
   cgoplot,wave[iplt],poly(wave[iplt],fitpars),color=cntcol
   cgoplot,[nad1,nad1],yran,color=lincol,linesty=2
   cgoplot,[nad2,nad2],yran,color=lincol,linesty=2
   cgoplot,[hei,hei],yran,color=lincol,linesty=2
   ytmp = yran[0]+(yran[1]-yran[0])*0.05
   cgoplot,fitranlo,[ytmp,ytmp],color=cntcol,linesty=2
   cgoplot,fitranhi,[ytmp,ytmp],color=cntcol,linesty=2
   if keyword_set(autoindices) then begin
      doabs=0
;     Absorption line region
      if autoindices[0] ne -1l AND autoindices[1] ne -1l then begin
         cgoplot,[wave[autoindices[0]],wave[autoindices[0]]],yran,$
                 color=autocol,linesty=1
         cgoplot,[wave[autoindices[1]],wave[autoindices[1]]],yran,$
                 color=autocol,linesty=1
         doabs=1
      endif
;     Emission line region, case of detection
      if autoindices[2] ne -1l AND autoindices[3] ne -1l then begin
         cgoplot,[wave[autoindices[2]],wave[autoindices[2]]],yran,$
                 color=autocol,linesty=1
         cgoplot,[wave[autoindices[3]],wave[autoindices[3]]],yran,$
                 color=autocol,linesty=1
;     Emission line region, case of upper limit
      endif else if doabs then begin
         cgoplot,[wave[autoindices[1]+iabsoff],wave[autoindices[1]+iabsoff]],yran,$
                 color=autocol,linesty=1
         iup = value_locate(wave,wave[autoindices[1]+iabsoff]+$
                            emwid)
         cgoplot,[wave[iup],wave[iup]],yran,$
                 color=autocol,linesty=1
      endif
   endif
   al_legend,['Normalization region','Galaxy rest frame','Abs/Em Ranges'],$
             col=[cntcol,lincol,autocol],linesty=[2,2,1],/left,/top


   yran = [0,1.3]
   cgplot,wave[iplt],nflux[iplt],xran=pltran,yran=yran,/xsty,/ysty,$
          color='White',axisc='White',layout=[1,2,2],ytit='Normalized Flux'
   cgoplot,pltran,[1,1],color=cntcol
   cgoplot,[nad1,nad1],yran,color=lincol,linesty=2
   cgoplot,[nad2,nad2],yran,color=lincol,linesty=2
   cgoplot,[hei,hei],yran,color=lincol,linesty=2
   if keyword_set(autoindices) then begin
      doabs=0
;     Absorption line region
      if autoindices[0] ne -1l AND autoindices[1] ne -1l then begin
         cgoplot,[wave[autoindices[0]],wave[autoindices[0]]],yran,$
                 color=autocol,linesty=1
         cgoplot,[wave[autoindices[1]],wave[autoindices[1]]],yran,$
                 color=autocol,linesty=1
         doabs=1
      endif
;     Emission line region, case of detection
      if autoindices[2] ne -1l AND autoindices[3] ne -1l then begin
         cgoplot,[wave[autoindices[2]],wave[autoindices[2]]],yran,$
                 color=autocol,linesty=1
         cgoplot,[wave[autoindices[3]],wave[autoindices[3]]],yran,$
                 color=autocol,linesty=1
;     Emission line region, case of upper limit
      endif else if doabs then begin
         cgoplot,[wave[autoindices[1]+iabsoff],wave[autoindices[1]+iabsoff]],yran,$
                 color=autocol,linesty=1
         iup = value_locate(wave,wave[autoindices[1]+iabsoff]+$
                            emwid)
         cgoplot,[wave[iup],wave[iup]],yran,$
                 color=autocol,linesty=1
      endif
   endif
   al_legend,['Continuum fit','Galaxy rest frame','Abs/Em Ranges'],$
             col=[cntcol,lincol,autocol],linesty=[0,2,1],/left,/bottom

   img = cgsnapshot(filename=outfile,/jpeg,/nodialog,quality=100)
 
end
