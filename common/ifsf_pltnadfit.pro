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
;    cmodcomp: in, optional, type=string, default='Cyan'
;    cvlab: in, optional, type=string, default='Cyan'
;    czstarline: in, optional, type=string, default='Orange'
;    czstarval: in, optional, type=string, default='Orange'
;    czsysline: in, optional, type=string, default='Green'
;    czsysrval: in, optional, type=string, default='Green'
;      Color for various things ...
;    ps: in, optional, type=byte
;      Select .eps output.
;    specres: in, required, type=double
;    upsample: in, optional, type=integer, def=0
;      Factor by which to "upsample" the spectra (opposite of bin!) for the
;      model computation.
;    xran: in, optional, type=dblarr(2), default=all
;      Wavlength range of output plot.
;    yran: in, optional, type=dblarr(2), default=[0\,1.5]
;      Normalized flux range of output plot.
;    zstar: in, optional, type=double
;      Stellar redshift.
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
;      2020jun25, DSNR, now also plot v(NaD em) and only plot zsys if different 
;                       from zstar; option to output .eps file; color options
;      2021oct21, DSNR, made upsample parameter
;      2022jul29, DSNR, changed a bunch of labelling; possibly not backwards compatible ...
;
; :Copyright:
;    Copyright (C) 2013--2022 David S. N. Rupke
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
pro ifsf_pltnadfit,wave,flux,err,param,outfile,zsys,xran=xran,yran=yran,$
                   ps=ps,specres=specres,zstar=zstar,$
                   upsample=upsample,vacuum=vacuum,$
                   ; colors
                   cvlab=cvlab,cmodcomp=cmodcomp,czstarline=czstarline,$
                   czstarval=czstarval,czsysval=czsysval,czsysline=czsysline
                   
  if ~keyword_set(specres) then specres=0b
  if ~keyword_set(upsample) then upsample=0b
  if ~keyword_set(vacuum) then vacuum=0b
  if keyword_set(ps) then dops=1 else dops=0
  if ~keyword_set(cvlab) then cvlab='Cyan'
  if ~keyword_set(cmodcomp) then cmodcomp='Cyan'
  if ~keyword_set(czstarline) then czstarline='Orange'
  if ~keyword_set(czstarval) then czstarval='Orange'
  if ~keyword_set(czsysline) then czsysline='Green'
  if ~keyword_set(czsysval) then czsysval='Green'

; Default wavelengths
  c = 299792.458d
  linelist = ifsf_linelist(['NaD1','NaD2','HeI5876'],/quiet,vacuum=vacuum)

; Observed-frame wavelengths
  w_nad1 = linelist['NaD1']*(1d +zsys)
  w_nad2 = linelist['NaD2']*(1d +zsys)
  w_hei = linelist['HeI5876']*(1d +zsys)
  if keyword_set(zstar) then begin
     wst_nad1 = linelist['NaD1']*(1d +zstar)
     wst_nad2 = linelist['NaD2']*(1d +zstar)
  endif
  
  modhei=1
  modnadabs=1
  modnadem=1
  modflux = ifsf_nadfcn(wave,param,modhei=modhei,modnadabs=modnadabs,$
                        modnadem=modnadem,specres=specres,$
                        wavnadabs=wavnadabs,wavhei=wavhei,wavnadem=wavnadem,$
                        upsample=upsample)
  size_hei = size(modhei)
  if size_hei[0] eq 1 then nhei = 1 $
  else if size_hei[0] eq 2 then nhei = fix(size_hei[2]) $
  else nhei=0l
  size_nadabs = size(modnadabs)
  if size_nadabs[0] eq 1 then nnadabs = 1 $
  else if size_nadabs[0] eq 2 then nnadabs = fix(size_nadabs[2]) $
  else nnadabs=0l
  size_nadem = size(modnadem)
  if size_nadem[0] eq 1 then nnadem = 1 $
  else if size_nadem[0] eq 2 then nnadem = fix(size_nadem[2]) $
  else nnadem=0l

; PLOT
  
  if dops then begin
     cgps_open,filename=outfile+'.eps',/encapsulated,xsize=10,ysize=7.5,$
               bits_per_pixel=8,/color,/inches,/nomatch
     !P.charsize=1
     !P.charthick=2
     !P.thick=2
     cbackg = 'White'
     caxis='Black'
     cdef='Black'
   endif else begin
     set_plot,'Z'
     device,decomposed=0,set_resolution=[1280,960],set_pixel_depth=24
     !P.charsize=1.5
     !P.charthick=1
     cgerase,'Black'
     cbackg = 'Black'
     caxis='White'
     cdef='White'
   endelse

  if ~ keyword_set(xran) then xran = [min(wave),max(wave)]
  if ~ keyword_set(yran) then yran = [0d,1.5d]

  cgplot,wave,flux,xran=xran,yran=yran,/xsty,/ysty,$
         backg=cbackg,axiscolor=caxis,color=cdef,$
         xtit='Wavelength ($\Angstrom$)',ytit='Normalized F$\down$$\lambda$'
  cgoplot,wave,err,color=cdef,linesty=2
  if keyword_set(zstar) then begin
     cgoplot,[wst_nad1,wst_nad1],yran,color=czstarline,linesty=2
     cgoplot,[wst_nad2,wst_nad2],yran,color=czstarline,linesty=2
  endif
  if zsys ne zstar then begin
     cgoplot,[w_nad1,w_nad1],yran,color=czsysline,linesty=2
     cgoplot,[w_nad2,w_nad2],yran,color=czsysline,linesty=2
     cgoplot,[w_hei,w_hei],yran,color=czsysline,linesty=2
  endif else begin
     cgoplot,[w_hei,w_hei],yran,color=czstarline,linesty=2
  endelse
  for i=0,nhei-1 do cgoplot,wave,modhei[*,i],color=cmodcomp,thick=2
  for i=0,nnadabs-1 do cgoplot,wave,modnadabs[*,i],color=cmodcomp,thick=2
  for i=0,nnadem-1 do cgoplot,wave,modnadem[*,i],color=cmodcomp,thick=2
  cgoplot,wave,modflux,color='Red',thick=4
  xtit = xran[1] - (xran[1]-xran[0])*0.05
  ytit = yran[1] - (yran[1]-yran[0])*0.05
  if keyword_set(zstar) then $
    cgtext,xtit,ytit,'z$\downstar$ = '+string(zstar,format='(D0.5)'),$
           color=czstarval,align=1
  if zsys ne zstar then begin
     ytit -= (yran[1]-yran[0])*0.05
     cgtext,xtit,ytit,'z$\downem$ = '+string(zsys,format='(D0.5)'),$
            color=czsysval,align=1
  endif
  for i=0,nnadabs-1 do begin
     ytit -= (yran[1]-yran[0])*0.05
     if keyword_set(zstar) then begin
        zdiff = wavnadabs[i]/wst_nad1 - 1d
        vel = c * ((zdiff+1d)^2d - 1d) / ((zdiff+1d)^2d + 1d)
        cgtext,xtit,ytit,$
           '$\Delta$v'+string('(NaDabs',i+1,'-star)=',vel,' km/s',$
           format='(A0,I0,A0,I0,A0)'),color=cvlab,align=1
     endif else begin
        zdiff = wavnadabs[i]/w_nad1 - 1d
        vel = c * ((zdiff+1d)^2d - 1d) / ((zdiff+1d)^2d + 1d)
        cgtext,xtit,ytit,$
           '$\Delta$v'+string('(NaDabs',i+1,'-em)=',vel,' km/s',$
           format='(A0,I0,A0,I0,A0)'),color=cvlab,align=1
     endelse
;     cgtext,xtit,ytit,$
;            string('NaD abs., comp. ',i+1,' v = ',vel,' km/s',$
;                   format='(A0,I0,A0,I0,A0)'),$
;            color=cvlab,align=1
  endfor
  for i=0,nnadem-1 do begin
    ytit -= (yran[1]-yran[0])*0.05
    zdiff = wavnadem[i]/w_nad1 - 1d
    vel = c * ((zdiff+1d)^2d - 1d) / ((zdiff+1d)^2d + 1d)
;    cgtext,xtit,ytit,$
;           string('NaD em., comp. ',i+1,' v = ',vel,' km/s',$
;           format='(A0,I0,A0,I0,A0)'),$
;           color=cvlab,align=1
   cgtext,xtit,ytit,$
      '$\Delta$v'+string('(NaDem',i+1,'-em)=',vel,' km/s',$
      format='(A0,I0,A0,I0,A0)'),color=cvlab,align=1

  endfor
  for i=0,nhei-1 do begin
     ytit -= (yran[1]-yran[0])*0.05
     zdiff = wavhei[i]/w_hei - 1d
     vel = c * ((zdiff+1d)^2d - 1d) / ((zdiff+1d)^2d + 1d)
     if zsys ne zstar then begin
        cgtext,xtit,ytit,$
           '$\Delta$v'+string('(HeI',i+1,'-em)=',vel,' km/s',$
           format='(A0,I0,A0,I0,A0)'),color=cvlab,align=1
     endif else begin
        cgtext,xtit,ytit,$
           '$\Delta$v'+string('(HeI',i+1,'-star)=',vel,' km/s',$
           format='(A0,I0,A0,I0,A0)'),color=cvlab,align=1
     endelse
  endfor

  if dops then cgps_close $
  else img = cgsnapshot(filename=outfile,/jpeg,/nodialog,quality=100)

end
