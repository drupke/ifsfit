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
;    instr: in, required, type=structure
;      Results of fit.
;    outfile: in, required, type=string
;      Full path and name of output plot.
;    wavenorm: in, required, type=dblarr(N)
;      Wavelength array of normalized data.
;    fluxnorm: in, required, type=dblarr(N)
;      Flux array of normalized data.
;    parnorm: in, required, type=dblarr(M)
;      Polynomial coefficients of continuum normalization.
;
; :Keywords:
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
pro ifsf_pltnaddat,instr,outfile,wavenorm,fluxnorm,parnorm

; Plot defaults
  xran_rest = [5793,5993]
  nad1_rest = 5895.92d
  nad2_rest = 5889.95d
  he_rest   = 5875.661d
  cfit1ran_rest = [5800d,5865d]
  cfit2ran_rest = [5905d,5980d]

  cleanplot,/silent
  set_plot,'Z'
  device,decomposed=0,set_resolution=[1280,960],set_pixel_depth=24
  !P.charsize=1
  !P.charthick=1

  wave = instr.wave
  spectot = instr.spec
  specstars = instr.spec - instr.specfit
  speclines = instr.spec_nocnt
  modtot = instr.specfit + (instr.spec - instr.spec_nocnt)
  modstars = instr.spec - instr.spec_nocnt
  modresid = modstars
  specresid = specstars
  
  norm = max(modstars)
  spectot /= norm
  specstars /= norm
  speclines /= norm
  specresid /= norm
  modtot /= norm
  modstars /= norm
  modresid /= norm

  cfit1ran = (1d + instr.z)*cfit1ran_rest
  cfit2ran = (1d + instr.z)*cfit2ran_rest
  xran1 = (1d + instr.z)*xran_rest
  nad1 = (1d + instr.z)*nad1_rest
  nad2 = (1d + instr.z)*nad2_rest
  he = (1d + instr.z)*he_rest
  i1 = where(wave gt xran1[0] AND wave lt xran1[1],ct1)
  
  maxthresh=0.2
  ntop = 20
  ntop++
  nbottom = 20
  nbottom--

  ydat = specresid
  ymod = modresid
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
  yran[0] = 0d

  cgplot,wave,specresid,xran=xran1,yran=yran,/xsty,/ysty,$
         backg='Black',color='White',layout=[1,2,1]
  cgoplot,wave,poly(wave,parnorm)/norm,color='Red'
  cgoplot,[nad1,nad1],yran,color='Green',linesty=2
  cgoplot,[nad2,nad2],yran,color='Green',linesty=2
  cgoplot,[he,he],yran,color='Green',linesty=2
  ytmp = yran[0]+(yran[1]-yran[0])*0.05
  cgoplot,cfit1ran,[ytmp,ytmp],color='Blue',linesty=2
  cgoplot,cfit2ran,[ytmp,ytmp],color='Blue',linesty=2

  yran = [0,1.3]
  cgplot,wavenorm,fluxnorm,xran=xran1,yran=yran,/xsty,/ysty,$
         layout=[1,2,2]
  cgoplot,xran1,[1,1],color='Red'
  cgoplot,[nad1,nad1],yran,color='Green',linesty=2
  cgoplot,[nad2,nad2],yran,color='Green',linesty=2
  cgoplot,[he,he],yran,color='Green',linesty=2

  img = cgsnapshot(filename=outfile,/jpeg,/nodialog,quality=100)
 
end
