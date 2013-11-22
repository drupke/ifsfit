; docformat = 'rst'
;
;+
;
; Normalize region around Na D using a polynomial.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    type=dblarr(Nwave,2): Normalized flux and error vectors.
;
; :Params:
;    wave: in, required, type=dblarr(N)
;      Wavelengths of data to fit.
;    flux: in, required, type=dblarr(N)
;      Fluxes of data to fit.
;    err: in, required, type=dblarr(N)
;      Flux errors of data to fit.
;    z: in, required, type=double
;      Redshift to use in computing default fit ranges.
;
; :Keywords:
;    fitord: in, optional, type=double, default=2
;      Polynomial order to use in fitting continuum.
;    fitordlo: in, optional, type=dblarr(2), default=[5810\,5865]*(1+z)
;      Wavelength limits of region below Na D to use in fit.
;    fitordhi: in, optional, type=dblarr(2), default=[5905\,5960]*(1+z)
;      Wavelength limits of region above Na D to use in fit.
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
;      2013nov22, DSNR, created (copied normalization bits from old
;                       routine 'printnadspec')
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
function ifsf_normnadspec,wave,flux,err,z,fitord=fitord,$
                     fitranlo=fitranlo,fitranhi=fitranhi

  if ~ keyword_set(fitord) then fitord=2
  if ~ keyword_set(fitranlo) then fitranlo = (1d +z)*[5810d,5865d]
  if ~ keyword_set(fitranhi) then fitranhi = (1d +z)*[5905d,5960d]

  ifit = where((wave ge fitranlo[0] AND wave le fitranlo[1]) OR $
               (wave ge fitranhi[0] AND wave le fitranhi[1]))

  parinfo = replicate({value:0d},fitord)
  pars = mpfitfun('poly',wave[ifit],flux[ifit],err[ifit],parinfo=parinfo,/quiet)

  nflux = flux / poly(wave,pars)
  nerr = err / poly(wave,pars)

  return,[[nflux],[nerr]]

end
