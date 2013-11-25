; docformat = 'rst'
;
;+
;
; Refine initial guess for stellar redshift by fitting CaII K.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Best-fit redshift.
;
; :Params:
;    zinit: in, required, type=double
;      Initial redshift guess.
;    lambda: in, required, type=dblarr(N)
;      Wavelengths of spectrum.
;    flux: in, required, type=dblarr(N)
;      Fluxes of spectrum.
;    flux: in, required, type=dblarr(N)
;      Flux errors of spectrum.
;
; :Keywords:
;    sigma: out, optional, type=double
;      Best-fit sigma.
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
;      2011mar, Joshua Fuchs, created
;      2011may09, DSNR, ported to generic function
;      2013nov25, DSNR, documented, renamed, added license and copyright 
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
function uhsf_refinestz,zinit,wave,flux,err,sigma=sigma

  carest = 3934.777d

; Determine Ca continuum
  w1 = (1d +zinit)*carest -20d
  w2 = (1d +zinit)*carest -12d
  w3 = (1d +zinit)*carest +12d
  w4 = (1d +zinit)*carest +20d
  caindex = where((wave GT w1 AND wave LT w2) OR $
                  (wave GT w3 AND wave LT w4),ct)
  if caindex[0] ne -1 then begin
     p = poly_fit(wave[caindex],flux[caindex],1,$
                  measure_errors=1d/err[caindex]^2d )

; Subtract continuum
     caab = where(wave GT w2 AND wave LT w3)
     cacontflux = poly(wave[caab],p)
     castart = flux[caab] - cacontflux
  
; Fit Gaussian
     cafit = mpfitpeak(wave[caab],castart,a,$
                       error=err[caab],/NEGATIVE,nterms=3)

     caredshift = (a[1] / carest) -1d
     sigma = a[2]

  endif else begin

     print,'IFSF_REFINESTZ: Not enough data to fit continuum, or'
     print,'CaII K line not present. Returning initial guess.'
     caredshift = zinit

  endelse

  return,caredshift

end
