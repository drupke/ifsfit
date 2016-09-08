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
;    A structure of normalized flux and error vectors.
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
;    pars: out, required, type=dblarr(M)
;      Polynomial coefficients of continuum normalization.
;
; :Keywords:
;    fitord: in, optional, type=double, default=2
;      Polynomial order to use in fitting continuum.
;    fitranlo: in, optional, type=dblarr(2), default=[5810\,5865]*(1+z)
;      Wavelength limits of region below Na D to use in fit.
;    fitranhi: in, optional, type=dblarr(2), default=[5905\,5960]*(1+z)
;      Wavelength limits of region above Na D to use in fit.
;    subtract: in, optional, type=byte
;      Subtract fitted continuum instead of dividing.
;    snavg_thresh: in, optional, type=double
;      Average flux/err in fit region below which spaxel is rejected.
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
;      2014may08, DSNR, tweaked for clarity
;      2014may27, DSNR, added output of un-normalized error and 
;                       subtraction-normalized flux; converted output to a 
;                       structure
;      2014jun10, DSNR, now outputs indices into original arrays, as well
;      2016may03, DSNR, added trigger to deal with very noisy data or data
;                       where continuum goes below 0
;      2016sep02, DSNR, fixed bug where fitting order was too small by 1
;    
; :Copyright:
;    Copyright (C) 2013--2016 David S. N. Rupke
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
function ifsf_normnad,wave,flux,err,z,pars,fitord=fitord,$
                      fitranlo=fitranlo,fitranhi=fitranhi,subtract=subtract,$
                      snavg_thresh=snavg_thresh


   if ~ keyword_set(fitord) then fitord=2
   if ~ keyword_set(fitranlo) then fitranlo = (1d +z)*[5810d,5865d]
   if ~ keyword_set(fitranhi) then fitranhi = (1d +z)*[5905d,5960d]

   ifit = where((wave ge fitranlo[0] AND wave le fitranlo[1]) OR $
                (wave ge fitranhi[0] AND wave le fitranhi[1]))
   igd = where(wave ge fitranlo[0] AND wave le fitranhi[1],ctgd)

;  Check to make sure data is OK
   snavg = mean(flux[ifit]/err[ifit])
   if ~ keyword_set(snavg_thresh) then snavg_thresh = 1d
   if snavg ge snavg_thresh then begin

      parinfo = replicate({value:0d},fitord+1)
      pars = mpfitfun('poly',wave[ifit],flux[ifit],err[ifit],parinfo=parinfo,/quiet)

      if keyword_set(subtract) then begin
         nflux = flux[igd] - poly(wave[igd],pars)
         nerr = err[igd]
      endif else begin
         nflux = flux[igd] / poly(wave[igd],pars)
         nerr = err[igd] / poly(wave[igd],pars)
      endelse
      
   endif else begin
  
      pars = dblarr(n_elements(fitord))
      nflux = dblarr(ctgd)
      nerr = dblarr(ctgd)
   
   endelse

   nwave = wave[igd]
   unflux = flux[igd]
   unerr = err[igd]
   out = {wave: nwave,$
          flux: unflux,$
          err: unerr,$
          nflux: nflux,$
          nerr: nerr,$
          ind: igd}
   
   return,out

end
