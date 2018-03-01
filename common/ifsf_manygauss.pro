; docformat = 'rst'
;
;+
; This function is input to MPFIT, which uses it to compute the
; emission line spectrum of multiple Gaussian emission lines
; simultaneously. This routine assumes constant dispersion (in A/pix),
; and uses this fact to optimize the computation by working with
; smaller sub-arrays. Use MANYGAUSS_SLOW (which works with large
; arrays) for data without constant dispersion.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    An N-element array containing an emission-line spectrum.
;
; :Params:
;    wave: in, required, type=dblarr(N)
;      Wavelengths
;    param: in, required, type=dblarr
;      Best-fit parameter array output by MPFIT.
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
;      2009, DSNR, copied from manygauss_slow.pro and rewritten
;      2013sep, DSNR, switch sigma from wavelength to velocity space
;      2013nov13, DSNR, documented, renamed, added license and copyright
;      2014apr10, DSNR, fixed cases of floating underflow
;      2016sep26, DSNR, switched to deconvolving resolution in situ by adding
;                       resolution in quadrature in wavelength space to each
;                       velocity component
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
function ifsf_manygauss, wave, param, specresarr=specresarr

  c = 299792.458d

  ppoff = param[0]
  nwave = n_elements(wave)
  nline = (N_ELEMENTS(param)-ppoff)/3

  find = ppoff+indgen(nline)*3
  wind = find + 1
  sind = find + 2

  dispersion = wave[1] - wave[0]
  if keyword_set(specresarr) then begin
     wsr = value_locate(specresarr[*,0],param[wind])
     srsigslam = param[wind]/specresarr[wsr,1]/2.35d
  endif else begin
     srsigslam = dblarr(nline)+param[2]
  endelse
; resolution in wavelength space [sigma] assumed to be in third element of PARAM
  sigs = sqrt((param[sind]/c * param[wind])^2d + srsigslam^2d)
  maxsig = max(sigs)

  nsubwave = round(10d * maxsig / dispersion)
  halfnsubwave = round(nsubwave / 2)
  nsubwave = halfnsubwave*2+1
  indsubwaves = rebin(transpose(fix(indgen(nsubwave)-halfnsubwave)),$
                      nline,nsubwave)

  fluxes = param[find]
  refwaves = param[wind]
  indrefwaves_real = (refwaves - wave[0]) / dispersion
  indrefwaves = fix(indrefwaves_real)
  indrefwaves_frac = indrefwaves_real - double(indrefwaves)
  dwaves = (indsubwaves - rebin(indrefwaves_frac,nline,nsubwave))*$
           dispersion
  indsubwaves += rebin(indrefwaves,nline,nsubwave)
  
  yvals = dblarr(nwave)
  for i=0,nline-1 do begin
     gind = where(indsubwaves[i,*] ge 0 AND indsubwaves[i,*] le nwave-1,$
                  count)
;    The "mask" parameter eliminates floating underflow by removing very large
;    negative exponents. See http://www.idlcoyote.com/math_tips/underflow.html
;    for more details.
     if count gt 0 then begin
        exparg = -(dwaves[i,gind]/sigs[i])^2d/2d
        mask = (abs(exparg) lt 80)
        yvals[indsubwaves[i,gind]] += transpose(fluxes[i]*mask*EXP(exparg*mask))      
     endif
  endfor

  RETURN, yvals

end
