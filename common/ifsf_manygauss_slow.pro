; docformat = 'rst'
;
;+
;
; This function is input to MPFIT, which uses it to compute the
; emission line spectrum of multiple Gaussian emission lines
; simultaneously. This routine assumes nothing about dispersion, and
; thus is slower than MANYGAUSS because it uses larger arrays.
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
;      2009, DSNR, copied base code from Harus Jabran Zahid and re-wrote
;      2010nov05, DSNR, added sigma in velocity space
;      2013nov13, DSNR, switched default sigma to velocity space
;      2013nov13, DSNR, documented, renamed, added license and copyright 
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
function ifsf_manygauss_slow, wave, param

  c = 299792.458d
  
  ppoff = param[0]
  nwave = n_elements(wave)
  nline = (N_ELEMENTS(param)-ppoff)/3

  find = ppoff+indgen(nline)*3
  wind = find + 1
  sind = find + 2

  rrwave = rebin(wave,nwave,nline)
  fluxes = rebin(reform(param[find],1,nline),nwave,nline)
  refwaves = rebin(reform(param[wind],1,nline),nwave,nline)
  sigs = rebin(reform(param[sind],1,nline),nwave,nline)

  dwave = rrwave-refwaves
  sigslam = sigs/c*refwaves
  yvals = fluxes*EXP(-dwave*dwave/sigslam/sigslam/2d)

  ysum = total(yvals,2)

  RETURN, ysum

end
