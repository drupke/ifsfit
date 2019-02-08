; docformat = 'rst'
;
;+
;
;  Compute E(B-V) and its error from the Cardelli, Clayton, & Mathis
;  (1989) extinction curve using two fluxes.
;
;  Root equation is
;      m_i = m_o - (A/A_V * E(B-V) * R_V)
;  Taking two wavelengths, differencing, and rearranging:
;      E(B-V) = [(m1 - m2)_o - (m1 - m2)_i] / [R_V (A1/AV - A2/AV)]
;
;  Substitute fluxes for mags. to get the equations in the source code.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    N x 1 or N x 2 array, where N is the number of input spectra and the other
;    dimension is the calculated E(B-V) (and its error, if fluxerr is specified).
;
; :Params:
;    lambdas: in, required, type=dblarr(2)
;      Rest wavelengths of lines from which to compute extinction
;    fluxes: in, required, type=dblarr(N,2)
;      Fluxes of two lines; first dimension allows more than one spectrum.
;    predfluxrat: in, required, type=double
;      Predicted flux ratio for E(B-V) = 0
;      
; :Keywords:
;    rv: in, optional, type=double, default=3.1
;      R_V = A_V / E(B-V)
;    fluxerr: in, optional, type=dblarr(N,2)
;      Flux errors of two lines; first dimension allows more than one spectrum.
;      Error in E(B-V) is calculated only if this keyword is present.
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
;      2009sep08, DSNR, created
;      2015apr17, DSNR, copied to IFSFIT package; added documentation, license,
;                       and copyright
;
; :Copyright:
;    Copyright (C) 2014 David S. N. Rupke
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
function ifsf_ebv_ccm,lambdas,fluxes,predfluxrat,rv=rv,fluxerr=fluxerr

  if ~ keyword_set(rv) then rv=3.1d

; Find A_lambda/A_V for each input lambda
  alamav1 = ifsf_extcurve_ccm(lambdas[0],rv=rv)
  alamav2 = ifsf_extcurve_ccm(lambdas[1],rv=rv)
  alamdiff = alamav1 - alamav2

; Arrange arrays 
  sfluxes = size(fluxes)
  if (sfluxes[0] eq 1) then begin
     fluxes = transpose(fluxes)
     if keyword_set(fluxerr) then fluxerr = transpose(fluxerr)
  endif

; Compute E(B-V) and error
  fluxrat = fluxes[*,0] / fluxes[*,1]
  ebv = -2.5d / rv / alamdiff[0] * alog10(fluxrat/predfluxrat)
  if keyword_set(fluxerr) then begin
     fluxraterr = fluxrat*sqrt((fluxerr[*,0]/fluxes[*,0])^2d + $
                               (fluxerr[*,1]/fluxes[*,1])^2d)
     ebverr = abs(-2.5d / rv / alamdiff[0] / alog(10d) * $
                  fluxraterr/fluxrat)
     ebv = [[ebv],[ebverr]]
  endif

; Rearrange arrays if necessary
  if (sfluxes[0] eq 1) then begin
     fluxes = transpose(fluxes)
     ebv = transpose(ebv)
     if keyword_set(fluxerr) then fluxerr = transpose(fluxerr)
  endif

  return,ebv

end
