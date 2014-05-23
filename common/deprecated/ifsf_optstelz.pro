; docformat = 'rst'
;
;+
;
; Find best-fit stellar redshift via cross-correlation with a stellar
; model. Data and model are assumed to be either both in the rest
; frame or both in the observed frame.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Best fit redshift, as a redshift difference between the data and
;    model.
;
; :Params:
;    wave: in, required, type=dblarr(N)
;      Wavelengths.
;    data: in, required, type=dblarr(N)
;      Flux of spectrum to fit.
;    model: in, required, type=dblarr(N)
;      Stellar model.
;
; :Keywords:
;    lrange: in, optional, type=dblarr(2), default=[0\,Inf]
;      Wavelength range over which to do x-correlation.
;    dzrange: in, optional, type=dblarr(2), default=[-0.005\,0.005]
;      Range of lags, in redshift space.
;    dzstep: in, optional, type=double, default=0.00002
;      Lag step, in redshift space.
;    quiet: in, optional, type=byte, default=0
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
;      2013sep25, DSNR, created
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
function ifsf_optstelz,wave,data,model,lrange=lrange,$
                       dzrange=dzrange,dzstep=dzstep,$
                       quiet=quiet

  if ~ keyword_set(dzrange) then dzrange = [-0.005,0.005]
  if ~ keyword_set(dzstep) then dzstep = 0.00002

  if keyword_set(lrange) then begin
     ilrange = where(wave gt lrange[0] AND wave lt lrange[1])
     data_use = data[ilrange]
     model_use = model[ilrange]
     wave_use = wave[ilrange]
  endif else begin
     data_use = data
     model_use = model
     wave_use = wave
  endelse

; grab the dispersion from the middle of the wavelength array
  disp = wave_use[fix(n_elements(wave_use)/2) + 1] - $
         wave_use[fix(n_elements(wave_use)/2)]
  meanwave = mean(wave_use)

; Make lag array (in units of dispersion elements)
;
; # of lags
  nlags = fix((dzrange[1] - dzrange[0])/dzstep + 1)
; Make # of lags odd
  if not nlags then nlags++
; lag array in z space
  lags_z = (dindgen(nlags) - floor(nlags/2)) * dzstep
; lag array in space where 1 unit = 1 dispersion element
  lags = lags_z * meanwave / disp

; Run cross-correlation
  rslt = xcor(data_use,model_use,lags)
  
; Find lag at max. CC in z space
  dzstar = lags_z[rslt[1]]
  
  if not quiet then print,'IFSF_OPTSTELZ: Adding ',$
                          string(dzstar,'(D0.5)'),' to z_star.'

  return,dzstar

end
