; docformat = 'rst'
;
;+
;
; Convert output of MPFIT, with best-fit line parameters in a single
; array, into a structure with separate arrays for different line
; parameters. Compute total line fluxes from the best-fit line
; parameters.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    A structure with separate arrays for different line
;    parameters. Each array has NxM elements, where N is the number of
;    emission lines fit and M is the number of velocity
;    components. Tags: flux, fluxerr, fluxpk, fluxpkerr, wave, and
;    sigma.
;
; :Params:
;    param: in, required, type=dblarr
;      Best-fit parameter array output by MPFIT.
;    perror: in, optional, type=dblarr
;      Errors in best fit parameters, output by MPFIT.
;
; :Keywords:
;    waveran: in, optional, type=dblarr(2)
;      Set to upper and lower limits to return line parameters only
;      for lines within the given wavelength range. Lines outside this
;      range have fluxes set to 0.
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
;      2009may26, DSNR, created
;      2009jun07, DSNR, added error propagation and rewrote
;      2013nov01, DSNR, added documentation
;      2013nov25, DSNR, renamed, added copyright and license
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
function ifsf_sepfitpars,param,perror,waveran=waveran
  
  c = 299792.458d

  ppoff = param[0]
  ncomp = param[1]

  nlines = (n_elements(param)-ppoff)/3/ncomp
  fluxind = (indgen(nlines*ncomp)*3 + ppoff)
  waveind = fluxind+1
  sigind = fluxind+2

; Errors
  fluxpkerr = perror[fluxind]
  sigerr = dblarr(nlines*ncomp)

; Gaussian fluxes
  fluxpk = param[fluxind]
  gflux = GAUSSFLUX(fluxpk,param[sigind]/c * param[waveind],$
                    normerr=fluxpkerr,sigerr=sigerr)

; Set fluxes to 0 outside of wavelength range
  if keyword_set(waveran) then begin
     inoflux = where(waveran[0] gt param[waveind]*(1 - 3d*param[sigind]/c) OR $
                     waveran[1] lt param[waveind]*(1 + 3d*param[sigind]/c),ct)
     if ct gt 0 then begin
        gflux.flux[inoflux] = 0d
        gflux.flux_err[inoflux] = 0d
        fluxpk[inoflux] = 0d
        fluxpkerr[inoflux] = 0d
     endif
  endif
; Set fluxes to 0 if NaN or infinite errors
  inans = where(finite(gflux.flux_err) eq 0,ctnan)
  if ctnan gt 0 then begin
     gflux.flux[inans] = 0d
     gflux.flux_err[inans] = 0d
  endif
  inans = where(finite(fluxpkerr) eq 0,ctnan)
  if ctnan gt 0 then begin
     fluxpk[inans] = 0d
     fluxpkerr[inans] = 0d
  endif

  outstr = { $
           flux:reform(gflux.flux,nlines,ncomp),$
           fluxerr:reform(gflux.flux_err,nlines,ncomp),$
           fluxpk:reform(fluxpk,nlines,ncomp),$
           fluxpkerr:reform(fluxpkerr,nlines,ncomp),$
           wave:reform(param[waveind],nlines,ncomp),$
           sigma:reform(param[sigind],nlines,ncomp) $
           }
  
  return,outstr

end
