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
;    A structure with separate hashes for different line parameters. The hashes 
;    are indexed by line, and each value is an array over components. 
;    Tags: flux, fluxerr, fluxpk, fluxpkerr, nolines, wave, and sigma.
;
; :Params:
;    linelist: in, required, type=hash(lines)
;      List of emission line rest-frame wavelengths.
;    param: in, required, type=dblarr(N)
;      Best-fit parameter array output by MPFIT.
;    perror: in, optional, type=dblarr(N)
;      Errors in best fit parameters, output by MPFIT.
;    parinfo: in, required, type=structure
;      Structure input into MPFIT. Each tag has N values, one per parameter. 
;      Used to sort param and perror arrays.
;
; :Keywords:
;    waveran: in, optional, type=dblarr(2)
;      Set to upper and lower limits to return line parameters only
;      for lines within the given wavelength range. Lines outside this
;      range have fluxes set to 0.
;    tflux: out, optional, type=structure
;      A structure with separate hashes for total flux and error. The hashes 
;      are indexed by line. Tags: tflux, tfluxerr
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
;      2014jan13, DSNR, re-written to use hashes rather than arrays
;      2014feb26, DSNR, replaced ordered hashes with hashes
;      2015sep20, DSNR, compute total line flux and error
;    
; :Copyright:
;    Copyright (C) 2013-2015 David S. N. Rupke
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
function ifsf_sepfitpars,linelist,param,perror,parinfo,waveran=waveran,$
                         tflux=tflux

; Return 0 if no lines were fit
  if n_elements(param) eq 1 then begin
    outstr = {nolines:0}
    goto,nolines
  endif

  c = 299792.458d
  maxncomp = param[1]

; If this turns out to be slow, can try creating the hashes by concatenation 
; instead.
  flux = hash(linelist->keys())
  fluxerr = hash(linelist->keys())
  fluxpk = hash(linelist->keys())
  fluxpkerr = hash(linelist->keys())
  wave = hash(linelist->keys())
  sigma = hash(linelist->keys())

  if keyword_set(tflux) then begin
     tf = hash(linelist->keys())
     tfe = hash(linelist->keys())
  endif

  foreach line,linelist->keys() do begin
    fluxpk[line] = param[where(parinfo.line eq line AND $
                               parinfo.parname eq 'flux_peak')]
    fluxpkerr[line] = perror[where(parinfo.line eq line AND $
                                   parinfo.parname eq 'flux_peak')]
    wave[line] = param[where(parinfo.line eq line AND $
                             parinfo.parname eq 'wavelength')]
    sigma[line] = param[where(parinfo.line eq line AND $
                             parinfo.parname eq 'sigma')]
;   Compute total Gaussian flux
    gflux = IFSF_GAUSSFLUX(fluxpk[line],sigma[line]/c * wave[line],$
                           normerr=fluxpkerr[line],sigerr=dblarr(maxncomp))
    flux[line] = gflux.flux
    fluxerr[line] = gflux.flux_err

;   Set fluxes to 0 outside of wavelength range, or if NaNs or infinite errors
    if keyword_set(waveran) then begin
      inoflux = where(waveran[0] gt wave[line]*(1 - 3d*sigma[line]/c) OR $
                      waveran[1] lt wave[line]*(1 + 3d*sigma[line]/c) OR $
                      finite(fluxerr[line]) eq 0 OR $
                      finite(fluxpkerr[line]) eq 0,ct)
      if ct gt 0 then begin
        flux[line,inoflux] = 0d
        fluxerr[line,inoflux] = 0d
        fluxpk[line,inoflux] = 0d
        fluxpkerr[line,inoflux] = 0d
      endif
    endif

    igd = where(flux[line] gt 0d,ctgd)
    if keyword_set(tflux) then begin
       if ctgd gt 0 then begin
         tf[line] = total(flux[line,igd])
         tfe[line] = sqrt(total(fluxerr[line,igd]^2d))
       endif else begin
         tf[line] = 0d
         tfe[line] = 0d
       endelse
    endif
  endforeach


  outstr = {nolines:0,flux:flux,fluxerr:fluxerr,fluxpk:fluxpk,$
            fluxpkerr:fluxpkerr,wave:wave,sigma:sigma}
  if keyword_set(tflux) then tflux = {tflux:tf,tfluxerr:tfe}

nolines:

  return,outstr

end
