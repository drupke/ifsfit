;
; History
;   09may26  DSNR  created
;   09jun07  DSNR  added error propagation and rewrote
;
function uhsf_sepfitpars,param,perror,waveran=waveran,fluxpkerr=fluxpkerr

  c = 299792.458d

  ppoff = param[0]
  ncomp = param[1]

  nlines = (n_elements(param)-ppoff)/3/ncomp
  fluxind = (indgen(nlines*ncomp)*3 + ppoff)
  waveind = fluxind+1
  sigind = fluxind+2

; Errors
  if ~ keyword_set(fluxpkerr) then fluxpkerr = perror[fluxind]
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
; Set fluxes to 0 in NaN or infinite errors
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
