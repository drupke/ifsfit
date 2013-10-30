function lris_feb09_sepfitpars,param,perror,stitchwave,$
                               waveran=waveran,fluxpkerr=fluxpkerr
;
; History
;   09may26  DSNR  created
;   09jun07  DSNR  added error propagation and rewrote
;

; Estimate for fractional flux contribution from broad line component
  broadfrac = 0.1d

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
  gflux = gaussflux(fluxpk,param[sigind],$
                    normerr=fluxpkerr,sigerr=sigerr)
  ired = where(param[waveind] gt stitchwave,ctred)
  gflux.flux[ired] *= 1d + broadfrac
  gflux.flux_err[ired] *= 1d + broadfrac
  gfluxred = gaussflux(fluxpk[ired]*param[11],replicate(param[12],ctred),$
                       normerr=fluxpkerr[ired]*param[11],$
                       sigerr=sigerr[ired])
  gflux.flux[ired] = gflux.flux[ired] + gfluxred.flux
  gflux.flux_err[ired] = sqrt(gflux.flux_err[ired]^2+gfluxred.flux_err^2)
; Set fluxes to 0 outside of wavelength range
  if keyword_set(waveran) then begin
     inoflux = where(waveran[0] gt param[waveind]-3d*param[sigind] OR $
                    waveran[1] lt param[waveind]+3d*param[sigind],ct)
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
