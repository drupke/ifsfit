;
; History
;   10may27  DSNR  copied from LRIS version
;
; Re-order emission lines so that first component is one with largest
; Halpha flux.  Presently, stellar emission redshift is determined
; based on component that is first in the redshift and fit arrays.
;
pro gmos_orderlines,struct

  ppoff = struct.param[0]
  ncomp = struct.param[1]
  ppoff0 = ppoff - (ncomp-1)

; Extract line fit parameters
  linepars = sepfitpars(struct.param,struct.perror)
; Number of fitted lines
  nlines = n_elements(linepars.flux[*,0])

; Get Halpha fluxes for each component, and sort in descending order
  haind = where(struct.linelabel eq 'Halpha')
  hafluxes = dblarr(ncomp)
  for i=0,ncomp-1 do $
     hafluxes[i] = linepars.fluxpk[haind,i]
  sortha = reverse(sort(hafluxes))

; Reorder emission-line redshifts in fit structure
  ptmparr = [struct.param[2],struct.param[ppoff0:ppoff-1]]
  ptmparr = ptmparr[sortha]
  struct.param[2] = ptmparr[0]
  struct.param[ppoff0:ppoff-1] = ptmparr[1:ncomp-1]
  
; Reorder fitted line parameters.  Note that wavelength and sigma
; errors are not reordered, but we don't use the wavelength error and
; the sigma error is fixed in sepfitpars.
  for i=0,ncomp-1 do begin
     foff = ppoff + i*nlines*3 + indgen(nlines)*3
     woff = foff + 1
     soff = foff + 2
     struct.param[foff]  = linepars.fluxpk[*,sortha[i]]
     struct.perror[foff] = linepars.fluxpkerr[*,sortha[i]]
     struct.param[woff]  = linepars.wave[*,sortha[i]]
     struct.param[soff]  = linepars.sigma[*,sortha[i]]
  endfor

end
