;
; History
;   13mar05  DSNR  created
;
; Re-order emission lines so that first component is one with largest
; peak flux.
;
pro nifs_orderlines,struct,doubleline=doubleline

  ppoff = struct.param[0]
  ncomp = struct.param[1]
  ppoff0 = ppoff - (ncomp-1)

; Extract line fit parameters
  linepars = sepfitpars(struct.param,struct.perror)
; Number of fitted lines
  nlines = n_elements(linepars.flux[*,0])

  if doubleline eq 'paa' then begin

; Get Paa fluxes for each component, and sort in descending order
     paaind = where(struct.linelabel eq 'paa')
     paafluxes = dblarr(ncomp-1)
     paafluxes[0] = linepars.fluxpk[paaind,0]
     paafluxes[1] = linepars.fluxpk[paaind,2]
     sortpaa = reverse(sort(paafluxes))

; Reorder emission-line redshifts in fit structure
     ptmparr = [struct.param[ppoff0-1],$
                struct.param[ppoff-1]]
     ptmparr = ptmparr[sortpaa]
     struct.param[ppoff0-1] = ptmparr[0]
     struct.param[ppoff-1] = ptmparr[1]

; Reorder fitted line parameters.  Note that wavelength and sigma
; errors are not reordered, but we don't use the wavelength error and
; the sigma error is fixed in sepfitpars.
     for i=0,1 do begin
        if i eq 0 then foff = ppoff + indgen(nlines)*3 $
        else foff = ppoff + 2*nlines*3 + indgen(nlines)*3
        woff = foff + 1
        soff = foff + 2
        if sortpaa[i] eq 1 then sortpaa[i] = 2
        struct.param[foff]  = linepars.fluxpk[*,sortpaa[i]]
        struct.perror[foff] = linepars.fluxpkerr[*,sortpaa[i]]
        struct.param[woff]  = linepars.wave[*,sortpaa[i]]
        struct.param[soff]  = linepars.sigma[*,sortpaa[i]]
     endfor

  endif else begin

; Get H2 S1 fluxes for each component, and sort in descending order
     h2s1ind = where(struct.linelabel eq 'H2_10_S1')
     h2s1fluxes = dblarr(ncomp-1)
     for i=1,ncomp-1 do $
        h2s1fluxes[i-1] = linepars.fluxpk[h2s1ind,i]
     sorth2s1 = reverse(sort(h2s1fluxes))

; Reorder emission-line redshifts in fit structure
     ptmparr = [struct.param[ppoff0:ppoff-1]]
     ptmparr = ptmparr[sorth2s1]
     struct.param[ppoff0:ppoff-1] = ptmparr[0:ncomp-2]
  
; Reorder fitted line parameters.  Note that wavelength and sigma
; errors are not reordered, but we don't use the wavelength error and
; the sigma error is fixed in sepfitpars.
     for i=1,ncomp-1 do begin
        foff = ppoff + i*nlines*3 + indgen(nlines)*3
        woff = foff + 1
        soff = foff + 2
        struct.param[foff]  = linepars.fluxpk[*,sorth2s1[i-1]+1]
        struct.perror[foff] = linepars.fluxpkerr[*,sorth2s1[i-1]+1]
        struct.param[woff]  = linepars.wave[*,sorth2s1[i-1]+1]
        struct.param[soff]  = linepars.sigma[*,sorth2s1[i-1]+1]
     endfor

  endelse

end
