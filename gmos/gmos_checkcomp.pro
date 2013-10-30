function gmos_checkcomp,struct,sigcut=sigcut

;  fluxthresh = 1d-8 ; GMOS
  fluxthresh = 1d-18 ; KPNO
; 0.8d is approximately resolution for GMOS R600 grating.  10 A in
; sigma is ~1000 km/s FWHM.
  sigthreshlo = 0.8d
  sigthreshhi = 15d
  comppctthresh = 0.05d
  if ~ keyword_set(sigcut) then sigcut = 3d

  ppoff = struct.param[0]
  ncomp = struct.param[1]

  linepars = sepfitpars(struct.param,struct.perror)
  nlines = n_elements(linepars.flux[*,0])

  haind = where(struct.linelabel eq 'Halpha')
  n2ind = where(struct.linelabel eq '[NII]6548')
  goodcomp = -1

check:

  for i=0,ncomp-1 do begin
     haflux = linepars.fluxpk[haind,i]
     hasig = linepars.sigma[haind,i]
     hafluxerr = linepars.fluxpkerr[haind,i]
     n2flux = linepars.fluxpk[n2ind,i]
;;      if (haflux gt 0 AND n2flux gt 0 AND haflux gt fluxthresh $
;;          AND haflux gt sigcut*hafluxerr) then begin
     if (haflux gt 0 AND haflux gt fluxthresh AND $
         haflux gt sigcut*hafluxerr AND $
         haflux gt comppctthresh*linepars.fluxpk[haind,0] AND $
         hasig lt sigthreshhi AND $
         hasig gt sigthreshlo) then begin        
        if goodcomp[0] eq -1 then goodcomp = i $
        else goodcomp = [goodcomp,i]
     endif
  endfor

;;   if goodcomp[0] eq -1 then begin
;;      sigcut -= 0.5d
;;      print,'GMOS_CHECKCOMP: Lowering sigma threshold to ',sigcut,$
;;            format='(A,D0.1)'
;;      goto,check
;;   endif

  return,goodcomp

end
