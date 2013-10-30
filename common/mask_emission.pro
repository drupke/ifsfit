;+
; NAME:
;     MASK_EMISSION
;
; PURPOSE:
;     Remove data near emission lines.
;
; EXPLANATION:
;
; CALLING SEQUENCE
;     mask_emission,lambda,linelambda,halfwidth
;
; INPUTS:
;     lambda - array of wavelengths
;     linelambda - array of central wavelengths of lines to mask
;     halfwidth - array of half widths of masking region around each line
;
; OUTPUT:
;     Array of lambda-array indices indicating non-masked wavelengths.
;
; METHOD:
;
; REVISION HISTORY:
;     08fall  HJZ  created
;-

function mask_emission, lambda, linelambda, halfwidth

  good = indgen(n_elements(lambda))
  nn = n_elements(linelambda)

  for i = 0, nn - 1 do begin
     gg = where(lambda[good] lt linelambda[i] - halfwidth[i] or $
                lambda[good] gt linelambda[i] + halfwidth[i],ct)
     if ct gt 0 then good = good[gg]
  endfor

  return, good

end
