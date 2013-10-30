; docformat = 'rst'
;
;+
;
; Remove data near emission lines.
;
; :Returns:
;    Array of lambda-array indices indicating non-masked wavelengths.
;
; :Params:
;    lambda: in, required, type=dblarr
;      Wavelengths of spectrum
;    linelambda: in, required, type=dblarr(nlines)
;      Central wavelengths of lines to mask
;    halfwidth: in, required, type=dblarr(nlines)
;      Half width (in km/s) of masking region around each line
;
; :Keywords:
;    nomaskran: in, optional, type=dblarr(2\,nreg)
;      Set of lower and upper wavelength limits of regions *not* to
;      mask.
;
; :Author:
;    Jabran Zahid and David Rupke
;
; :History:
;    ChangeHistory::
;      2008, HJZ, created
;      2013oct, DSNR, added nomaskran keyword
;-
function uhsf_masklin, lambda, linelambda, halfwidth, nomaskran=nomaskran

  c = 299792.458d

  good = indgen(n_elements(lambda))
  nn = n_elements(linelambda)


  for i = 0, nn - 1 do begin
     gg = where(lambda[good] lt linelambda[i]*(1 - halfwidth[i]/c) or $
                lambda[good] gt linelambda[i]*(1 + halfwidth[i]/c),ct)
     if ct gt 0 then good = good[gg]
  endfor

  if keyword_set(nomaskran) then begin
     for j=0,n_elements(nomaskran)/2 - 1 do begin
        goodadd = where(lambda ge nomaskran[0,j] AND $
                        lambda le nomaskran[1,j],ctadd)
        if ctadd gt 0 then begin
           if ct gt 0 then good = [good,goodadd] else good = goodadd
           ct += ctadd
        endif
     endfor
  endif

  good = good[sort(good)]
  good = good[uniq(good)]

  return, good

end
