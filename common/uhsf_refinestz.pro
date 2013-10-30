;
; History
;  11marXY  JTF   created
;  11may09  DSNR  ported to generic function
;

function uhsf_refinestz,zinit,wave,flux,err,sigma=sigma

  carest = 3934.777

; Determine Ca continuum
  w1 = (1d +zinit)*carest -20d
  w2 = (1d +zinit)*carest -12d
  w3 = (1d +zinit)*carest +12d
  w4 = (1d +zinit)*carest +20d
  caindex = where((wave GT w1 AND wave LT w2) OR $
                  (wave GT w3 AND wave LT w4),ct)
  if caindex[0] ne -1 then begin
     p = poly_fit(wave[caindex],flux[caindex],1,$
                  measure_errors=1d/err[caindex]^2d )

; Subtract continuum
     caab = where(wave GT w2 AND wave LT w3)
     cacontflux = poly(wave[caab],p)
     castart = flux[caab] - cacontflux
  
; Fit Gaussian to CaII line
     cafit = mpfitpeak(wave[caab],castart,a,$
                       error=err[caab],/NEGATIVE,nterms=3)

; Calculate redshift based on Ca fitting
     caredshift = (a[1] / carest) -1d
     sigma = a[2]

  endif else begin

     caredshift = zinit

  endelse

  return,caredshift

end
