;
; History
;   08xxx??  HJZ   created
;   09yyy??  DSNR  re-written
;   10nov05  DSNR  added sigma in velocity space
;
function manygauss_slow, wave, param, velsig=velsig

  ppoff = param[0]
  nwave = n_elements(wave)
  nline = (N_ELEMENTS(param)-ppoff)/3

  find = ppoff+indgen(nline)*3
  wind = find + 1
  sind = find + 2

  rrwave = rebin(wave,nwave,nline)
  fluxes = rebin(reform(param[find],1,nline),nwave,nline)
  refwaves = rebin(reform(param[wind],1,nline),nwave,nline)
  sigs = rebin(reform(param[sind],1,nline),nwave,nline)

; Sigma in velocity space?  Then convert to A.
  if keyword_set(velsig) then sigs *=  refwaves/299792d

  yvals = fluxes*EXP(-(rrwave-refwaves)*(rrwave-refwaves)/sigs/sigs/2d)
  ysum = total(yvals,2)

  RETURN, ysum

end
