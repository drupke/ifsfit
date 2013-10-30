;
; History
;   09yyy??  DSNR  copied from manygauss_slow.pro
;   10dec07  DSNR  copied from manygauss.pro
;                  edited to include polynomial baseline
;

; NOTE!  This routine assumes constant dispersion (in A/pix).  If you
; don't have it, either modify this routine to your needs or use
; MANYGAUSS_SLOW.

function gmos_manygauss, wave, param

  ppoff = param[0]
  nwave = n_elements(wave)
  nline = (N_ELEMENTS(param)-ppoff)/3

; output array
  yvals = dblarr(nwave)

; polynomial
  npoly = param[4]
  wavelo = param[5]
  wavehi = param[6]
  n2hawave = where(wave ge wavelo AND wave le wavehi)
  yvals[n2hawave] = poly(wave[n2hawave]-mean(wave[n2hawave]),param[7:7+npoly])

  find = ppoff+indgen(nline)*3
  wind = find + 1
  sind = find + 2

  dispersion = wave[1] - wave[0]
  maxsig = max(param[sind])

  nsubwave = round(10d * maxsig / dispersion)
  halfnsubwave = round(nsubwave / 2)
  nsubwave = halfnsubwave*2+1
  indsubwaves = rebin(transpose(fix(indgen(nsubwave)-halfnsubwave)),$
                      nline,nsubwave)

  fluxes = param[find]
  sigs = param[sind]
  refwaves = param[wind]
  indrefwaves_real = (refwaves - wave[0]) / dispersion
  indrefwaves = fix(indrefwaves_real)
  indrefwaves_frac = indrefwaves_real - double(indrefwaves)
  dwaves = (indsubwaves - rebin(indrefwaves_frac,nline,nsubwave))*$
           dispersion
  indsubwaves += rebin(indrefwaves,nline,nsubwave)
  
  for i=0,nline-1 do begin
     gind = where(indsubwaves[i,*] ge 0 AND indsubwaves[i,*] le nwave-1,$
                  count)
     if count gt 0 then $
        yvals[indsubwaves[i,gind]] += transpose(fluxes[i]*$
                                      EXP(-(dwaves[i,gind]/sigs[i])^2d/2d))
  endfor

  RETURN, yvals

end
