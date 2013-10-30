function lrisbr_manygauss, wave, param, floatcomp=floatcomp

  ppoff = param[0]
  nwave = n_elements(wave)
  nline = (N_ELEMENTS(param)-ppoff)/3

  find = ppoff+indgen(nline)*3
  wind = find + 1
  sind = find + 2

  fluxes = param[find]
  refwaves = param[wind]
  sigs = param[sind]
  if keyword_set(floatcomp) then begin
     refwaves = [refwaves,param[7]]
  endif

  dispersion = wave[1] - wave[0]
  maxsig = max([param[10],param[sind]])
  if keyword_set(floatcomp) then begin
     nline+=1
     maxsig = max([maxsig,param[8]])
  endif

  nsubwave = round(10d * maxsig / dispersion)
  halfnsubwave = round(nsubwave / 2)
  nsubwave = halfnsubwave*2+1
  indsubwaves = rebin(transpose(fix(indgen(nsubwave)-halfnsubwave)),$
                      nline,nsubwave)
  indrefwaves_real = (refwaves - wave[0]) / dispersion
  indrefwaves = fix(indrefwaves_real)
  indrefwaves_frac = indrefwaves_real - double(indrefwaves)
  dwaves = (indsubwaves - rebin(indrefwaves_frac,nline,nsubwave))*$
           dispersion
  indsubwaves += rebin(indrefwaves,nline,nsubwave)

  yvals = dblarr(nwave)
  for i=0,nline-1 do begin
     gind = where(indsubwaves[i,*] ge 0 AND indsubwaves[i,*] le nwave-1,$
                  count)
     if count gt 0 then begin
        if ~ keyword_set(floatcomp) OR i ne nline-1 then $
           yvals[indsubwaves[i,gind]] += transpose(fluxes[i]*$
                                                   EXP(-(dwaves[i,gind]/sigs[i])^2d/2d)) $
        else $
           yvals[indsubwaves[i,gind]] += transpose(param[6]*$
                                                   EXP(-(dwaves[i,gind]/param[8])^2d/2d))
     endif
  endfor

  RETURN, yvals

end
