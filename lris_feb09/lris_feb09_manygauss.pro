function lris_feb09_manygauss,wave,param, $
                              floatcomp=floatcomp, $
                              stitchwave=stitchwave

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
;;     ifluxmax = where(fluxes eq max(fluxes))
     refwaves = [refwaves,param[14]]
  endif

  dispersion = wave[1] - wave[0]
  maxsig = max([param[9],param[sind]])
  if keyword_set(floatcomp) then begin
     nline+=1
     maxsig = max([maxsig,param[15]])
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
        if ~ keyword_set(floatcomp) OR i ne nline-1 then begin
           yvals[indsubwaves[i,gind]] += transpose(fluxes[i]*$
                                                   EXP(-(dwaves[i,gind]/sigs[i])^2d/2d))
           if refwaves[i] gt stitchwave then begin
              yvals[indsubwaves[i,gind]] += transpose(fluxes[i]*param[8]*$
                                                      EXP(-(dwaves[i,gind]/param[9])^2d/2d))
              yvals[indsubwaves[i,gind]] += transpose(fluxes[i]*param[11]*$
                                                      EXP(-((dwaves[i,gind]-param[10])/param[12])^2d/2d))
           endif
        endif else begin
           yvals[indsubwaves[i,gind]] += transpose(param[13]*$
                                                   EXP(-(dwaves[i,gind]/param[15])^2d/2d))
        endelse
     endif
  endfor

  RETURN, yvals

end
