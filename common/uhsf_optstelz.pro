;
; History
;  13sep25  DSNR  created
;
function uhsf_optstelz,wave,data,model,lrange=lrange,$
                       dzrange=dzrange,dzstep=dzstep,$
                       quiet=quiet

  if ~ keyword_set(dzrange) then dzrange = [-0.005,0.005]
  if ~ keyword_set(dzstep) then dzstep = 0.00002

  if keyword_set(lrange) then begin
     ilrange = where(wave gt lrange[0] AND wave lt lrange[1])
     data_use = data[ilrange]
     model_use = model[ilrange]
     wave_use = wave[ilrange]
  endif else begin
     data_use = data
     model_use = model
     wave_use = wave
  endelse

; grab the dispersion from the middle of the wavelength array
  disp = wave_use[fix(n_elements(wave_use)/2) + 1] - $
         wave_use[fix(n_elements(wave_use)/2)]
  meanwave = mean(wave_use)

; Make lag array (in units of dispersion elements)
;
; # of lags
  nlags = fix((dzrange[1] - dzrange[0])/dzstep + 1)
; Make # of lags odd
  if not nlags then nlags++
; lag array in z space
  lags_z = (dindgen(nlags) - floor(nlags/2)) * dzstep
; lag array in space where 1 unit = 1 dispersion element
  lags = lags_z * meanwave / disp

; Run cross-correlation
  rslt = xcor(data_use,model_use,lags)
  
; Find lag at max. CC in z space
  dzstar = lags_z[rslt[1]]
  
  if not quiet then print,'UHSF_OPTSTELZ: Adding ',$
                          string(dzstar,'(D0.5)'),' to z_star.'

  return,dzstar

end
