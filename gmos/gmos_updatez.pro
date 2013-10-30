pro gmos_updatez,param,wave,zarr

  ppoff = param[0]
  ncomp = param[1]
  ppoff0 = ppoff - 2*(ncomp-1)

; Fitted updates to redshifts
  if ncomp gt 1 then $
     for i=0,ncomp-2 do zarr[i+1]=zarr[0]+param[ppoff0+2*i+1]
  zarr[0] += param[2]

  wave /= 1d + param[2]

  nlines = (n_elements(param)-ppoff)/3
  waveind = ppoff + indgen(nlines)*3 + 1
  param[waveind] /= 1d + param[2]
     
end
