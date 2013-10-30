;
; History
;  09nov24  DSNR  added log(lambda) treatment
;

function disperse_template, template, lambda, sigma, loglam=loglam

  npoints = size(template)
  new_temp = dblarr(npoints[1], npoints[2])

  if keyword_set(loglam) then begin

     uselam = alog(lambda)
     disp = uselam[1] - uselam[0]
     sigma_pix = (sigma / 299792d) / disp

  endif else begin

     uselam = lambda
     disp = uselam[1] - uselam[0]
     sigma_pix = sigma / disp

  endelse

  kernelsize = 10d*sigma_pix
  
; Properly center the kernel!  Otherwise one gets a wavelength shift 
; in the resulting convolved array.
;  kernel = dindgen(kernelsize) - kernelsize/2d
  kernel = [dindgen(kernelsize)-fix(kernelsize),$
            dindgen(kernelsize +1)]
  kernel = exp( -kernel^2d/(2d*sigma_pix^2d) )
  kernel /= total(kernel)

  for i = 0, npoints[2] - 1 do $
     new_temp[0,i] = convol(template[*,i], kernel, /edge_wrap)

  return, new_temp

end
