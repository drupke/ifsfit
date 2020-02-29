function ifsf_extcurve_ccm,lambda,rv=rv
;
; History
;  09sep08  DSNR  created
;  20feb28  DSNR  ported to IFSFIT
;
; Compute A(lambda)/A(V) from the Cardelli, Clayton, & Mathis (1989)
; extinction curve.
;
  if ~ keyword_set(rv) then rv=3.1d

  invlam = 1d/(lambda/10000d)

  lran_ir     = [0.3d,1.1d]
  lran_niropt = [1.1d,3.3d]
  lran_uv     = [3.3d,8.0d]

  acoeff_niropt = double([1,0.17699,-0.50447,-0.02427,0.72085,0.01979,-0.77530,0.32999])
  bcoeff_niropt = double([0,1.41338,2.28305,1.07233,-5.38434,-0.62251,5.30260,-2.09002])

  a = dblarr(n_elements(lambda))
  b = dblarr(n_elements(lambda))
  for i=0,n_elements(lambda)-1 do begin

     if (invlam[i] ge lran_ir[0] AND invlam[i] le lran_ir[1]) then begin
        x = invlam[i]
        a[i] = 0.574d * x^1.61d
        b[i] =-0.527d * x^1.61d
     endif else if (invlam[i] ge lran_niropt[0] AND invlam[i] le lran_niropt[1]) then begin
        y = invlam[i] - 1.82d
        a[i] = poly(y,acoeff_niropt)
        b[i] = poly(y,bcoeff_niropt)
     endif else if (invlam[i] ge lran_uv[0] AND invlam[i] le lran_uv[1]) then begin
        x = invlam[i]
        if x lt 5.9 then begin
           fa = 0d
           fb = 0d
        endif else begin
           y = x - 5.9d
           fa =-0.04473d*y^2d - 0.009779d*y^3d
           fb = 0.2130d *y^2d + 0.1207d  *y^3d
        endelse
        a[i] = 1.752d - 0.316d*x - 0.104d/((x-4.67d)^2d + 0.341d) + fa
        b[i] =-3.090d + 1.825d*x + 1.206d/((x-4.62d)^2d + 0.263d) + fb
     endif else begin
        print,'IFSF_EXTCURVE_CCM: Out of wavelength range of CCM89.  Returning 0.'
        return,0d
     endelse

  endfor

  alamav = a + b/rv

  return,alamav

end
