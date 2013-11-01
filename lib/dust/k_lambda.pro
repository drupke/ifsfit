;+
; NAME:
;       K_LAMBDA()
;
; PURPOSE:
;       Return a variety of reddening curves.
;
; CALLING SEQUENCE:
;       k_lambda = k_lambda(wave,[r_v=],calzetti=calzetti,$
;                     ccm=ccm,lmc=lmc,smc=smc)
;
; INPUTS:
;       wave  - wavelength vector [Angstroms]
;
; OPTIONAL INPUTS:
;       r_v   - total to selective extinction ratio (default 3.1)
;
;
; KEYWORD PARAMETERS:
;       calzetti - return the Calzetti (1994) starburst extinction
;                  curve [NB: continuum extinction curve]
;       charlot  - return the Charlot & Fall (2000) extinction curve
;       ccm      - return the Cardelli, Clayton, & Mathis 1989
;                  Galactic extinction curve
;       lmc      - return the Large Magellanic Cloud extinction curve
;       smc      - return the Small Magellanic Cloud extinction curve
;
; OUTPUTS:
;       k_lambda - defined as A(lambda)/E(B-V)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       Much of the code for this routine has been excised from the
;       Goddard routines CCM_UNRED(), FM_UNRED(), and CALZ_UNRED() for
;       my own nefarious purposes.
;
; PROCEDURES USED:
;       READ_SMC()
;
; EXAMPLE:
;       Plot the calzetti and the CCM extinction curves in the
;       optical and UV.
;
;       IDL> wave = findgen(6000)+1000.0 ; Angstrom
;       IDL> plot, wave, k_lambda(wave,/calzetti), xsty=3, ysty=3
;       IDL> oplot, wave, k_lambda(wave,/ccm)
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 September 19, U of A
;       jm03mar2uofa - added the Charlot & Fall (2000) curve
;-

function k_lambda, wave, r_v=r_v, calzetti=calzetti, charlot=charlot, $
  ccm=ccm, lmc=lmc, smc=smc

    if n_elements(wave) eq 0L then begin
       print, 'Syntax - k_lambda = k_lambda(wave,[r_v=],calzetti=calzetti,$'
       print, '   charlot=charlot,ccm=ccm,lmc=lmc,smc=smc)'
       return, -1
    endif

    k_lambda = -1.0

; select the default extinction curve (CCM)
    
    if (n_elements(calzetti) eq 0L) and (n_elements(charlot) eq 0L) and $
      (n_elements(ccm) eq 0L) and (n_elements(lmc) eq 0L) and $
      (n_elements(smc) eq 0L) then ccm = 1L
    
    if keyword_set(calzetti) then begin ; Calzetti

       if n_elements(r_v) eq 0 then r_v = 4.05
       w1 = where((wave ge 6300) and (wave le 22000),c1)
       w2 = where((wave ge  912) and (wave lt  6300),c2)
       x  = 10000.0/wave        ; wavelength in inverse microns

       if (c1 + c2) ne n_elements(wave) then message, /info, $
         'Warning - some elements of wavelength vector outside valid domain.'

       k_lambda = 0.0*wave

       if c1 gt 0 then $
         k_lambda[w1] = 2.659*(-1.857 + 1.040*x[w1]) + R_V
       
       if c2 gt 0 then $
         k_lambda[w2] = 2.659*(poly(x[w2], [-2.156, 1.509d0, -0.198d0, 0.011d0])) + R_V
       
    endif

    if keyword_set(charlot) then begin ; Charlot & Fall (2000)

       if n_elements(r_v) eq 0L then r_v = 3.1
       k_lambda = r_v*(wave/5500.0)^(-0.7)
       
    endif
    
    if keyword_set(ccm) then begin ; Cardelli, Clayton, & Mathis 1989

       if n_elements(r_v) eq 0L then r_v = 3.1

       x = 10000./ wave         ; convert to inverse microns 
       npts = n_elements(x)
       a = fltarr(npts)  
       b = fltarr(npts)

       good = where( (x gt 0.3) and (x lt 1.1),ngood ) ; infrared
       if ngood gt 0 then begin
          a[good] =  0.574 * x[good]^(1.61)
          b[good] = -0.527 * x[good]^(1.61)
       endif

       good = where( (x ge 1.1) and (x lt 3.3),ngood) ; optical/NIR
       if ngood GT 0 then begin

          y = x[good] - 1.82
;         c1 = [ 1. , 0.17699, -0.50447, -0.02427,  0.72085,    $ ;Original
;           0.01979, -0.77530,  0.32999 ] ;coefficients
;         c2 = [ 0.,  1.41338,  2.28305,  1.07233, -5.38434,    $ ;from CCM89
;           -0.62251,  5.30260, -2.09002 ]
          c1 = [ 1. , 0.104,   -0.609,    0.701,  1.137,    $ ;New coefficients
            -1.718,   -0.827,    1.647, -0.505 ] ;from O'Donnell (1994)
          c2 = [ 0.,  1.952,    2.908,   -3.989, -7.985,    $
            11.102,    5.491,  -10.805,  3.347 ]

          a[good] = poly( y, c1)
          b[good] = poly( y, c2)

       endif

       good = where( (x ge 3.3) and (x lt 8),ngood) ; mid-UV
       if ngood gt 0 then begin

          y = x[good]
          F_a = fltarr(Ngood)    & F_b = fltarr(Ngood)
          good1 = where( (y GT 5.9), Ngood1 )
          if Ngood1 GT 0 then begin
             y1 = y[good1] - 5.9
             F_a[ good1] = -0.04473 * y1^2 - 0.009779 * y1^3
             F_b[ good1] =   0.2130 * y1^2  +  0.1207 * y1^3
          endif
          
          a[good] =  1.752 - 0.316*y - (0.104 / ( (y-4.67)^2 + 0.341 )) + F_a
          b[good] = -3.090 + 1.825*y + (1.206 / ( (y-4.62)^2 + 0.263 )) + F_b

       endif

       good = where( (x GE 8) and (x LE 11),ngood ) ; far-UV
       if ngood gt 0 then begin

          y = x[good] - 8.
          c1 = [ -1.073, -0.628,  0.137, -0.070 ]
          c2 = [ 13.670,  4.257, -0.420,  0.374 ]
          a[good] = poly(y, c1)
          b[good] = poly(y, c2)

       endif

       k_lambda = r_v * (a + b/r_v)

    endif

    if keyword_set(lmc) then begin

       if n_elements(r_v) eq 0L then r_v = 3.1

       x = 10000./ wave         ; convert to inverse microns 
       k_lambda = x*0.

       if n_elements(x0) EQ 0 then x0 = 4.596  
       if n_elements(gamma) EQ 0 then gamma = 0.91
       if n_elements(c4) EQ 0 then c4   =  0.64  
       if n_elements(c3) EQ 0 then c3    =  2.73	
       if n_elements(c2) EQ 0 then c2    = 1.11
       if n_elements(c1) EQ 0 then c1    =  -1.28

; compute UV portion of A(lambda)/E(B-V) curve using FM fitting
; function and R-dependent coefficients 
       
       xcutuv = 10000.0/2700.0
       xspluv = 10000.0/[2700.0,2600.0]
       iuv = where(x ge xcutuv, N_UV)
       IF (N_UV GT 0) THEN xuv = [xspluv,x[iuv]] ELSE  xuv = xspluv

       yuv = c1  + c2*xuv
       yuv = yuv + c3*xuv^2/((xuv^2-x0^2)^2 +(xuv*gamma)^2)
       yuv = yuv + c4*(0.5392*((xuv>5.9)-5.9)^2+0.05644*((xuv>5.9)-5.9)^3)
       yuv = yuv + R_V
       yspluv  = yuv[0:1]       ; save spline points

       IF (N_UV GT 0) THEN k_lambda[iuv] = yuv[2:*] ; remove spline points
       
; compute optical portion of A(lambda)/E(B-V) curve using cubic spline
; anchored in UV, optical, and IR 

       xsplopir = [0,10000.0/[26500.0,12200.0,6000.0,5470.0,4670.0,4110.0]]
       ysplir   = [0.0,0.26469,0.82925]*R_V/3.1 
       ysplop   = [poly(R_V, [-4.22809e-01, 1.00270, 2.13572e-04] ), $
         poly(R_V, [-5.13540e-02, 1.00216, -7.35778e-05] ), $
         poly(R_V, [ 7.00127e-01, 1.00184, -3.32598e-05] ), $
         poly(R_V, [ 1.19456, 1.01707, -5.46959e-03, 7.97809e-04, $ 
         -4.45636e-05] ) ]
       
       ysplopir = [ysplir,ysplop]

       iopir = where(x lt xcutuv, Nopir)
       if (Nopir GT 0) then $
         k_lambda[iopir] = CSPLINE([xsplopir,xspluv],[ysplopir,yspluv],x[iopir])

    endif

    if keyword_set(smc) then begin

       if n_elements(r_v) eq 0L then r_v = 3.1
       k_lambda = read_smc(wave,r_v=r_v)
    
    endif

    if n_elements(k_lambda) eq 1L then k_lambda = k_lambda[0]
    
return, k_lambda
end
