;+
; NAME:
;       GET_EBV()
;
; PURPOSE:
;       Get E(B-V) from the Balmer decrement and an extinction curve.
;
; CALLING SEQUENCE:
;       ebv = get_ebv(decrement,[decrement_err=,ebv_err=],_extra=extra)
;
; INPUTS:
;       decrement - the observed Balmer decrement F(H-alpha)/F(H-beta)
;
; OPTIONAL INPUTS:
;       decrement_err - error in DECREMENT
;
; KEYWORD PARAMETERS:
;       hgamma - de-redden using the H-alpha/H-gamma line ratio
;       extra  - keywords for K_LAMBDA()
;
; OUTPUTS:
;       ebv     - color excess E(B-V) [mag]
;
; OPTIONAL OUTPUTS:
;       ebv_err - error in EBV [mag]
;
; PROCEDURES USED:
;       K_LAMBDA()
;
; COMMENTS:
;       EBV and EBV_ERR are returned as type double. 
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 May 28, U of A
;       jm02nov6uofa - added K_LAMBDA() compatibility
;       jm03jun13uofa - added HGAMMA keyword
;-

function get_ebv, decrement, decrement_err=decrement_err, ebv_err=ebv_err, $
  hgamma=hgamma, _extra=extra

    ndec = n_elements(decrement)
    if ndec eq 0L then begin
       print, 'Syntax - ebv = get_ebv(decrement,[decrement_err=,ebv_err=],_extra=extra)'
       return, -1
    endif

    ndecerr = n_elements(decrement_err)
    if ndecerr eq 0L then decrement_err = decrement*0.0 else begin
       if ndecerr ne ndec then begin
          print, 'DECREMENT and DECREMENT_ERR have incompatible dimensions.'
          return, -1L
       endif
    endelse

    if keyword_set(hgamma) then begin
       
       wave1 = 4340.464
       wave2 = 6562.80
       R_int = double(2.86/0.468) ; intrinsic H-alpha/H-gamma ratio for type B recombination
    
    endif else begin

       wave1 = 4861.33         ; 4861.0
       wave2 = 6562.80         ; 6563.0
       R_int = double(2.86)    ; intrinsic H-alpha/H-beta ratio for type B recombination

    endelse
    
;   R_int_err = 0.05*R_int      ; 5% error (between 5000-20000 K; Caplan & Deharveng 1986)
    R_int_err = 0.0D
    
    rcurve = k_lambda(wave1,_extra=extra)-k_lambda(wave2,_extra=extra)
    ebv = (alog10(decrement)-alog10(R_int))/(0.4D*rcurve)

    ebv_err = 1.0D/(alog(10)*0.4D*rcurve) * $
      sqrt( (decrement_err/decrement)^2.0 + $
      (R_int_err/R_int)^2.0 )
       
return, ebv
end
