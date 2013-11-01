;+
; NAME:
;       DUST_CORRECT()
;
; PURPOSE:
;       Correct a line flux measurement for dust extinction.  
;
; CALLING SEQUENCE:
;       Syntax - newflux = dust_correct(lineflux,linewave,$
;          [lineflux_err=,decrement=,decrement_err=,ebv=],$
;          [ebv_err=,newflux_err=],_extra=extra
;
; INPUTS:
;       lineflux  - emission line flux [NFLUX]
;       linewave  - emission line wavelength (scalar or [NFLUX] array)
;
; OPTIONAL INPUTS:
;       decrement     - Balmer decrement (H-alpha/H-beta or
;                       H-alpha/H-gamma) [NFLUX]
;       errdecrement  - corresponding error in DECREMENT [NFLUX]
;       ebv           - color excess E(B-V) [NFLUX]
;       errebv        - corresponding error in EBV [NFLUX]
;       
;
; KEYWORD PARAMETERS:
;       extra         - keywords for GET_EBV() and K_LAMBDA() 
;
; OUTPUTS:
;       newflux       - extinction-corrected LINEFLUX [NFLUX]
;
; OPTIONAL OUTPUTS:
;       newflux_err   - if either ERRDECREMENT or EBV_ERR are passed
;                       then compute the error in LINEFLUX [NFLUX]
;
; COMMENTS:
;       If both DECREMENT and EBV are passed then DECREMENT is used to
;       calculate and overwrite EBV.
;
; PROCEDURES USED:
;       GET_EBV(), K_LAMBDA()
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 May 28, U of A
;-

function dust_correct, lineflux, linewave, lineflux_err=lineflux_err, $
  decrement=decrement, errdecrement=errdecrement, ebv=ebv, errebv=errebv, $
  newflux_err=newflux_err, _extra=extra

    nflux = n_elements(lineflux)
    nlines = n_elements(linewave)

    if (nflux eq 0L) or (nlines eq 0L) then begin
       print, 'Syntax - newflux = dust_correct(lineflux,linewave,$'
       print, '   [lineflux_err=,decrement=,errdecrement=,ebv=],$'
       print, '   [errebv=,newflux_err=],_extra=extra'
       return, -1L
    endif
    
    if (nlines ne 1L) and (nlines ne nflux) then begin
       print, 'LINEWAVE must be a scalar quantity or an [NFLUX] array.'
       return, -1L
    endif

    nfluxerr = n_elements(lineflux_err)
    if nfluxerr eq 0L then lineflux_err = lineflux*0.0 else begin
       if nfluxerr ne nflux then begin
          print, 'LINEFLUX and LINEFLUX_ERR have incompatible dimensions.'
          return, -1L
       endif 
    endelse
       
    ndec = n_elements(decrement)
    nebv = n_elements(ebv)

    if (ndec eq 0L) and (nebv eq 0L) then begin
       print, 'Either DECREMENT or EBV must be specified.'
       return, -1L
    endif

    if (nebv ne 0L) then begin
       if nebv ne nflux then begin
          print, 'EBV and LINEFLUX have incompatible dimensions.'
          return, -1L
       endif
    endif

    if (ndec ne 0L) then begin
       if ndec ne nflux then begin
          print, 'DECREMENT and LINEFLUX have incompatible dimensions.'
          return, -1L
       endif
; calculate the color excess and error
       ebv = get_ebv(decrement,decrement_err=errdecrement,errebv=errebv,_extra=extra)
    endif

; boost the flux appropriately and compute the error

    kl = k_lambda(linewave,_extra=extra)
    newflux = lineflux*10D0^(0.4*ebv*kl)

; include the error in E(B-V)

;   newflux_err = sqrt( (lineflux_err*10.0^(0.4*ebv*kl))^2.0 + $
;     (lineflux*0.4*kl*alog(10)*10.0^(0.4*ebv*kl)*errebv)^2.0 )

; do not include the error in E(B-V)

    newflux_err = lineflux_err*10.0^(0.4*ebv*kl)

return, newflux
end


