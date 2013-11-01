;+
; NAME:
;       EMISSION_MASK()
;
; PURPOSE:
;       Mask out optical nebular emission lines and telluric
;       absorption features from a spectrum. 
;
; CALLING SEQUENCE:
;       mask = emission_mask(wave,[z=,width=,spectrum=],$
;          [good=,bad=],/telluric,/cosmic)
;
; INPUTS:
;       wave  - wavelength vector [Angstrom]
;
; OPTIONAL INPUTS:
;       z        - dimensionless redshift
;       width    - width of the mask around each emission line
;                  [Angstrom] (default 15.0)
;       spectrum - spectrum corresponding to WAVE.  see COSMIC
;       extra    - keywords for DJS_ITERSTAT
;
; KEYWORD PARAMETERS:
;       telluric - flag wavelengths corresponding to telluric
;                  absorption lines in MASK
;       cosmic   - attempt to identify and flag cosmic rays (deviant
;                  pixels) in SPECTRUM and add them to MASK 
;
; OUTPUTS:
;       mask  - pixel mask for WAVE (1B is good and 0B is bad)
;
; OPTIONAL OUTPUTS:
;       good  - indices in WAVE that have not been masked
;       bad   - indices in WAVE that have been masked
;
; COMMENTS:
;       See Figure 1 in Bessell (1999) for the strongest telluric
;       feature wavelengths.
;
; PROCEDURES USED:
;       DJS_ITERSTAT
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 April 8, U of A
;-

function emission_mask, wave, z=z, width=width, spectrum=spectrum, $
  good=good, bad=bad, telluric=telluric, cosmic=cosmic

    nwave = n_elements(wave)
    if nwave eq 0L then begin
       print, 'Syntax - mask = emission_mask(wave,[z=,width=,spectrum=],$'
       print, '   [good=,bad=],/telluric,/cosmic)'
       return, -1
    endif

    if n_elements(z) eq 0L then z = 0.0
    if n_elements(width) eq 0L then width = 15.0 ; mask width [Angstrom]
    
    mask = make_array(nwave,/byte,value=1)

    linelist = (1+z)*[$
                       3726.032,$ ; OII
                       3728.815,$ ; OII
                       3869.06, $ ; NeIII
;                      3797.898,$ ; H8
;                      3835.384,$ ; H7
;                      3889.049,$ ; H6
;                      3970.072,$ ; H-epsilon
                       4101.734,$ ; H-delta
                       4340.46 ,$ ; H-gamma
                       4363.21 ,$ ; OIII
                       4861.33 ,$ ; H-beta
                       4958.91 ,$ ; OIII
                       5006.84 ,$ ; OIII
                       5875.96 ,$ ; HeI
                       5890.0  ,$ ; Na D doublet
                       5896.0  ,$ ; Na D doublet
                       6300.30 ,$ ; OI
                       6548.04 ,$ ; NII
                       6562.80 ,$ ; H-alpha
                       6583.46 ,$ ; NII
                       6716.14,$  ; SII
                       6730.81$   ; SII
    ]

    for i = 0L, n_elements(linelist)-1L do $
      mask = mask and ((wave lt linelist[i]-width) or (wave gt linelist[i]+width))

; flag telluric features
    
    if keyword_set(telluric) then begin

       bband = where((wave gt 6860.0) and (wave lt 6890.0),nbband)
       if nbband ne 0L then mask[bband] = 0B

       aband = where((wave gt 7600.0) and (wave lt 7630.0),naband)
       if naband ne 0L then mask[aband] = 0B

       tell = where((wave gt 7170.0) and (wave lt 7350.0),ntell)
       if (ntell ne 0L) then mask[tell] = 0B

    endif

; attempt to remove cosmic-rays if requested    

    nspec = n_elements(spectrum)
    if keyword_set(cosmic) and (nspec ne 0L) then begin

       if nspec ne nwave then begin
          print, 'SPECTRUM and WAVE must have the same number of elements.'
       endif else begin

;         res = djs_avsigclip(reform(spectrum,nspec,1),0,inmask=(mask eq 0B),$
;           outmask=crmask,sigrej=3.0)
;         crmask = crmask eq 0B

          noelines = where(mask eq 1B)
          djs_iterstat, spectrum[noelines], sigrej=3.0, mask=crmask, _extra=extra

          mask[noelines] = (mask[noelines] + crmask) gt 1B

       endelse
    endif

    good = where(mask eq 1B,ngood,comp=bad,ncomp=nbad)

return, mask
end
