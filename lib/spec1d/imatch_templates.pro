;+
; NAME:
;       IMATCH_TEMPLATES()
;
; PURPOSE:
;       Match a set of templates to a data spectrum in terms of
;       wavelength spacing and spectral resolution.
;
; CALLING SEQUENCE:
;          newflux = imatch_templates(starflux,starwave,$
;             newwave, [starres=,newres=])
;
; INPUTS:
;       starflux - stellar template flux array [NSTARFLUX,NSTAR] 
;       starwave - wavelength vector corresponding to STARFLUX
;                  [NSTARFLUX]
;       newwave  - desired output wavelength spacing (i.e., of the
;                  data spectrum)
;
; OPTIONAL INPUTS:
;       starres  - FWHM spectral resolution of STARFLUX [Angstrom] 
;       newres   - FWHM spectral resolution of NEWFLUX [Angstrom]
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;       newflux  - resampled and broadened STARFLUX
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;       The templates will be padded with zeros where there is no
;       wavelength overlap.  
;
;       This routine currently only supports a single Gaussian sigma
;       as the broadening kernel.  The slow-down is significant for a
;       wavelength-varying kernel. 
;
; PROCEDURES USED:
;       COMBINE1FIBER, GCONVOLVE()
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 May 8, U of A
;       jm03jun6uofa, excised from IFITSPEC()
;-

function imatch_templates, starflux, starwave, newwave, starres=starres, newres=newres

    npix = n_elements(newwave)

    if (n_elements(starflux) eq 0L) or (n_elements(starwave) eq 0L) or (npix eq 0L) then begin
       print, 'Syntax - newflux = imatch_templates(starflux,starwave,$'
       print, '   newwave, [starres=,newres=])'
       return, -1L
    endif
    
    sdim = size(starflux,/n_dimension)
    if sdim ne 2L then begin
       print, 'STARFLUX must have two dimensions.'
       return, -1L
    endif
    
    ssize = size(starflux,/dimension)
    nstar = ssize[1]
    nstarflux = ssize[0]

    if n_elements(starwave) ne nstarflux then begin
       print, 'STARFLUX and STARWAVE must have the same number of pixels.'
       return, -1L
    endif
    
; conversion factor between Gaussian FWHM and one-sigma width 
    
    fwhm2sig = 2.0*sqrt(2.0*alog(2.0)) ; = 2.35
    
; re-sample the stellar templates onto the new wavelength spacing 

    tempflux = fltarr(npix,nstar)

    for istar = 0L, nstar-1L do begin

       combine1fiber, alog10(starwave), starflux[*,istar], starflux[*,istar]*0.0+1.0, $
         newloglam=alog10(newwave), newflux=outflux
       tempflux[*,istar] = outflux

    endfor

; broaden - we should have a more sophisticated broadening option
    
    disp = starwave[1]-starwave[0]                           ; [Angstrom/pixel]
    sigma = sqrt(newres^2.0 - starres^2.0) / fwhm2sig / disp ; [pixels]

    newflux  = tempflux*0.0

    for j = 0L, nstar-1L do newflux[*,j] = gconvolve(reform(tempflux[*,j]),sigma[0],/edge_truncate)
;   for j = 0L, nstar-1L do for k = 0L, npix-1L do flux[k,j] = (gconvolve(reform(tempflux[*,j]),sigma[k]))[k]

return, newflux
end

