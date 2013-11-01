;+
; NAME:
;       IUPPER_LIMITS()
;
; PURPOSE:
;       Compute upper limits on emission line(s).
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;       flux     - spectrum [erg/s/cm2/Angstrom]
;       linefit  - linefit structure from IFITSPEC for every line
;                  requiring an upper limit calculation
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;       linefit  - (modified)
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;       We assume that the line widths (LINESIGMA) of the undetected
;       emission lines have been tied to 
;
;
; PROCEDURES USED:
;
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 June 11, U of A, excised from IFITSPEC() 
;-

function iupper_limits, flux, linefit, snrcut=snrcut

    nline = n_elements(linefit)
    if nline eq 0L then begin
       print, 'Syntax - linefit = iupper_limits(flux,linefit,[snrcut=])'
       return, -1L
    endif

    light = 2.99792458D5            ; speed of light [km/s]
    zline = djs_mean(linefit.linez) ; mean redshift
    
    for j = 0L, nline-1L do begin
       
; compute the Gaussian width of the undetected emission line.  convert
; from [km/s] to [Angstrom]

       sigma = linefit[j].linewave*(1+zline)*linefit[j].linesigma/light ; [Angstrom]

; compute the continuum and the continuum error
       
       djs_iterstat, flux, sigrej=2.5, maxiter=50, median=linecontlevel, sigma=linecontlevel_err
       linearea = snrcut*linecontlevel_err*sqrt(2.0*!pi)*sigma
       
       linefit[j].linesigma_err = -3.0
       linefit[j].linearea = linearea
       linefit[j].linearea_err = -3.0
       linefit[j].linebox = 0.0
       linefit[j].linebox_err = -3.0
       linefit[j].linecontlevel = linecontlevel
       linefit[j].linecontlevel_err = -3.0
       linefit[j].linechi2 = -3.0

;      print, linefit[j].linename, linecontlevel, linecontlevel_err
       
    endfor

return, linefit
end
