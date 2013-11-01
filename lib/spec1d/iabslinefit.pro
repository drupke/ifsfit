;+
; NAME:
;       IABSLINEFIT()
;
; PURPOSE:
;       Fit absorption line fluxes and equivalent widths in a galaxy
;       spectrum.
;
; CALLING SEQUENCE:
;       absfit = iabslinefit(wave,flux,ferr,abslines,$
;          [lineres=],/doplot)
;
; INPUTS:
;       wave     - wavelength vector [NPIX]
;       flux     - spectrum flux in erg/s/cm2/Angstrom [NPIX]
;       ferr     - corresponding error spectrum [NPIX]
;       abslines - fit the absorption lines specified in this
;                  structure [NLINE, see READ_LINEPARS()]
;
; OPTIONAL INPUTS:
;       lineres    - Gaussian sigma spectral resolution line width in 
;                    km/s [NLINE]
;
; KEYWORD PARAMETERS:
;       doplot   - generate a plot of the fit to each absorption line 
;
; OUTPUTS:
;       absfit   - absorption-line fitting results [see ILINEFIT()]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; PROCEDURES USED:
;       ICREATE_LINESTRUCT(), MEASURE_LINEEW()
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 July 4, U of A, excised from IBACKFIT()
;       jm03jul30uofa - verify that the absorption lines are within
;                       the given wavelength range
;-

function iabslinefit, wave, flux, ferr, abslines, lineres=lineres, doplot=doplot

    nwave = n_elements(wave)
    nflux = n_elements(flux)
    nferr = n_elements(ferr)
    nline = n_elements(abslines)

    if (nwave eq 0L) or (nflux eq 0L) or (nferr eq 0L) then begin
       print, 'Syntax - absfit = iabslinefit(wave,flux,ferr,abslines,$'
       print, '   [lineres=],/doplot)'
       return, -1L
    endif

    if (nwave ne nflux) or (nwave ne nferr) then begin
       print, 'Dimensions of WAVE, FLUX, and FERR do not agree.'
       return, -1L
    endif
    
    if n_elements(lineres) ne 0L then begin
       if n_elements(lineres) ne nline then begin
          print, 'Dimensions of LINERES and ABSLINES do not agree.'
          return, -1L
       endif
    endif

    absfit = icreate_linestruct(nline)
    absfit.linename = abslines.line
    absfit.linewave = abslines.wave
    absfit.linesigma = 200.0      ; [km/s]
    absfit.linearea_err = 100.0   ; this is to trick MEASURE_LINEEW()
    absfit.line_blend = abslines.blend

    bad = where((absfit.linewave gt max(wave)) or (absfit.linewave lt min(wave)),$
      nbad,comp=good,ncomp=ngood)

    if nbad ne 0L then begin

       splog, 'The absorption line(s) '+absfit[bad].linename+$
          ' are outside the data wavelength range.'

       if n_elements(lineres) ne 0L then absfit[bad].linesigma = lineres[bad]
       absfit[bad].linesigma_err     = -2.0
       absfit[bad].linearea_err      = -2.0
       absfit[bad].linebox_err       = -2.0
       absfit[bad].linecontlevel_err = -2.0
       absfit[bad].lineew_area_err   = -2.0
       absfit[bad].lineew_box_err    = -2.0
       absfit[bad].linechi2          = -2.0
       
    endif

    if ngood ne 0L then begin

       splog, 'Fitting '+strn(ngood)+' absorption lines.'

       fitthem = absfit[good]

       fitthem = measure_lineew(fitthem,wave,flux,ferr,snrcut=0.0,$
         lineres=lineres,sigmax=700.0,/absorption,/overwrite,$
         doplot=doplot)

       absfit[good] = fitthem

    endif
    
    good = where(absfit.linearea_err gt 0.0,ngood,comp=bad,ncomp=nbad)
    if nbad ne 0L then splog, 'Absorption line(s) '+$
      absfit[bad].linename+' not successfully fitted.'

return, reform(absfit)
end    
