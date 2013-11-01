;+
; NAME:
;       ILINEFIT()
;
; PURPOSE:
;       Fit an emission-line spectrum.
;
; CALLING SEQUENCE:
;
; INPUTS:
;	flux       - emission-line spectrum to fit [NPIX] 
;       wave       - corresponding wavelength vector [NPIX] 
;       linelambda - rest-frame wavelengths of emission lines to fit
;                    in Angstroms [NLINE] 
;       lineres    - Gaussian sigma spectral resolution line width in
;                    km/s [NLINE] 
;
; OPTIONAL INPUTS:
;       invvar     - inverse variance spectrum [NPIX].  if not
;                    provided then the emission line flux errors will
;                    be incorrect and the fitting may even fail
;       linename   - String name(s) of line(s) to be copied to output
;                    structure.
;       zindex     - Lines with the same ZINDEX are constrained to
;                    have the same redshift; default to a unique list
;                    of indices [NLINE].
;       windex     - Lines with the same WINDEX are constrained to
;                    have the same width [NLINE]; default to a unique
;                    list of indices.
;       findex     - Lines with the same FINDEX are constrained to
;                    have the same flux ratio as input in FVALUE;
;                    defult to a unique list of indices [NLINE].
;       fvalue     - If FINDEX is specified, then constrain the lines
;                    to have the same flux ratios as in FVALUE
;                    [NLINE]; default to all values of 1.  These
;                    values are also used as the initial guesses for
;                    the line strengths.
;       sigmin     - minimum emission-line width [km/s] (default 0)
;       sigmax     - maximum emission-line width [km/s] (default 400)
;       background - Background vector(s) to fit simultaneously with
;                    the lines, where the scaling of each vector is
;                    fit [NPIX,NBACK]. The redshift of these
;                    background vectors is not fit, but rather we
;                    maintain the one-to-one correspondence of each
;                    BACKGROUND pixel to each FLUX pixel.  The initial
;                    guess for the scaling of each background vector
;                    is unity.
;       sigguess   - Initial guess for *intrinsic* sigmas of all lines in log-10
;                    Angstroms (scalar or vector with one entry per
;                    line); default to 1.5d-4 (105 km/sec).
;       zguess     - redshift guess (needs to be pretty close, or set
;                    to zero and pass the rest wavelength vector).
;                    can be a scalar or a vector with one entry per
;                    emission line
;
; KEYWORD PARAMETERS:
;       fixzguess  - do not allow the emission-line redshifts to vary
;                    by more than +/- 5% relative to ZGUESS
;       allowneg   - allow negative emission lines (to fit absorption
;                    lines).  the default is to constrain the lines to
;                    be positive
;
; OUTPUTS:
;       linestruct - Output structure with result of line fits [NLINE] 
;          LINENAME        - String name of line copied from input
;                            parameter by the same name, or else
;                            constructed from the rounded-down wavelength
;                            of each line in air
;          LINEWAVE        - rest-frame wavelength [Angstrom]
;                            (copy of LINELAMBDA)
;          LINEZ           - Redshift [dimensionless]
;          LINEZ_ERR       - Error in above
;          LINESIGMA       - Sigma of gaussian [km/s]
;          LINESIGMA_ERR   - Error in above
;          LINEAREA        - Area of line [flux-units * Ang]
;          LINEAREA_ERR    - Error in above
;          LINECONTLEVEL   - Continuum level at line center [flux-units];
;                            if the line center is outside the wavelength
;                            range, then return the nearest value (either
;                            the first or last value)
;          LINECONTLEVEL_ERR - Error in above, or -1L if the line center
;                            is outside the wavelength range.
;          LINENPIX        - Number of pixels within +/- 3 sigma of
;                            the line center that have INVVAR > 0.
;          LINEDOF         - LINENPIX minus the number of terms fit
;                            for that line, which could be fractional
;                            (if one parameter is fixed between N lines,
;                            then we say only 1/N-th of that parameter
;                            is fit in each of those lines).  This can
;                            be zero or negative.
;          LINECHI2       -  Chi^2 for all points within +/- 3 sigma of
;                            the line center; -1L if no such points.
;
; OPTIONAL OUTPUTS:
;       specfit - Fit spectrum including lines and background terms [NPIX].
;       bfit    - Fit spectrum including only background terms [NPIX].
;       bterms  - Coefficients for background terms [NBACK].
;
; COMMENTS:
;	This program assumes that a pure emission-line spectrum is
;	being fitted.  
;
;       All wavelengths must be in Angstroms.
;
;       If a line was dropped from the fit (for example, no points to
;       fit), then set the LINEAREA to 0 and the LINEAREA_ERR to -1L.
;
;       Also, if LINENPIX=0 for a line, then remove that line from the
;       fit. 
;
;       Possible bug:  Do not use lines with no points to fit in the
;       computation of degrees of freedom for other lines.
;
; EXAMPLES:
;
; INTERNAL ROUTINES:
;       ONEGAUSS(), MANYGAUSS()
;
; PROCEDURES USED:
;       ICREATE_LINESTRUCT()
;
; MODIFICATION HISTORY:
;       05-Feb-2002  Written by D. Schlegel, Princeton
;       J. Moustakas, 2002 March 15, U of A
;-

function onegauss, xval, pp, sigmares=sigmares
; the sigma line-width is comprised of the fixed spectral resolution
; width and the variable intrinsic line width
    
    sigma_squared = sigmares^2.0 + pp[2]^2.0 ; total Gaussian sigma width [log-Angstrom]
    sigma = sqrt(sigma_squared)
    
    term1 = exp( - (xval - pp[1])^2 / (2. * sigma_squared ) )
;    yval = pp[0] * term1 / (sqrt(2.*!pi) * sigma)
    yval = pp[0] * term1 

;   term1 = exp( - (xval - pp[1])^2 / (2. * pp[2]^2) )
;   yval = pp[0] * term1 / (sqrt(2.*!pi) * pp[2])

return, yval
end

function manygauss, xindx, pp, nline=nline, nback=nback, loglam=loglam, $
  sigmares=sigmares, background=background

    yval = 0.d

    for iline=0, nline-1 do $
      yval = yval + onegauss(loglam[xindx], pp[iline*3:iline*3+2], sigmares=sigmares[iline])
    for iback=0, nback-1 do $
      yval = yval + background[xindx,iback] * pp[nline*3+iback]

return, yval
end

function ilinefit, flux, wave, linelambda, lineres, invvar=invvar, $
  linename=linename, zindex=zindex, windex=windex, findex=findex, $
  fvalue=fvalue, sigmin=sigmin, sigmax=sigmax, background=background, $
  zguess=zguess1, sigguess=sigguess1, specfit=specfit, bfit=bfit, $
  bterms=bterms, fixzguess=fixzguess, allowneg=allowneg, $
  perror = perror, lfit = lfit

    light = 2.99792458D5 ; speed of light [km/s]

    npix = n_elements(flux)
    nwave = n_elements(wave)
    nline = n_elements(linelambda)

    if (npix eq 0L) or (nwave eq 0L) or (nline eq 0L) then begin
       print, 'Syntax - linestruct = ilinefit(flux,wave,linelambda,lineres,$'
       print, '   [invvar=,linename=,zindex=,windex=,findex=,fvalue=,sigmin=,$'
       print, '   sigmax=,background=,zguess=,sigguess=,specfit=,bfit=,bterms=],$'
       print, '   /fixzguess,/allowneg)'
       return, -1
    endif

    if (nwave ne npix) then begin
       splog, 'WAVE and FLUX do not have the same number of elements.'
       return, -1L
    endif

    nivar = n_elements(invvar)
    if (nivar eq 0L) then invvar = make_array(npix,value=1.0,/float) else begin
       if nivar ne npix then begin
          splog, 'FLUX and INVVAR do not have the same number of elements.'
          return, -1L
       endif
    endelse

    ndim = size(background,/n_dimension)
    if (ndim eq 0L) then begin
       background = flux*0.0
       nback = 0L
    endif else begin
       dims = size(background,/dimension)
       if (ndim eq 1) then begin
          nback = 1L
          if n_elements(background) ne npix then begin
             splog, 'FLUX and BACKGROUND do not have the same number of elements.'
             return, -1L
          endif
       endif
       if (ndim eq 2) then begin
          nback = dims[1]
          if n_elements(background[*,0]) ne npix then begin
             splog, 'FLUX and BACKGROUND[*,0] do not have the same number of elements.'
             return, -1L
          endif
       endif
    endelse 

; convert the spectral resolution from km/s to log-Angstrom units

    nlineres = n_elements(lineres)
    if (nlineres eq 0L) then begin
       splog, 'LINERES must be provided.'
       return, -1L
    endif

    if (nlineres ne nline) then begin
       splog, 'LINERES must have the same number of elements as LINELAMBDA.'
       return, -1L
    endif

    sigmares = lineres / light / alog(10.0)

; set defaults
    
    if (NOT keyword_set(zindex)) then zindex = lindgen(nline)
    if (NOT keyword_set(windex)) then windex = lindgen(nline)
    if (NOT keyword_set(findex)) then findex = lindgen(nline)
    if (NOT keyword_set(fvalue)) then fvalue = fltarr(nline) + 1.0

    if (keyword_set(zguess1)) then begin
       if (n_elements(zguess1) EQ 1) then begin
          zguess = replicate(zguess1[0], nline)
       endif else if (n_elements(zguess1) EQ nline) then begin
          zguess = zguess1
       endif else begin
          splog, 'Wrong number of elements for ZGUESS.'
          return, -1L
       endelse
    endif else begin
       zguess = replicate(0D,nline)
    endelse

    if (keyword_set(sigguess1)) then begin
       if (n_elements(sigguess1) EQ 1) then begin
          sigguess = replicate(sigguess1[0], nline)
       endif else if (n_elements(sigguess1) EQ nline) then begin
          sigguess = sigguess1
       endif else begin
          splog, 'Wrong number of elements for SIGGUESS.'
          return, -1L
       endelse
    endif else begin
       sigguess = replicate(105.0/light/alog(10.0),nline) ; default sigma width [log-Angstrom]
    endelse 
    
    loglam = alog10(wave) ; log-Angstrom

; initialize the LINESTRUCT output structure

    linestruct = icreate_linestruct(nline)

    linestruct.linewave = linelambda
    if (n_elements(linename) eq nline) then linestruct.linename = linename else $
      linestruct.linename = strtrim(string(long(linelambda),format='(I0)'),2)

; initialize the structures to be passed to the fitting routine

    parinfo = replicate({value:0.D, fixed:0, limited:[0,0], tied:'', $
      limits:[0.D,0]}, nline*3+nback)
    functargs = {nline: nline, nback: nback, loglam: loglam, $
      sigmares: sigmares, background: background}

 ; set the initial guesses of the fitting parameters

    if n_elements(sigmin) eq 0L then sigmin = 0.0D  ; constrain the sigma widths [km/s]
    if n_elements(sigmax) eq 0L then sigmax = 400.0
    
    for iline=0, nline-1 do begin

       parinfo[0+iline*3].value = fvalue[iline]
       parinfo[1+iline*3].value = alog10(linelambda[iline]*(1+zguess[iline]))
       parinfo[2+iline*3].value = sigguess[iline]

; unless ALLOWNEG has been set, constrain the emission lines to be
; positive

       if (n_elements(allowneg) eq 0L) then begin
;          parinfo[0+iline*3].limited[0] = 1B
;          parinfo[0+iline*3].limits[0] = 0.0D
          parinfo[0+iline*3].limited = [1B, 1B]
          parinfo[0+iline*3].limits  = [0.0D, 1d]
     
       endif

; constrain the sigma widths between [SIGMIN,SIGMAX] (log-Angstrom)

       parinfo[2+iline*3].limited = 1B
       parinfo[2+iline*3].limits[0] = sigmin/alog(10.0)/light
       parinfo[2+iline*3].limits[1] = sigmax/alog(10.0)/light

; if the redshift is relatively well known then constrain it such that
; it cannot change by more than 5%

       if keyword_set(fixzguess) then begin
          parinfo[1+iline*3].limited = 1B
          parinfo[1+iline*3].limits[0] = alog10(linelambda[iline]*(1-1.05*zguess[iline]))
          parinfo[1+iline*3].limits[1] = alog10(linelambda[iline]*(1+1.05*zguess[iline]))
       endif

    endfor

; Set the initial guess for the initial background level to be the
; median flux

    if (nback ne 0L) then parinfo[nline*3+0].value = median(flux)
    
;;   ; Set the initial guess for the scaling of each background vector to unity.
;;   for iback=0, nback-1 do begin
;;      parinfo[nline*3+iback].value = 1.0
;;   endfor

; Make a list of the number of fitting terms per line.  If a parameter
; is constrained between N lines, then we say each of those lines is
; only fitting 1/N-th of that parameter

    nfitterms = fltarr(nline)

; Apply constraints to peak flux values (not integrated flux in a line)

    allindex = findex[ uniq(findex,sort(findex)) ]
    for iall=0, n_elements(allindex)-1 do begin
       ii = where(findex EQ allindex[iall], ct)
       nfitterms[ii] = nfitterms[ii] + 1.0/ct
       if (ct GT 1) then begin
          for jj=1, ct-1 do begin
             fratio = fvalue[ii[jj]] / fvalue[ii[0]]
             parinfo[0+ii[jj]*3].tied = $
               string(fratio, 0+ii[0]*3, $
               format='(e12.5," * P(",i,")")')
          endfor
       endif
    endfor

; Apply constraints to couple redshifts

    allindex = zindex[ uniq(zindex,sort(zindex)) ]
    for iall=0, n_elements(allindex)-1 do begin
       ii = where(zindex EQ allindex[iall], ct)
       nfitterms[ii] = nfitterms[ii] + 1.0/ct
       if (ct GT 1) then begin
          for jj=1, ct-1 do begin
             lamshift = alog10(linelambda[ii[jj]] / linelambda[ii[0]])
             parinfo[1+ii[jj]*3].tied = $
               string(lamshift, 1+ii[0]*3, $
               format='(e12.5," + P(",i,")")')
          endfor
       endif
    endfor

; Apply constraints to couple widths

    allindex = windex[ uniq(windex,sort(windex)) ]
    for iall=0, n_elements(allindex)-1 do begin
       ii = where(windex EQ allindex[iall], ct)
       nfitterms[ii] = nfitterms[ii] + 1.0/ct
       if (ct GT 1) then begin
          for jj=1, ct-1 do begin
             parinfo[2+ii[jj]*3].tied = $
               string(2+ii[0]*3, format='("P(",i,")")')
          endfor
       endif
    endfor

; Do the fit!
    
    specfit = fltarr(npix)

    igood = where(invvar GT 0, ngood)
    status = 0
    if (ngood GT 0) then begin
       lfit = mpfitfun('manygauss', igood, double(flux[igood]), $
         1./sqrt(invvar[igood]), parinfo=parinfo, $
         covar=covar, perror=perror, yfit=specfit1, functargs=functargs, $
         nfev=nfev, niter=niter, status=status, /quiet)
;     splog, 'MPFIT number of function evaluations=', nfev
;     splog, 'MPFIT number of iterations=', niter
;     splog, 'MPFIT exit status=', status
       if (status EQ 5) then $
         splog, 'Warning: Maximum number of iterations reached: ', niter
       specfit[igood] = specfit1
    endif


    if (ngood EQ 0 OR status EQ 0) then begin
       splog, 'Too few points to fit ', ngood
       nparam = 3 * nline + nback
       lfit = fltarr(nparam)
       perror = fltarr(nparam)
       covar = fltarr(nparam,nparam)
    endif

; convert -0.0 to 0.0

    negzero = where_negzero(perror,negcount)
    if negcount ne 0L then perror[negzero] = 0.0
    
; For parameters that are fixed, assign them the same errors as those
; parameters to which they are tied.  For the flux area, scale those
; errors according to the ratio of the flux levels.

    allindex = findex[ uniq(findex,sort(findex)) ]
    for iall=0, n_elements(allindex)-1 do begin
       ii = where(findex EQ allindex[iall], ct)
       if (ct GT 1) then perror[ii*3+0] = perror[ii[0]*3+0] $
         * fvalue[ii] / fvalue[ii[0]]
    endfor

    allindex = zindex[ uniq(zindex,sort(zindex)) ]
    for iall=0, n_elements(allindex)-1 do begin
       ii = where(zindex EQ allindex[iall], ct)
       if (ct GT 1) then perror[ii*3+1] = perror[ii[0]*3+1]
    endfor

    allindex = windex[ uniq(windex,sort(windex)) ]
    for iall=0, n_elements(allindex)-1 do begin
       ii = where(windex EQ allindex[iall], ct)
       if (ct GT 1) then perror[ii*3+2] = perror[ii[0]*3+2]
    endfor

; Construct the line-measure outputs (and their errors)

    linestruct.linearea      = lfit[lindgen(nline)*3+0] $
      * alog(10.) * 10.^lfit[lindgen(nline)*3+1]
    linestruct.linez         =  (10.^lfit[lindgen(nline)*3+1] / linelambda - 1) $
      * (lfit[lindgen(nline)*3+1] GT 0) ; Set to to zero if the parameter is zero
    linestruct.linesigma     = abs(lfit[lindgen(nline)*3+2]) * alog(10.) * light
    linestruct.linearea_err  = perror[lindgen(nline)*3+0] $
      * alog(10.) * 10.^lfit[lindgen(nline)*3+1]
    linestruct.linez_err     = perror[lindgen(nline)*3+1] $
      * alog(10.) * (linestruct.linez + 1)
    linestruct.linesigma_err = perror[lindgen(nline)*3+2] * alog(10.) * light

; If a line was dropped from the fit by MPFIT (I'm not sure under what
; circumstances this occurs) then set the LINEAREA to 0 and the
; LINEAREA_ERR to -1L.  we will compute an upper limit on these lines
; in IFITSPEC()

    ibad = where(perror[lindgen(nline)*3+0] LE 0)
    if (ibad[0] NE -1) then begin
       linestruct[ibad].linearea = 0.0
       linestruct[ibad].linearea_err = -1.0
    endif

; assign a linewidth of LINERES to lines whose width is SIGMIN
; (unphysically narrow lines).  we will compute upper limits on these
; lines in IFITSPEC()
    
    ibad = where(perror[lindgen(nline)*3+2] eq SIGMIN)
    if (ibad[0] NE -1) then begin
       lfit[ibad*3+2] = sigmares[ibad]
       linestruct[ibad].linesigma = lineres[ibad]
       linestruct[ibad].linesigma_err = -1.0
    endif

; Find the background levels only

    if (nback EQ 0) then begin
       bterms = 0
       bfit = fltarr(npix)
       berr = bfit*0.0
    endif else begin
       bterms = lfit[nline*3:nline*3+nback-1]
       bcovar = covar[nline*3:nline*3+nback-1,nline*3:nline*3+nback-1]

; The following two methods for evaluating bfit are equivalent:

;      bfit = manygauss(lindgen(npix), bterms, nline=0, nback=nback, $
;      loglam=loglam, background=background)
       bfit = bterms ## background

       berr = fltarr(npix)
       for ipix=0, npix-1 do $
         berr[ipix] = sqrt( (transpose(background[ipix,*]) $
         # bcovar # transpose(background[ipix,*])) )
    endelse

; For each line, determine the background level at the line center and
; the number of pixels and chi^2 of each line fit

;;    logmin = min(loglam)
;;    logmax = max(loglam)
;;    for iline=0, nline-1 do begin
;;       if (lfit[iline*3+1] LT logmin) then begin
;;; Case where the line center is blueward of the entire spectrum.
;;          linestruct[iline].linecontlevel = 0.0      ; bfit[0]
;;          linestruct[iline].linecontlevel_err = -1.0
;;       endif else if (lfit[iline*3+1] GT logmax) then begin
;;; Case where the line center is redward of the entire spectrum.
;;          linestruct[iline].linecontlevel = 0.0      ; bfit[npix-1]
;;          linestruct[iline].linecontlevel_err = -1.0
;;       endif else begin
;;; Select the nearest pixel for evaluating the background
;;; level at this line center.
;;          junk = min(abs(loglam - lfit[iline*3+1]), ipix)
;;          linestruct[iline].linecontlevel = bfit[ipix]
;;          linestruct[iline].linecontlevel_err = berr[ipix]
;;       endelse
;;    endfor

; Find the pixels that are within +/- 3 sigma of the line center.
; Note that if the line is very (unphysically) narrow, it is possible
; to have no lines within this domain.  Reject those fits.  lines that
; floor against SIGMIN or SIGMAX will have very large errors, so don't
; flag those cases as special rejection

    dlam = wave[1]-wave[0] ; [Angstrom/pixel]
    var = 1.0/invvar       ; variance
    
    for iline=0, nline-1 do begin

       indx = where(loglam ge lfit[iline*3+1] - 3*lfit[iline*3+2] $ ; +/- 3-sigma
         and loglam le lfit[iline*3+1] + 3*lfit[iline*3+2] )

       if (indx[0] NE -1) then begin

          linestruct[iline].linenpix = total(invvar[indx] GT 0)
          linestruct[iline].linedof = linestruct[iline].linenpix-nfitterms[iline]

          if (linestruct[iline].linedof GT 0) then begin
             linestruct[iline].linechi2 = total( (flux[indx]-specfit[indx])^2*invvar[indx])
          endif 

; box-car flux measurement and error (jm02apr18uofa); measure the
; total flux contained within +/- 3-sigma of line center.  this
; measurement is meaningless for blended lines

          linestruct[iline].linebox = total(flux[indx]-bfit[indx])*dlam                ; [erg/s/cm2]
          linestruct[iline].linebox_err = sqrt(total(var[indx] + berr[indx]^2.0))*dlam ; [erg/s/cm2]

; Select the nearest pixel for evaluating the background
; level at this line center.

          junk = min(abs(loglam-lfit[iline*3+1]),ipix)
          linestruct[iline].linecontlevel = bfit[ipix]
          linestruct[iline].linecontlevel_err = berr[ipix]

       endif 

; Special-case rejection.  set AREA=0,AREA_ERR=-2 if there are no data
; points within the line-fitting region.  also reject if the edges of
; the line encroach on the endpoints of the wavelength range.  also
; set the fitted redshift to the redshift guess
; (jm03mar10uofa)

;      if (linestruct[iline].linenpix EQ 0) then begin
       if (linestruct[iline].linenpix EQ 0) or (indx[0] eq 0L) or $
         (indx[n_elements(indx)-1L] eq (npix-1L)) then begin
;      if (indx[0] eq 0L) or (indx[n_elements(indx)-1L] eq (npix-1L)) then begin

          linestruct[iline].linez             = zguess[iline]
          linestruct[iline].linez_err         = -2.0
          linestruct[iline].linesigma         =  0.0
          linestruct[iline].linesigma_err     = -2.0
          linestruct[iline].linearea          =  0.0
          linestruct[iline].linearea_err      = -2.0
          linestruct[iline].linebox           =  0.0
          linestruct[iline].linebox_err       = -2.0
          linestruct[iline].lineew_area       =  0.0
          linestruct[iline].lineew_area_err   = -2.0
          linestruct[iline].lineew_box        =  0.0
          linestruct[iline].lineew_box_err    = -2.0
          linestruct[iline].linecontlevel     =  0.0
          linestruct[iline].linecontlevel_err = -2.0
          linestruct[iline].linenpix          =  0
          linestruct[iline].linedof           =  0.0
          linestruct[iline].linechi2          = -2.0

; Set these line-fit coefficients equal to zero for when we
; re-evaluate SPECFIT.

          lfit[iline*3+0] = 0
          
       endif

    endfor

; Re-evaluate such that we get the functional fit at the rejected
; wavelengths too.  This also re-evaluates the fit for lines that have
; been rejected and removed due to no valid data points within the
; line-fitting region.

    if (arg_present(specfit)) then specfit = manygauss(lindgen(npix),lfit,_extra=functargs)

return, linestruct
end
