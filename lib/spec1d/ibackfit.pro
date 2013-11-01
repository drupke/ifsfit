;+
; NAME:
;	IBACKFIT()
;
; PURPOSE:
;	Fit a linear combination of spectral templates to a galaxy
;	spectrum.
;
; CALLING SEQUENCE:
;       backfit = ibackfit(flux,wave,starflux,[invvar=,mask=],$
;          [nback=,nmonte=,bsorder=,redindx=],/nodust,_extra=extra)
;
; INPUTS:
;	flux     - input spectrum [NPIX]
;	wave     - corresponding rest wavelength vector [NPIX]
;	starflux - stellar template spectra [NPIX,NSTAR] with the same
;                  wavelength spacing as the data
;
; OPTIONAL INPUTS:
;	invvar   - inverse variance spectrum [NPIX]
;	nback    - number of background polynomial terms to append to
;                  STARFLUX (default: 0)
;       nmonte   - number of Monte Carlo realizations to get the error
;                  in the continuum fit (default: 30)
;	bsorder  - order of the quartic b-spline to fit to the
;                  residuals of the background-subtracted data
;                  (default: 0)
;	redindx  - integer array [NSTAR]; stars with the same REDINDX
;                  are constrained to have the same reddening
;                  (default: one reddening for all templates)
;       extra    - keywords for K_LAMBDA() to specify the reddening 
;                  curve (default Charlot & Fall 2001)
;	
; KEYWORD PARAMETERS:
;	nodust   - do not incorporate dust reddening into the fit
;                  (ignore REDINDX)
;
; OUTPUTS:
;	backfit      - output data structure with the following fields 
;          continuum       - best-fitting continuum (best fit + RESIDFIT) 
;          continuum_sigma - error in CONTINUUM
;          chi2            - reduced chi2 of the fit, excluding masked
;                            pixels
;          nstar           - number of stellar templates
;          nback           - number of polynomial background terms 
;          maskpix         - masked pixels that were masked - see
;                            EMISSION_MASK() 
;          nmonte          - number of Monte Carlo realizations
;          ebv             - best-fitting color excess [NSTAR]
;          redindx         - INPUT
;          starcoeff       - best-fitting coefficients [NSTAR+NBACK] 
;
; OPTIONAL OUTPUTS:
;
; PROCEDURE:
;       We use a non-negative least squares matrix solution to
;       determine an initial guess for the coefficients of the fit.
;       We then improve upon these coefficients using MPFIT() and an
;       optional extinction curve.
;
; COMMENTS:
;	The amount of dust reddening is constrained to be between 0
;	and 1.0 in E(B-V).  Some of the MPFIT code is based on code in
;	D. Schlegel's LINEBACKFIT routine, especially constraining the
;	reddening parameters.
;
;       The reddening curve is normalized at 5500 Angstrom.
; 
; EXAMPLE:
;
; INTERNAL SUPPORT ROUTINES:
;       BVLSBACKFIT(), BACKMODEL(), MPBACKFIT()
;
; PROCEDURES USED:
;	POLY_ARRAY(), BSPLINE_ITERFIT(), MPFITFUN(), EMISSION_MASK(), 
;	SPLOG, K_LAMBDA(), BVLS
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2002 May 14-15, U of A - major rewrite of
;          existing code
;       jm03jun9uofa - various bug fixes
;-

function bvlsbackfit, fitflux, fitsigma, fitstarflux, nstar=nstar, nmonte=nmonte, quiet=quiet
; ---------------------------------------------------------------------------    
; BVLS continuum fitting
; ---------------------------------------------------------------------------    

    npix = n_elements(fitflux)
    
    bnd = fltarr(2,nstar) ; boundaries
    bnd[1,*] = 1.0

    if (~ keyword_set(quiet)) then $
      if (nmonte gt 0L) then splog, 'Monte Carlo fitting the continuum '+strn(nmonte)+' times with BVLS.' $
      else splog, 'Fitting the continuum with BVLS.'
;   t1 = systime(/seconds)
    for j = 0L, nmonte do begin
        
       Abvls = fitstarflux
       if j eq 0L then Bbvls = fitflux else Bbvls = fitflux + randomn(seed,npix)*fitsigma ; perturbed flux
       
       bvls, Abvls, Bbvls, bnd, coeff, rnorm, nsetp, w, $
         index, ierr, eps=eps, itmax=itmax, iter=iter

       if j eq 0L then eigencoeff = coeff else eigencoeff = [ [eigencoeff], [coeff] ]
       
    endfor 
;   print, systime(/seconds)-t1
    
return, eigencoeff
end

function backmodel, x, p, wave=wave, redcurve=redcurve, nstar=nstar, $
  nback=nback, nodust=nodust
; jm02may15uofa
; we assume that there are as many reddening terms as there are
; stellar templates.  ie, nred = nstar = n_elements(p)-nback.

; normalize the reddening curve at 5500 Angstrom

    ntemplate = nstar+nback

    if keyword_set(nodust) then begin

       newflux = x

    endif else begin

; redden each stellar template by the normalized reddening curve

       newflux = x[*,0L:nstar-1L]*0.0
       for k = 0L, nstar-1L do newflux[*,k] = x[*,k] * $
         10.0^(0.4*p[ntemplate+k]*(redcurve-interpol(redcurve,wave,5500.0)))

;      plotflux = newflux       ; x
;      djs_plot, wave, plotflux[*,0], xsty=3, ysty=3, ps=4
;      for k = 1L, nstar-1L do djs_oplot, wave, plotflux[*,k], ps=4
    
    endelse

; make the model

    model = newflux # p[0L:ntemplate-1L]

return, model
end

function mpbackfit, fitflux, fitwave, fitsigma, fitstarflux, $
                    redindx=redindx, redcurve=redcurve, $
                    coeffguess=coeffguess, nstar=nstar, nback=nback, $
                    nmonte=nmonte, nodust=nodust, quiet=quiet
; ---------------------------------------------------------------------------    
; MPFIT continuum fitting    
; ---------------------------------------------------------------------------    

    npix = n_elements(fitflux)
    
    if keyword_set(nodust) then ntemplate = nstar+nback else ntemplate = 2*nstar+nback
    
; initialize the fitting parameters.  constrain the stellar
; contribution to be positive and also constrain the range of the
; reddening

    parinfo = {value: 0.0,     $
               fixed: 0L,      $
               tied: '',       $
               limited: [0,0], $
               limits: [0.0,0.0]}
    parinfo = replicate(parinfo,ntemplate)
    
    parinfo[0:nstar+nback-1L].limited[0] = 1L    ; stellar contribution must be positive

    if n_elements(coeffguess) eq 0L then $
      parinfo[0:nstar+nback-1L].value = 1.0/nstar else $  ; initial guess (equal contributions)    
      parinfo[0:nstar+nback-1L].value = coeffguess        ; input coefficients
    
    if not keyword_set(nodust) then begin 

       parinfo[nstar+nback:ntemplate-1L].limited = [1L,1L] ; constrain the reddening range
       parinfo[nstar+nback:ntemplate-1L].limits = [0.0,1.5]
       parinfo[nstar+nback:ntemplate-1L].value = 0.1       ; initial guess

; tie reddening values to one another if requested
    
       allindx = redindx[uniq(redindx,sort(redindx))]
       for i = 0L, n_elements(allindx)-1L do begin
          ii = where(redindx eq allindx[i],count)
          if count gt 1L then begin
             for j = 1L, count-1L do begin
                parinfo[nstar+nback+ii[j]].tied = string(nstar+nback+ii[0],format='("P[",i,"]")')
             endfor
          endif
       endfor

    endif 
       
; reddening

    functargs = {wave: fitwave, redcurve: redcurve, nstar: nstar, $
      nback: nback, nodust: keyword_set(nodust)}

; Monte Carlo the fit to find the error in the continuum

    if (~ keyword_set(quiet)) then $
      if nmonte gt 0L then splog, 'Monte Carlo fitting the continuum '+strn(nmonte)+' times with MPFIT.' $
      else splog, 'Fitting the continuum with MPFIT.'

;   t1 = systime(/seconds)
    oflux = fitflux 
    for j = 0L, nmonte do begin

       if j gt 0L then fitflux = oflux + randomn(seed,npix)*fitsigma ; perturbed flux

       coeff = mpfitfun('backmodel',fitstarflux,fitflux,fitsigma,$
         parinfo=parinfo,functargs=functargs,covar=covar,perror=perror,yfit=yfit1,$
         nfev=nfev,niter=niter,status=status,/quiet)

       if j eq 0L then eigencoeff = coeff else eigencoeff = [ [eigencoeff], [coeff] ]
       
    endfor 
;   print, systime(/seconds)-t1

return, eigencoeff
end

function ibackfit, flux, wave, starflux, invvar=invvar, mask=mask, $
  nback=nback, nmonte=nmonte, bsorder=bsorder, redindx=redindx, $
  nodust=nodust, quiet=quiet, _extra=extra

    if (~ keyword_set(quiet)) then quiet=0

    light = 2.99792458D5

    if n_elements(flux) eq 0L then begin
       print, 'Syntax - backfit = ibackfit(flux,wave,starflux,[invvar=,mask=],$'
       print, '   [nback=,nmonte=,bsorder=,redindx=],/nodust,_extra=extra)'
       return, -1
    endif 
       
    ndim = size(flux,/n_dimension) 
    if ndim ne 1L then message, 'FLUX vector is not one dimensional.' 
    fsize = size(flux,/dimension)
    npix = fsize[0]

    if (n_elements(wave) eq 0L) or (n_elements(wave) ne npix) then $
      message, 'Either WAVE is not defined or FLUX and WAVE have incompatible dimensions.'
    
    if n_elements(invvar) eq 0L then invvar = replicate(1.0,npix) else $
      if n_elements(invvar) ne npix then message, 'FLUX and INVVAR have incompatible dimensions.'

    if n_elements(mask) eq 0L then mask = replicate(1,npix) else $
      if n_elements(mask) ne npix then message, 'FLUX and MASK have incompatible dimensions.'

    if n_elements(nback) eq 0L then nback = 0L
    if n_elements(nmonte) eq 0L then nmonte = 30L
    if n_elements(bsorder) eq 0L then bsorder = 0.0 ; 30.0

    stardim = size(starflux,/n_dimension)
    starsize = size(starflux,/dimension)
    if starsize[0] ne npix then message, 'FLUX and STARFLUX have incompatible dimensions.'
    if stardim eq 1L then nstar = 1L else nstar = starsize[1]

    if n_elements(redindx) eq 0L then redindx = lonarr(nstar) else $
      if n_elements(redindx) ne nstar then message, 'STARFLUX and REDINDX have incompatible dimensions.'
    
; add polynomial templates to the stellar templates

    if nback gt 0L then begin

       polyflux = poly_array(npix,nback)
       starflux = [ [starflux],[polyflux] ]       

    endif

    if keyword_set(nodust) then ntemplate = nstar+nback else $
      ntemplate = $
      nstar + $                 ; number of stellar templates
      nback + $                 ; number of polynomial terms
      nstar                     ; number of reddening terms (constrained by REDINDX)

    fitinvvar = invvar*mask
    good = where(fitinvvar gt float(0),ngood,comp=maskpix,ncomp=nmask)
    
    fitflux = flux[good]
    fitwave = wave[good]
    fitsigma = 1.0/sqrt(fitinvvar[good])

    fitstarflux = starflux[good,*]
    
; fit the continuum.  find the best-fitting solution with BVLS and, in
; general, Monte Carlo this solution to obtain the error in the
; coefficients.  then input this guess into MPFIT to refine the
; solution and to solve for the reddening
    
; BVLS

    bvlscoeff = bvlsbackfit(fitflux,fitsigma,fitstarflux,nstar=(nstar+nback),nmonte=nmonte,quiet=quiet)
    if nmonte gt 0L then coeffguess = bvlscoeff[*,0] else coeffguess = bvlscoeff

; MPFIT    
    
    if keyword_set(nodust) then redcurve = fltarr(npix) else $ ; reddening curve
      redcurve = k_lambda(wave,/calzetti,_extra=extra)

    mpfitcoeff = mpbackfit(fitflux,fitwave,fitsigma,fitstarflux,redindx=redindx,$
      redcurve=redcurve[good],coeffguess=coeffguess,nstar=nstar,nback=nback,$
      nmonte=0,nodust=nodust,quiet=quiet)

    continuum = backmodel(starflux,mpfitcoeff,wave=wave,redcurve=redcurve,nstar=nstar,$
      nback=nback,nodust=nodust)

    if not keyword_set(nodust) then $
      reddening = mpfitcoeff[nstar+nback:ntemplate-1L] else $ ; E(B-V)
      reddening = fltarr(nstar)
    starcoeff = mpfitcoeff[0L:nstar+nback-1L]                 ; final coefficients

    percent = 100.0*starcoeff/total(starcoeff) ; percentage contribution
    
; evaluate the error in the continuum if NMONTE > 0

    if nmonte gt 0L then begin
    
       continuum_array = continuum*0.0 # (fltarr(nmonte)+1)
       for k = 0L, nmonte-1L do continuum_array[*,k] = backmodel(starflux,[bvlscoeff[*,k+1],reddening],$
         wave=wave,redcurve=redcurve,nstar=nstar,nback=nback,nodust=nodust)
       continuum_sigma = continuum*0.0

       for i = 0L, npix-1L do continuum_sigma[i] = stddev(continuum_array[i,*])
       
    endif else continuum_sigma = continuum*0.0
    
;   plot, wave, flux, xsty=3, ysty=3, ps=10
;   for j = 0L, nmask-1L do plots, wave[maskpix[j]], flux[maskpix[j]], ps=4, color=djs_icolor('red')
;   djs_oplot, wave, continuum, color='green', thick=2.0
;   plot, wave, flux-continuum, xsty=3, ysty=3, ps=10
    
; fit the residuals with a high-order b-spline (probably never appropriate)

    residuals = flux-continuum
    residvar = residuals*0.0+1.0

    if bsorder gt 0 then begin

       sset = bspline_iterfit(wave,residuals,invvar=fitinvvar,bkpt=0,$
         everyn=npix/bsorder,yfit=residfit,nord=nord,lower=1.0,upper=1.0)

    endif else residfit = continuum*0.0

    continuum = continuum + residfit

    chi2 = total(fitinvvar[good]*(flux[good]-continuum[good])^2.0) / $
      (ngood-ntemplate-1L)
    if (~ keyword_set(quiet)) then begin
      splog, 'Reduced chi2 of the continuum fit = ', chi2
      if (~ keyword_set(nodust)) then splog, 'Best-fitting E(B-V) = ', reddening[uniq(reddening)]
    endif

; NB: we are not returning RESIDFIT and BSORDER
    
    backfit = {$
      continuum:       float(continuum), $
      continuum_sigma: float(continuum_sigma), $
      continuum_chi2:  float(chi2), $
      nstar:           long(nstar), $
      nback:           long(nback), $
      nmonte:          long(nmonte), $
      maskpix:         long(maskpix), $
      ebv:             float(reddening), $
      redindx:         long(redindx), $
      starcoeff:       float(starcoeff)}
    
;  djs_plot, flux, ps=10, xsty=3, ysty=3
;  ploterror, flux, 1.0/sqrt(invvar), ps=10, xsty=3, ysty=3
;  djs_oplot, continuum, color='red', ps=10

return, backfit
end    
