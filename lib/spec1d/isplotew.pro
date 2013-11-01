function isplotew, wave, flux, ferr, linewave, linename, nmonte=nmonte, $
  boxwidth=boxwidth, doplot=doplot, silent=silent
    
; jm02jul9uofa
; measure the equivalent width of a line similarly to splot (fitting
; the continuum locally) and estimate the error in the continuum using
; Monte Carlo.  

    nline = n_elements(linewave)
    if nline eq 0L then begin
       print, 'Syntax - lineew = isplotew(wave,flux,ferr,linewave,linename,$'
       print, '           [nmonte=],doplot=doplot)'
       return, -1L
    endif

    if nline gt 1L then begin

       if (n_elements(boxwidth) eq 1L) then boxwidth = replicate(boxwidth,nline)
       if (n_elements(boxwidth) eq 0L) then boxwidth = replicate(35.0,nline)
       
       for i = 0L, nline-1L do begin

          lineew1 = isplotew(wave,flux,ferr,linewave[i],linename[i],nmonte=nmonte,$
            boxwidth=boxwidth[i],doplot=doplot,silent=silent)
          if i eq 0L then lineew = lineew1 else lineew = [ [lineew], [lineew1] ]

       endfor

       return, reform(lineew)
    endif

    if not keyword_set(silent) then splog, 'Fitting '+linename
    
; default settings
    
    if n_elements(nmonte) eq 0L then nmonte = 20L      ; number of Monte-Carlo realizations
    if n_elements(boxwidth) eq 0L then boxwidth = 35.0 ; Angstrom (should be about 20 Angstrom for H-gamma)

    nterms = 5L             ; Gaussian + line
    cd1_1 = wave[1]-wave[0] ; Angstrom/pixel

    lineew = lineew_init(linewave,linename) ; output structure
    
; check for an out-of-range line

    if ((linewave+boxwidth) gt max(wave)) or ((linewave-boxwidth) lt min(wave)) then return, lineew       
       
; zoom in on the line

    get_element, wave, linewave+[-boxwidth,+boxwidth], xx

    lwave = wave[xx[0]:xx[1]]
    nwave = n_elements(lwave)
    lineflux = flux[xx[0]:xx[1]]
    lineerr = ferr[xx[0]:xx[1]]

; fit a negative Gaussian + line to the line profile; Monte Carlo to
; get the error in the flux and the EW (especially from the continuum
; at line center)

    continuum_arr = fltarr(nmonte)
    boxflux_arr = fltarr(nmonte)
    gaussflux_arr = fltarr(nmonte)
    
    lflux = lineflux

; constrain the wavelength of the Gaussian

    parinfo = {value: 0D, limited: [0,0], limits: [0.0D,0.0D]}
    parinfo = replicate(parinfo,nterms)
    parinfo[1].value = linewave
    parinfo[1].limited = [1L,1L]
    parinfo[1].limits = linewave + [-5.0,+5.0] ; +/- 5 Angstrom

    get_element, lwave, parinfo[1].limits, xx
    parinfo[0].value = max(lflux[xx[0]:xx[1]])
    
    for j = 0L, nmonte-1L do begin

       yfit = mpfitpeak(lwave,lflux,fitcoeff,nterms=nterms,error=lineerr,$
         /gaussian,perror=perror,quiet=quiet,bestnorm=bestnorm,parinfo=parinfo)

       contfit = poly(lwave,fitcoeff[3:4])
       continuum_arr[j] = interpol(contfit,lwave,linewave)

       indx = where((lwave ge fitcoeff[1]-3.0*fitcoeff[2]) and $
                    (lwave le fitcoeff[1]+3.0*fitcoeff[2]),nindx)

       if nindx ne 0L then begin

          boxflux_arr[j] = total(lflux[indx]-contfit[indx])*cd1_1   ; box flux
          gaussflux_arr[j] = sqrt(2.0*!dpi)*fitcoeff[2]*fitcoeff[0] ; gaussian flux

          if keyword_set(doplot) and j eq 0L then begin
             
             scale = 1E15
             djs_plot, lwave, scale*lflux, xsty=3, ysty=3, ps=10, charsize=1.5, charthick=2.0, $
               xtitle='Wavelength ('+angstrom()+')', ytitle='f_{\lambda} (10^{15} '+flam_units()+')'
             djs_oplot, lwave, scale*yfit, ps=10, thick=2.0, color='green'
             djs_oplot, lwave, scale*contfit, thick=2.0, color='red'
             djs_oplot, [lwave[indx[0]],lwave[indx[0]]], !y.crange, line=2, thick=2.0
             djs_oplot, [lwave[indx[nindx-1L]],lwave[indx[nindx-1L]]], !y.crange, line=2, thick=2.0
             cc = get_kbrd(1)

          endif
       
       endif

       lflux = lineflux + randomn(seed,nwave)*lineerr

    endfor

; construct the equivalent width and error

    continuum = djs_mean(continuum_arr) & continuum_err = stddev(continuum_arr)
    boxflux = djs_mean(boxflux_arr) & boxflux_err = stddev(boxflux_arr)
    gaussflux = djs_mean(gaussflux_arr) & gaussflux_err = stddev(gaussflux_arr)

    ew = (boxflux/continuum + gaussflux/continuum) / 2.0
    ew_err = stddev([boxflux/continuum,gaussflux/continuum])

;   gauss_ew_err = gaussflux_err^2.0/continuum^2.0 + (gaussflux/continuum^2.0)^2.0*continuum_err^2.0
;   box_ew_err = boxflux_err^2.0/continuum^2.0 + (boxflux/continuum^2.0)^2.0*continuum_err^2.0
;   ew_err = (gauss_ew_err + box_ew_err) / 2.0
    
; fill the output structure    

    lineew.continuum = [continuum,continuum_err]
    lineew.boxflux = [boxflux,boxflux_err]
    lineew.gaussflux = [gaussflux,gaussflux_err]
    lineew.ew = [ew,ew_err]

    if keyword_set(doplot) then cc = get_kbrd(1)

return, lineew
end
