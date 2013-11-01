;+
; NAME:
;       ATLAS_FIGURE_FITSPEC
;
; PURPOSE:
;       Generate a plot of the spectral fitting. 
;
; CALLING SEQUENCE:
;       atlas_figure_fitspec, specdata, specfit, $
;          [wave=,flux=,toc=,/postscript
;
; INPUTS:
;       specdata - see ATLAS_FITSPEC()
;       specfit  - see ATLAS_FITSPEC()
;
; OPTIONAL INPUTS:
;       wave     - galaxy wavelength vector [NPIX]
;       flux     - galaxy spectrum [NPIX]
;       toc      - "table of contents" information on SPECFIT
;
; KEYWORD PARAMETERS:
;       postscript - if generating postscript then set this keyword to
;                    prevent initializing IDL windows
;
; OUTPUTS:
;       Three plots are generated.
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;      SPLOG, EMISSION_MASK(), MRDFITS(), PAGEMAKER, DJS_ITERSTAT,
;      DJS_PLOT, DJS_OPLOT, GET_ELEMENT, LEGEND, STRUCT_TRIMTAGS(),
;      PLOTSYM 
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 May 12, U of A, based on an earlier code
;-

pro atlas_figure_fitspec, specdata, specfit, wave=wave, flux=flux, toc=toc, $
  postscript=postscript

    if (n_elements(specdata) eq 0L) or (n_elements(specfit) eq 0L) then begin
       print, 'Syntax - atlas_figure_fitspec, specdata, specfit, $'
       print, '   [wave=,flux=,toc=,/postscript'
       return
    endif

    nwave = n_elements(wave)
    npix = n_elements(flux)

    if (nwave eq 0L) or (npix eq 0L) then begin
       splog, 'WAVE and FLUX must be passed.'
       return
    endif
    if nwave ne npix then begin
       splog, 'WAVE and FLUX must have the same number of elements.'
       return
    endif
    
; define some convenient internal variables

    continuum = specfit[*,0]
    speclinefit = specfit[*,2]
;   speclineewfit = specfit[*,3]
    espectrum = flux-continuum
    charsize = 1.1
    galaxy = strcompress(specdata.galaxy,/remove)
    specfile = repstr(repstr(strcompress(specdata.specfile,/remove),'.fits',''),'.ms','')

    mask = emission_mask(wave,z=specdata.z_abs,width=30.0,bad=maskpix,good=keeppix,/telluric)

;    starflux = mrdfits(toc.eigendir+toc.eigenfile,0,eigenhead,/silent)
;    starwave = make_wave(eigenhead)
;    starcoeff = specdata.starcoeff
;    keep = where(starcoeff gt 0.01,nstar)
;    starflux = starflux[*,keep]
;    starcoeff = starcoeff[keep]/total(starcoeff[keep])

    lineres = specdata.specres
    
    nline = n_elements(specdata.linename)
    nabsline = n_elements(specdata.abs_linename)

    plotsym, 0, 1, /fill
    
; ---------------------------------------------------------------------------
; 2-panel plot    
; ---------------------------------------------------------------------------

; initialize the page 1 postscript
    
    pagemaker, nx=1, ny=2, position=position, /normal, $
      xmargin=[0.6,0.2], ymargin=[0.5,1.0], yspace=0.5

    if not keyword_set(postscript) then window, 0, xs=450, ys=450    

; --------------------------------------------------    
; PANEL 1
; --------------------------------------------------    
    
    residuals = (flux-continuum)/flux
    
    x = wave[0:npix/2] & y = flux[0:npix/2]

;   djs_iterstat, y, median=md, sigma=sig, mean=norm, mask=mmask
    djs_iterstat, y, median=md, sigma=sig, mean=mn, mask=mmask
    norm = max(y[where(mmask)])
    smin = min(y[where(mmask)])/norm

    yrange = fltarr(2)
    yrange[0] = ((md-3*sig)/norm)>smin
    yrange[1] = ((md+6*sig)/norm)<1.1
    
    xrange = minmax(x)

;   norm = max(y-speclinefit[0:npix/2])
;   yrange = minmax(abs(y-speclinefit[0:npix/2])/norm)*[0.6,1.05]
        
    djs_plot, x, y/norm, ps=10, ysty=3, xsty=3, yrange=yrange, xthick=2.0, ythick=2.0, $
      ytickname=replicate(' ',10), thick=2.0, xrange=xrange, $
      ytitle=ytitle, charsize=charsize, charthick=2.0, position=position[*,0]
    djs_oplot, x, continuum[0:npix/2]/norm, ps=10, color='red', thick=5.0

    for i = 0L, n_elements(maskpix)-1L do $
      if (wave[maskpix[i]] gt !x.crange[0]) and (wave[maskpix[i]] lt !x.crange[1]) and $
      (flux[maskpix[i]]/norm gt !y.crange[0]) and (flux[maskpix[i]]/norm lt !y.crange[1]) then $
      plots, wave[maskpix[i]], flux[maskpix[i]]/norm, ps=8, color=djs_icolor('purple')

; information legend
    
    label = [toc.templates,'z = '+string(specdata.z_obj,format='(G0.0)')]
    legend, label, /left, /top, charsize=1.5, charthick=2.0, box=0
        
;; residual spectrum    
;
;    residposition = position[*,0]
;    residposition[3] = (residposition[3]-residposition[1])/5.0+residposition[1]
;
;    yresid = residuals[0:npix/2L]
;    yresidrange = max(abs(residuals[keeppix]))*[-0.5,0.5]
;    meanscatter = djs_mean(yresid)
;    
;    djs_plot, [0], [0], /nodata, /noerase, position=residposition, $
;      xthick=5.0, ythick=5.0, xtickname=replicate(' ',10), $ ; ytickname=replicate(' ',10), $
;      xsty=7, ysty=7, yrange=yresidrange, xrange=xrange
;    djs_oplot, x, yresid, ps=10, color='cyan', thick=2.0
;    djs_oplot, !x.crange, [0,0], line=0, thick=2.0
;;   legend, string(100*meanscatter,format='(F4.1)')+'%', /right, /bottom, /box, $
;;     charsize=1.5, charthick=2.0
;    
;    for i = 0L, n_elements(maskpix)-1L do $
;      if (wave[maskpix[i]] gt !x.crange[0]) and (wave[maskpix[i]] lt !x.crange[1]) and $
;      (residuals[maskpix[i]] gt !y.crange[0]) and (residuals[maskpix[i]] lt !y.crange[1]) then $
;      plots, wave[maskpix[i]], residuals[maskpix[i]], ps=8, color=djs_icolor('purple')
    
; --------------------------------------------------    
; PANEL 2    
; --------------------------------------------------    
        
    x = wave[npix/2+1:npix-1] & y = flux[npix/2+1:npix-1]

;   djs_iterstat, y, median=md, sigma=sig, mean=norm
    djs_iterstat, y, median=md, sigma=sig, mean=mn, mask=mmask
    norm = max(y[where(mmask)])
    smin = min(y[where(mmask)])/norm

    yrange = fltarr(2)
    yrange[0] = ((md-3*sig)/norm)>smin
    yrange[1] = ((md+6*sig)/norm)<1.1
;   yrange = (md+sig*[-3,6])/norm

    xrange = minmax(x)
    
;   norm = max(y-speclinefit[npix/2+1:npix-1])
;   yrange = minmax(abs(y-speclinefit[npix/2+1:npix-1])/norm)*[0.6,1.05]
    
    djs_plot, x, y/norm, ps=10, ysty=3, xsty=3, yrange=yrange, xthick=2.0, ythick=2.0, $
      ytickname=replicate(' ',10), thick=2.0, ytitle=ytitle, $
      charsize=charsize, charthick=2.0, position=position[*,1], /noerase, $
      xtitle='Observed Wavelength ('+angstrom()+')', xrange=xrange
    djs_oplot, x, continuum[npix/2+1:npix-1]/norm, ps=10, color='red', thick=5.0

    for i = 0L, n_elements(maskpix)-1L do $
      if (wave[maskpix[i]] gt !x.crange[0]) and (wave[maskpix[i]] lt !x.crange[1]) and $
      (flux[maskpix[i]]/norm gt !y.crange[0]) and (flux[maskpix[i]]/norm lt !y.crange[1]) then $
      plots, wave[maskpix[i]], flux[maskpix[i]]/norm, ps=8, color=djs_icolor('purple')
       
;; residual spectrum    
;
;    residposition = position[*,1]
;    residposition[3] = (residposition[3]-residposition[1])/5.0+residposition[1]
;
;    yresid = residuals[npix/2L+1:npix-1]
;    yresidrange = max(abs(residuals[keeppix]))*[-0.5,0.5]
;    
;    djs_plot, [0], [0], /nodata, /noerase, position=residposition, $
;      xthick=5.0, ythick=5.0, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
;      xsty=7, ysty=7, yrange=yresidrange, xrange=xrange
;    djs_oplot, x, yresid, ps=10, color='cyan', thick=2.0
;    djs_oplot, !x.crange, [0,0], line=0, thick=2.0
;
;    for i = 0L, n_elements(maskpix)-1L do $
;      if (wave[maskpix[i]] gt !x.crange[0]) and (wave[maskpix[i]] lt !x.crange[1]) and $
;      (residuals[maskpix[i]] gt !y.crange[0]) and (residuals[maskpix[i]] lt !y.crange[1]) then $
;      plots, wave[maskpix[i]], residuals[maskpix[i]], ps=8, color=djs_icolor('purple')

; chi2 of the continuum fit

    legend, textoidl('\chi^2 = '+string(specdata.continuum_chi2,format='(G0.0)')), $
      /left, /top, charsize=1.5, charthick=2.0, box=0
        
; title    
    
    refpos = reform(position[*,0])
    xpos = (refpos[2]-refpos[0])/2.0+refpos[0]
    ypos = refpos[3]*1.01
    
    xyouts, xpos, ypos, galaxy+' ['+specfile+'] Continuum', /normal, charsize=1.5, charthick=2.0, align=0.5

; ---------------------------------------------------------------------------
; 3-panel plot    
; ---------------------------------------------------------------------------

;;;; initialize the page 1 postscript
;;;    
;;;    pagemaker, nx=1, ny=3, position=position, /normal, $
;;;      xmargin=[0.2,0.2], ymargin=[0.5,1.0], yspace=0.5
;;;
;;;    if not keyword_set(postscript) then window, 0, xs=450, ys=450    
;;;
;;;; --------------------------------------------------    
;;;; PANEL 1
;;;; --------------------------------------------------    
;;;    
;;;    x = wave[0:npix/2] & y = flux[0:npix/2]
;;;    norm = max(y-speclinefit[0:npix/2])
;;;    yrange = minmax(abs(y-speclinefit[0:npix/2])/norm)*[0.95,1.05]
;;;    
;;;    djs_plot, x, y/norm, ps=10, ysty=3, xsty=3, yrange=yrange, xthick=2.0, ythick=2.0, $
;;;      ytickname=replicate(' ',10), thick=2.0, $
;;;      ytitle=ytitle, charsize=charsize, charthick=2.0, position=position[*,0]
;;;    djs_oplot, x, continuum[0:npix/2]/norm, ps=10, color='green', thick=2.0
;;;
;;;    for i = 0L, n_elements(maskpix)-1L do $
;;;      if (wave[maskpix[i]] gt !x.crange[0]) and (wave[maskpix[i]] lt !x.crange[1]) and $
;;;      (flux[maskpix[i]]/norm gt !y.crange[0]) and (flux[maskpix[i]]/norm lt !y.crange[1]) then $
;;;      plots, wave[maskpix[i]], flux[maskpix[i]]/norm, ps=8, color=djs_icolor('purple')
;;;       
;;;; --------------------------------------------------    
;;;; PANEL 2    
;;;; --------------------------------------------------    
;;;        
;;;    x = wave[npix/2+1:npix-1] & y = flux[npix/2+1:npix-1]
;;;    norm = max(y-speclinefit[npix/2+1:npix-1])
;;;    yrange = minmax(abs(y-speclinefit[npix/2+1:npix-1])/norm)*[0.95,1.05]
;;;    
;;;    djs_plot, x, y/norm, ps=10, ysty=3, xsty=3, yrange=yrange, xthick=2.0, ythick=2.0, $
;;;      ytickname=replicate(' ',10), thick=2.0, ytitle=ytitle, $
;;;      charsize=charsize, charthick=2.0, position=position[*,1], /noerase
;;;    djs_oplot, x, continuum[npix/2+1:npix-1]/norm, ps=10, color='green', thick=2.0
;;;
;;;    for i = 0L, n_elements(maskpix)-1L do $
;;;      if (wave[maskpix[i]] gt !x.crange[0]) and (wave[maskpix[i]] lt !x.crange[1]) and $
;;;      (flux[maskpix[i]]/norm gt !y.crange[0]) and (flux[maskpix[i]]/norm lt !y.crange[1]) then $
;;;      plots, wave[maskpix[i]], flux[maskpix[i]]/norm, ps=8, color=djs_icolor('purple')
;;;       
;;;; --------------------------------------------------    
;;;; PANEL 3
;;;; --------------------------------------------------    
;;;
;;;    djs_plot, starwave, starcoeff[0]*im_normalize(starflux[*,0],/max), ps=10, $
;;;      ysty=3, xsty=3, xthick=2.0, ythick=2.0, thick=1.0, $
;;;      ytitle=ytitle, charsize=charsize, charthick=2.0, position=position[*,2], $
;;;      color='cyan', yrange=[0.01,1], ytickname=replicate(' ',10), /ylog, /noerase, $
;;;      xtitle='Observed Wavelength ('+angstrom()+')'
;;;    for j = 1L, nstar-1L do djs_oplot, starwave, starcoeff[j]*im_normalize(starflux[*,j],/max), $
;;;      ps=10, thick=2.0, color='cyan'
;;;    
;;;    xyouts, 0.45, 0.97, galaxy, /normal, charsize=1.5, charthick=2.0, align=0.0

; initialize the page 2 postscript

    pagemaker, nx=4, ny=ceil(nline/4.0), position=position, /normal, $
      xmargin=[0.6,0.2], ymargin=[0.5,1.0], yspace=0.0, xspace=0.0

    if not keyword_set(postscript) then window, 2, xs=450, ys=450

    box = 100.0
    tags = tag_names(specdata[0])
    for k = 0L, nline-1L do begin

       if k eq 0L then noerase = 0L else noerase = 1L

       line = repstr(specdata.linename[k],'_',' ')

; emission-line wavelength
       
       match = where(specdata.linename[k]+'_WAVE' eq tags)
       linewave = specdata.(match)*(1.0+specdata.z_line)

       leftbox  = linewave-10.0*lineres/2.35
       rightbox = linewave+10.0*lineres/2.35
       get_element, wave, [leftbox,rightbox], xx

       yrange = minmax(espectrum[xx[0]:xx[1]])*[1.0,2]

       match = where(specdata.linename[k] eq tags)
       errflag = specdata.(match)[1]
       
       if errflag eq -2.0 then begin

          upper = 'Not Measured'
          label = [line,upper]

          djs_plot, [0], [0], /nodata, xsty=3, ysty=3, $
            xthick=2.0, ythick=2.0, thick=3.0, noerase=noerase, ps=10, $
            xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
            position=position[*,k], yrange=yrange

       endif else begin
          
          djs_plot, wave[xx[0]:xx[1]], espectrum[xx[0]:xx[1]], xsty=3, ysty=3, $
            xthick=2.0, ythick=2.0, thick=5.0, noerase=noerase, ps=10, $
            xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
            position=position[*,k], yrange=yrange
          djs_oplot, wave[xx[0]:xx[1]], speclinefit[xx[0]:xx[1]], thick=5.0, $
            color='red', ps=10
          djs_oplot, linewave*[1,1], !y.crange, line=2, thick=5.0, color='navy'

          case errflag of
             -3.0: label = [line,'Upper Limit']
             else: begin
                snr = 'S/N = '+string((specdata.(match))[0]/(specdata.(match))[1],format='(I0)')
                label = [line,snr]
             endelse
          endcase
       
       endelse

       legend, label, /left, /top, box=0, charsize=1.0, charthick=2.0       

    endfor

; title    
    
    xyouts, xpos, ypos, galaxy+' Line Fluxes', /normal, charsize=1.5, charthick=2.0, align=0.5
;   xyouts, xpos, ypos, galaxy+' ['+specfile+'] Line Fluxes', /normal, charsize=1.5, charthick=2.0, align=0.5
    
; initialize the page 3 postscript

;   nall = nline + nabsline
    
    pagemaker, nx=4, ny=ceil(nline/4.0), position=position, /normal, $
      xmargin=[0.6,0.1], ymargin=[0.5,0.2], yspace=0.0, xspace=0.0

    if not keyword_set(postscript) then window, 3, xs=450, ys=450

    minibox = 50.0
    tags = tag_names(specdata[0])

; plot the emission lines

    for k = 0L, nline-1L do begin

       if k eq 0L then noerase = 0L else noerase = 1L

       line = repstr(specdata.linename[k],'_',' ')
       
; emission-line wavelength
       
       match = where(specdata.linename[k]+'_WAVE' eq tags)
       linewave = specdata.(match)*(1.0+specdata.z_line)

       leftbox  = linewave-15.0*lineres/2.35
       rightbox = linewave+15.0*lineres/2.35
       get_element, wave, [leftbox,rightbox], xx

; error flag       

       match = where(specdata.linename[k]+'_EW' eq tags)
       errflag = specdata.(match)[1]

       if (errflag eq -2.0) then begin

          upper = 'Not Measured'
          label = [line,upper]

          djs_plot, [0], [0], /nodata, xsty=3, ysty=3, $
            xthick=2.0, ythick=2.0, thick=3.0, noerase=noerase, ps=10, $
            xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
            position=position[*,k], yrange=yrange

       endif else begin

; continuum measurement
       
          cmatch = where(specdata.linename[k]+'_CONTINUUM' eq tags)
          line_continuum = (specdata.(cmatch))[0]

          djs_iterstat, flux[xx[0]:xx[1]], sigma=csig, sigrej=2.5
;         djs_iterstat, continuum[ww[0]:ww[1]], sigma=csig, sigrej=5.0
          yrange = line_continuum+[-5,+4]*csig

          djs_plot, wave[xx[0]:xx[1]], flux[xx[0]:xx[1]], xsty=3, ysty=3, $
            xthick=2.0, ythick=2.0, thick=5.0, noerase=noerase, ps=10, $
            xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
            position=position[*,k], yrange=yrange

          djs_oplot, wave[xx[0]:xx[1]], continuum[xx[0]:xx[1]]+speclinefit[xx[0]:xx[1]], $
            color='red', ps=10, thick=5.0
;         djs_oplot, wave[xx[0]:xx[1]], continuum[xx[0]:xx[1]], color='green', ps=10, thick=2.0
          djs_oplot, [wave[xx[0]],wave[xx[1]]], line_continuum*[1,1], $
            line=0, thick=5.0, color='navy'
          djs_oplot, linewave*[1,1], !y.crange, line=2, thick=5.0, color='navy'

          case errflag of
             -3.0: label = [line,'Upper Limit']
             else: begin
                snr = 'S/N = '+string((specdata.(match))[0]/(specdata.(match))[1],format='(I0)')
                label = [line,snr]
             endelse
          endcase
       
       endelse
          
       legend, label, /left, /bottom, box=0, charsize=1.0, charthick=2.0       

    endfor

; title    
    
    xyouts, xpos, ypos, galaxy+' EWs', /normal, charsize=1.5, charthick=2.0, align=0.5

return
end    
