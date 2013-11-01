;+
; NAME:
;       MEASURE_LINEEW()
;
; PURPOSE:
;       Measure emission- and absorption-line equivalent widths. 
;
; CALLING SEQUENCE:
;       linefit = measure_lineew(linefit,wave,flux,ferr,$
;          [snrcut=,speclineewfit=],/absorption,/overwrite,/doplot)
;
; INPUTS:
;       linefit - fitting structure outputed from ILINEFIT() 
;       wave    - spectrum wavelength vector [NPIX]
;       flux    - observed galaxy spectrum [erg/s/cm2/Angstrom] [NPIX]
;       ferr    - corresponding flux error
;
; OPTIONAL INPUTS:
;       snrcut    - compute upper limits on lines with S/N < SNRCUT
;                   (default 3.0)
;       lineres   - see ILINEFIT()
;
; KEYWORD PARAMETERS:
;       absorption - measure absorption lines rather than the default
;                    emission lines
;       overwrite  - overwrite the initial guess fitting parameters 
;                    in LINEFIT 
;       doplot     - generate a plot of the line fitting
;
; OUTPUTS:
;       linefit    - (modified) upper limits and line equivalent
;                    widths have been computed
;
; OPTIONAL OUTPUTS:
;       speclineewfit - spectral fit to every emission/absorption line
;                       [NPIX] (used for plotting)
;
; COMMENTS:
;       Measure equivalent widths.  we use the results of the
;       simultaneous gaussian fitting done in ILINEFIT() as
;       constrained guesses to the fitting here, but now we
;       incorporate a local linear continuum term to measure
;       equivalent widths and errors.  NOTE: for blended lines the
;       LINEEW_BOX equivalent width measurement is not sensible, just
;       as the LINEAREA measurement is not sensible 
;
; PROCEDURES USED:
;       GET_ELEMENT, IUPPER_LIMITS(), ILINEFIT(), POLY_ARRAY(),
;       DJS_PLOT, DJS_OPLOT, LEGEND, DJS_MEDIAN()
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 June 11, U of A, excised from IFITSPEC()
;-

function measure_lineew, linefit, wave, flux, ferr, snrcut=snrcut, $
  lineres=lineres, speclineewfit=speclineewfit, absorption=absorption, $
  overwrite=overwrite, doplot=doplot, _extra=extra

    nline = n_elements(linefit)
    nwave = n_elements(wave)
    nflux = n_elements(flux)
    nferr = n_elements(ferr)

    if (nline eq 0L) or (nwave eq 0L) or (nflux eq 0L) or (nferr eq 0L) then begin
       print, 'Syntax - linefit = measure_lineew(linefit,wave,flux,ferr,$'
       print, '   [snrcut=,speclineewfit=],/absorption,/overwrite,/doplot)'
       return, -1L 
    endif
    
    if n_elements(snrcut) eq 0L then snrcut = 3.0
    if n_elements(speclineewfit) eq 0L then speclineewfit = flux*0.0D
    
    light = 2.99792458D5             ; speed of light [km/s]
    medz = djs_median(linefit.linez)

; the main FOR loop is on each emission line or set of blended lines
     
    blend = strcompress(linefit.line_blend,/remove)
    blendindx = blend[uniq(blend)]

    boxwidth = 13.0 ; units of Gaussian sigma
    
    for i = 0L, n_elements(blendindx)-1L do begin

       ii = where(blend eq blendindx[i],blendcount)
       
; isolate the emission line(s) and restrict the wavelength range 
       
       bline = linefit[ii]
       zline = bline.linez
       meanz = djs_mean(zline)
       meanwave = djs_mean(bline.linewave)*(1+meanz)

       if (blendcount eq 1L) then begin
          leftbox  = meanwave - boxwidth*meanwave*lineres[ii]/light ; [Angstrom]
          rightbox = meanwave + boxwidth*meanwave*lineres[ii]/light
       endif else begin
          leftbox  = bline[0].linewave*(1+meanz) - boxwidth*meanwave*lineres[ii[0]]/light
          rightbox = bline[blendcount-1L].linewave*(1+meanz) + boxwidth*meanwave*lineres[ii[blendcount-1L]]/light
       endelse

       get_element, wave, [leftbox,rightbox], xx

       fwave = wave[xx[0]:xx[1]]
       fflux = flux[xx[0]:xx[1]]
       fferr = ferr[xx[0]:xx[1]]
       fivar = 1.0/fferr^2.0
       nfpix = n_elements(fwave)

;      splog, '['+strjoin(bline.linename,',')+'] '+strjoin(string(minmax(fwave),format='(G0.0)'),'-')

; place upper limits on undetected emission lines (LINEAREA_ERR = -1.0
; or LINEAREA/LINEAREA_ERR < SNRCUT)

       ulimit = where((bline.linearea_err eq -1.0) or $
         (bline.linearea/bline.linearea_err lt snrcut),$
         nulimit,comp=goodfit,ncomp=ngoodfit)

; compute upper limits for single unblended lines, e.g., [OIII], or
; for blended lines, e.g., [OII] 3727 doublet, only if *all* the
; blended lines were not detected

       if (nulimit gt 0L) and (nulimit eq blendcount) then begin

          splog, 'Computing upper limit(s) on ['+strjoin(bline.linename,', ')+'].'
          linefit[ii[ulimit]] = iupper_limits(fflux,linefit[ii[ulimit]],snrcut=snrcut)

          plotnum = 1L
          
       endif else begin 

; fit single emission lines.  also fit blended lines if at least one
; line was detected.  there is an assumption in SIGGUESS that the line
; widths have been tied together, but hopefully MPFIT() will behave
; properly even if SIGGUESS is off
          
          splog, 'Fitting ['+strjoin(bline.linename,', ')+'].'

          if n_elements(lineres) ne 0L then lineres1 = lineres[ii]

          background = poly_array(nfpix,2) ; linear continuum
          sigguess = bline.linesigma/alog(10.0)/light
          zindex = replicate(1,blendcount)
          windex = replicate(1,blendcount) ; tie the widths together

          if keyword_set(absorption) then $
            fvalue = replicate(min(fflux),blendcount) else $
            fvalue = replicate(max(fflux),blendcount)

          result = ilinefit(double(fflux),fwave,bline.linewave,$
            lineres1,invvar=fivar,linename=bline.linename,zindex=zindex,$
            windex=windex,fvalue=fvalue,background=background,$
            zguess=zline,sigguess=sigguess,specfit=yfit,bfit=bfit,$
            allowneg=absorption,_extra=extra)

          speclineewfit[xx[0]:xx[1]] = speclineewfit[xx[0]:xx[1]] + yfit
          
;         niceprint, bline.linename, result.linecontlevel, result.linecontlevel_err

          if keyword_set(overwrite) then linefit[ii] = result else begin
          
             linefit[ii].linecontlevel = result.linecontlevel
             linefit[ii].linecontlevel_err = result.linecontlevel_err

          endelse

          plotnum = 2L

; finally compute upper limits on blended lines that were tied to a
; detected line.  to calculate the upper limits properly we have to
; subtract the Gaussian model from the data

          if (nulimit gt 0L) and (nulimit lt blendcount) then begin

             residuals = fflux - yfit + bfit
             linefit[ii[ulimit]] = iupper_limits(residuals,linefit[ii[ulimit]],snrcut=snrcut)

             plotnum = 3L
             
          endif
             
       endelse 

;      niceprint, bline.linename, linefit[ii].linecontlevel, linefit[ii].linecontlevel_err
       
; compute the equivalent widths and errors

       linefit[ii] = fill_ew(linefit[ii])
       
       if keyword_set(doplot) then begin

          if keyword_set(absorption) then scale = 1.0 else scale = 1E15
          djs_plot, wave, scale*flux, xsty=11, ysty=3, xrange=[leftbox,rightbox]*[0.99,1.01], $
            ps=10, charsize=2.0, charthick=2.0, xthick=2.0, ythick=2.0, ymargin=[4,3], $
            xtitle='Wavelength '+angstrom()+')', ytitle='f_{\lambda} (10^{15} '+flam_units()+')'
          axis, /xaxis, xrange=[leftbox,rightbox]*[0.99,1.01]/(1+meanz), xthick=2.0, $
            charsize=2.0, charthick=2.0, xtitle='Rest Wavelength ('+angstrom()+')', xsty=3

          case plotnum of
             1L: begin
                djs_oplot, !x.crange, scale*linefit[ii[ulimit[0]]].linecontlevel*[+1,+1], line=0, thick=2.0
                djs_oplot, !x.crange, scale*(linefit[ii[ulimit[0]]].linecontlevel+$
                  linefit[ii[ulimit[0]]].linecontlevel_err*[-1,-1]), line=2, thick=2.0
                djs_oplot, !x.crange, scale*(linefit[ii[ulimit[0]]].linecontlevel+$
                  linefit[ii[ulimit[0]]].linecontlevel_err*[+1,+1]), line=2, thick=2.0
                legend, 'Upper limit on ['+strjoin(bline[ulimit].linename,', ')+']', /left, /top, $
                  box=0, charsize=2.0, charthick=2.0
             end
             2L: begin
                djs_oplot, fwave, scale*yfit, color='red', thick=4.0, ps=10
                djs_oplot, fwave, scale*bfit, color='navy', thick=4.0
                legend, bline.linename, /left, /top, box=0, charsize=2.0, charthick=2.0
             end
             3L: begin
                djs_oplot, fwave, scale*yfit, color='red', thick=4.0, ps=10
                djs_oplot, fwave, scale*bfit, color='navy', thick=4.0
                djs_oplot, !x.crange, scale*linefit[ii[ulimit[0]]].linecontlevel*[+1,+1], line=0, thick=2.0
                djs_oplot, !x.crange, scale*(linefit[ii[ulimit[0]]].linecontlevel+$
                  linefit[ii[ulimit[0]]].linecontlevel_err*[-1,-1]), line=2, thick=2.0
                djs_oplot, !x.crange, scale*(linefit[ii[ulimit[0]]].linecontlevel+$
                  linefit[ii[ulimit[0]]].linecontlevel_err*[+1,+1]), line=2, thick=2.0
                legend, ['Fitted '+bline[goodfit].linename,$
                  'Upper limit on ['+strjoin(bline[ulimit].linename,', ')+']'], $
                  /left, /top, box=0, charsize=2.0, charthick=2.0
             end
             else: 
          endcase

          cc = get_kbrd(1)

       endif 

    endfor 

return, linefit
end
