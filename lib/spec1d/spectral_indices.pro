;+
; NAME:
;       SPECTRAL_INDICES()
;
; PURPOSE:
;       Measure spectral indices (D4000, 41-50, and the Lick indices)
;       in a galaxy spectrum. 
;
; CALLING SEQUENCE:
;       indices = spectral_indices(wave,flux,ferr,$
;          [indexpath=,indexfile=],/doplot)
;
; INPUTS:
;       wave     - wavelength vector [NPIX]
;       flux     - spectrum flux in erg/s/cm2/Angstrom [NPIX]
;       ferr     - corresponding error spectrum [NPIX]
;
; OPTIONAL INPUTS:
;       indexpath - path name to INDEXFILE
;       indexfile - name of the file listing the indices
;
; KEYWORD PARAMETERS:
;       doplot    - generate a plot of the index measurements 
;
; OUTPUTS:
;       indices   - data structure with all the results
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       The indices are uncorrected for emission.  We measure D4000
;       (Bruzual 1983, ApJ, 273, 105), D4000_narrow (Balogh et
;       al. 1999, ApJ, 527, 54), 41-50 (Kennicutt 1992, ApJS, 79,
;       255), and all the Lick indices (see the default INDEXFILE for
;       references)
;
; DATA FILES:
;	${ISPEC_DIR}/etc/indexlist.dat
;
; PROCEDURES USED:
;       DJS_PLOT, DJS_OPLOT, GET_ELEMENT, TSUM(), DJS_ITERSTAT,
;       READCOL, IM_SYMBOLS()
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 October 31, U of A
;-

function spectral_indices, wave, flux, ferr, indexpath=indexpath, $
  indexfile=indexfile, doplot=doplot

    npix = n_elements(wave)
    nflux = n_elements(flux)
    nferr = n_elements(ferr)
    
    if (npix eq 0L) or (nflux eq 0L) or (nferr eq 0L) then begin
       print, 'Syntax - indices = spectral_indices(wave,flux,ferr,$'
       print, '   [indexpath=,indexfile=],/doplot)'
       return, -1L
    endif

    if (npix ne nflux) or (npix ne nferr) then begin
       print, 'Dimensions of WAVE, FLUX, and FERR do not agree.'
       return, -1L
    endif
    
    if n_elements(indexpath) eq 0L then $
      indexpath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='etc')
    if n_elements(indexfile) eq 0L then indexfile = 'indexlist.dat'

    im_symbols, 'circle', /fill, psize=2

; convert to F_nu (erg/s/cm2/Hz)

    light = 2.99792458D18
    fnu = wave*wave*flux/light
    fnu_err = wave*wave*ferr/light
    nu = light/wave
    
    colornames = ['D4000','D4000_NARROW','C41_50'] ; required for the output structure
    colors = create_struct('D4000', [0.0,-1.0], 'D4000_narrow', [0.0,-1.0], 'c41_50', [0.0,-1.0])

; measure D(4000)

; Bruzual 1983, ApJ, 273, 105
; ---------------------------

    wr = [4050.0,4250.0]        ; red wavelength interval
    wb = [3750.0,3950.0]        ; blue wavelength interval
       
    if (min(wave) lt wb[0]) and (max(wave) gt wr[1]) then begin
    
       get_element, wave, wr, rr
       get_element, wave, wb, bb

       red = tsum(wave,fnu,rr[0],rr[1])
       red_err = sqrt(tsum(wave,fnu_err^2.0,rr[0],rr[1]))

       blue = tsum(wave,fnu,bb[0],bb[1])
       blue_err = sqrt(tsum(wave,fnu_err^2.0,bb[0],bb[1]))

       D4000 = (wb[1]-wb[0]) / (wr[1]-wr[0]) * red / blue
       D4000_err = (wb[1]-wb[0]) / (wr[1]-wr[0]) * $
         im_compute_error(red,red_err,blue,blue_err,/quotient)

       colors.D4000 = [D4000,D4000_err]
       
       if keyword_set(doplot) then begin

          djs_plot, wave[(bb[0]-20)>0L:(rr[1]+20)<(npix-1)], fnu[(bb[0]-20)>0L:(rr[1]+20)<(npix-1)], $
            color='yellow', xsty=3, ysty=3, ps=10, xtitle='Rest Wavelength ('+angstrom()+')', $
            ytitle='f_{\nu}', charsize=2, charthick=2.0, title='D4000 (Bruzual)', xmargin=[13,3]

          polyfill, [wave[bb[0]],wave[bb[1]],wave[bb[1]],wave[bb[0]]], $
            [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
            /line_fill, orientation=45, color=djs_icolor('blue')
          
          polyfill, [wave[rr[0]],wave[rr[1]],wave[rr[1]],wave[rr[0]]], $
            [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
            /line_fill, orientation=45, color=djs_icolor('red')

          cc = get_kbrd(1)
          
       endif
       
    endif
       
; Balogh et al. 1999, ApJ, 527, 54
; ---------------------------

    wb = [3850.0,3950.0]        ; blue wavelength interval
    wr = [4000.0,4100.0]        ; red wavelength interval
       
    if (min(wave) lt wb[0]) and (max(wave) gt wr[1]) then begin

       get_element, wave, wr, rr
       get_element, wave, wb, bb

       red = tsum(wave,fnu,rr[0],rr[1])
       red_err = sqrt(tsum(wave,fnu_err^2.0,rr[0],rr[1]))

       blue = tsum(wave,fnu,bb[0],bb[1])
       blue_err = sqrt(tsum(wave,fnu_err^2.0,bb[0],bb[1]))

;      red = integral(wave,fnu,wr[0],wr[1])
;      red_err = sqrt(integral(wave,fnu_err^2.0,wr[0],wr[1]))
;      blue = integral(wave,fnu,wb[0],wb[1])
;      blue_err = sqrt(integral(wave,fnu_err^2.0,wb[0],wb[1]))

       D4000_narrow = (wb[1]-wb[0]) / (wr[1]-wr[0]) * red / blue
       D4000_narrow_err = (wb[1]-wb[0]) / (wr[1]-wr[0]) * $
         im_compute_error(red,red_err,blue,blue_err,/quotient)

       colors.D4000_narrow = [D4000_narrow,D4000_narrow_err]
       
       if keyword_set(doplot) then begin

          djs_plot, wave[(bb[0]-20)>0L:(rr[1]+20)<(npix-1)], fnu[(bb[0]-20)>0L:(rr[1]+20)<(npix-1)], $
            color='yellow', xsty=3, ysty=3, ps=10, xtitle='Rest Wavelength ('+angstrom()+')', $
            ytitle='f_{\nu}', charsize=2, charthick=2.0, title='D4000 Narrow (Balogh)', xmargin=[13,3]
          
          polyfill, [wave[bb[0]],wave[bb[1]],wave[bb[1]],wave[bb[0]]], $
            [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
            /line_fill, orientation=45, color=djs_icolor('blue')
          
          polyfill, [wave[rr[0]],wave[rr[1]],wave[rr[1]],wave[rr[0]]], $
            [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
            /line_fill, orientation=45, color=djs_icolor('red')
          
          cc = get_kbrd(1)
          
       endif
       
    endif 
       
; measure the 41-50 color index (Kennicutt 1992, ApJS, 79, 255)

    wb = [4050.0,4250.0]        ; blue wavelength interval
    wr = [4900.0,5100.0]        ; red wavelength interval

    if (min(wave) lt wb[0]) and (max(wave) gt wr[1]) then begin

       get_element, wave, wr, rr
       get_element, wave, wb, bb

       djs_iterstat, fnu[rr[0]:rr[1]], sigrej=0.8, maxiter=25, mean=red, sigma=red_err
;      red = djs_mean(fnu[rr[0]:rr[1]])
;      red_err = stddev(fnu[rr[0]:rr[1]])/(rr[1]-rr[0]+1)

       djs_iterstat, fnu[bb[0]:bb[1]], sigrej=0.8, maxiter=25, mean=blue, sigma=blue_err

       c41_50 = 2.5*alog10(red/blue)
       c41_50_err = (2.5/alog(10))*im_compute_error(red,red_err,blue,blue_err,/quotient)/(red/blue)

       colors.c41_50 = [c41_50,c41_50_err]

       if keyword_set(doplot) then begin

          djs_plot, wave[(bb[0]-20)>0L:(rr[1]+20)<(npix-1)], fnu[(bb[0]-20)>0L:(rr[1]+20)<(npix-1)], $
            color='yellow', xsty=3, ysty=3, ps=10, xtitle='Rest Wavelength ('+angstrom()+')', $
            ytitle='f_{\nu}', charsize=2, charthick=2.0, title='41-50 Kennicutt Index', xmargin=[13,3]
          
          polyfill, [wave[bb[0]],wave[bb[1]],wave[bb[1]],wave[bb[0]]], $
            [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
            /line_fill, orientation=45, color=djs_icolor('blue')
          
          polyfill, [wave[rr[0]],wave[rr[1]],wave[rr[1]],wave[rr[0]]], $
            [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
            /line_fill, orientation=45, color=djs_icolor('red')

          bwave = djs_mean(wb)
          rwave = djs_mean(wr)
          
          plots, bwave, blue, ps=8, color=djs_icolor('white')
          plots, rwave, red, ps=8, color=djs_icolor('white')

          cc = get_kbrd(1)
          
       endif

    endif

; measure the Lick indices

    if file_test(indexpath+indexfile,/regular) eq 0L then begin

       splog, 'Lick index file '+indexpath+indexfile+' not found.'
       indices = colors

     endif else begin
       
       readcol, indexpath+indexfile, licknames, w1, w2, wb1, wb2, wr1, wr2, units, $
         format='A,F,F,F,F,F,F,A', /silent, comment='#'
       nlick = n_elements(licknames)

; initialize the output structure

       indices = create_struct(licknames[0],[0.0,-1.0])
       for j = 1L, nlick-1L do indices = create_struct(indices,licknames[j],[0.0,-1.0])
       
;      indices = create_struct('indices', licknames)
;      for j = 0L, nlick-1L do indices = create_struct(indices,licknames[j],[0.0,-1.0])
;      indices = create_struct(licknames[0],0.0,licknames[0]+'_err',-1.0)
;      for j = 1L, nlick-1L do indices = create_struct(indices,licknames[j],0.0,licknames[j]+'_err',-1.0)
       
       for k = 0L, nlick-1L do begin

          if (wb1[k] gt min(wave)) and (wr2[k] lt max(wave)) then begin
             
             get_element, wave, [w1[k],w2[k]], linexx ; line wavelengths
             localflux = flux[linexx[0]:linexx[1]]
             localferr = ferr[linexx[0]:linexx[1]]
             localwave = wave[linexx[0]:linexx[1]]
             nlocal = n_elements(localwave)
             
             get_element, wave, [wb1[k],wb2[k]], bluexx ; blue continuum
             blue = integral(wave,flux,wb1[k],wb2[k]) / (wb2[k] - wb1[k])
             blue_err = sqrt(integral(wave,ferr^2.0,wb1[k],wb2[k]) / (wb2[k] - wb1[k]))
             bwave = djs_mean([wb1[k],wb2[k]])
             
             get_element, wave, [wr1[k],wr2[k]], redxx ; red continuum
             red = integral(wave,flux,wr1[k],wr2[k]) / (wr2[k] - wr1[k])
             red_err = sqrt(integral(wave,ferr^2.0,wr1[k],wr2[k]) / (wr2[k] - wr1[k]))
             rwave = djs_mean([wr1[k],wr2[k]])

             plotwave = wave[bluexx[0]:redxx[1]]
             
; compute the continuum
             
             if strmatch(licknames[k],'*BH*',/fold) eq 1B then begin

                continuum = ((red+blue)/2.0)[0]
                continuum_err = (sqrt(red_err^2.0+blue_err^2.0)/2.0)[0]

                plotcont = replicate(continuum,n_elements(plotwave))
                
             endif else begin

; SVDFIT gives the same answer but with errors
                
                slope = (red-blue)/(rwave-bwave) 
                coeff = [red-slope*rwave,slope]
                sigma = coeff*0.05 ; 5% error
                
;               coeff = svdfit([bwave,rwave],[blue,red],2,weights=1.0/[blue_err,red_err],sigma=sigma)
                
                continuum = poly(localwave,coeff)
                continuum_err = sqrt(sigma[0]^2.0 + sigma[1]^2.0*(localwave)^2.0)

                plotcont = poly(plotwave,coeff)
                
             endelse
             
             case units[k] of

                'A': begin

                   index = integral(localwave,1.0-localflux/continuum,localwave[0],localwave[nlocal-1L])
                   err = (localflux/continuum) * sqrt((localferr/localflux)^2.0+(continuum_err/continuum)^2.0)
                   index_err = sqrt(integral(localwave,err^2.0,localwave[0],localwave[nlocal-1L]))

                end
                
                'mag': begin

                   nflux = integral(localwave,localflux/continuum,localwave[0],localwave[nlocal-1L])
                   nerri = (localflux/continuum) * sqrt((localferr/localflux)^2.0 + (continuum_err/continuum)^2.0)
                   nerr = sqrt(integral(localwave,nerri^2.0,localwave[0],localwave[nlocal-1L]))
                   
                   index = -2.5 * alog10(nflux/(w2[k]-w1[k]))
                   index_err = 2.5/alog(10) * nerr / nflux
                   
                end
                
             endcase

             indices.(k) = [index[0],index_err[0]]
;            indices.(k+1) = [index[0],index_err[0]] ; (k+1) offsets from INDICES.INDICES
;            indices.(2*k) = index[0]
;            indices.(2*k+1L) = index_err[0]

             if keyword_set(doplot) then begin
                
                djs_plot, wave[(bluexx[0]-20)>0L:(redxx[1]+20)<(npix-1)], $
                  flux[(bluexx[0]-20)>0L:(redxx[1]+20)<(npix-1)], xmargin=[13,3], $
                  color='yellow', xsty=3, ysty=3, ps=10, xtitle='Rest Wavelength ('+angstrom()+')', $
                  ytitle='f_{\lambda}', charsize=2, charthick=2.0, title=repstr(licknames[k],'_',' ')
                
                polyfill, [wave[bluexx[0]],wave[bluexx[1]],wave[bluexx[1]],wave[bluexx[0]]], $ ; blue continuum
                  [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
                  /line_fill, orientation=45, color=djs_icolor('blue')
                
                polyfill, [wave[redxx[0]],wave[redxx[1]],wave[redxx[1]],wave[redxx[0]]], $ ; red continuum
                  [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
                  /line_fill, orientation=45, color=djs_icolor('red')
                
                plots, bwave, blue, ps=8, color=djs_icolor('white')
                plots, rwave, red, ps=8, color=djs_icolor('white')
;               djs_oplot, wave[bluexx[0]:redxx[1]], poly(wave[bluexx[0]:redxx[1]],coeff), line=0, thick=2, color='cyan'
                djs_oplot, plotwave, plotcont, line=0, thick=2, color='cyan'
                
                polyfill, [wave[linexx[0]],wave[linexx[1]],wave[linexx[1]],wave[linexx[0]]], $ ; line region
                  [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
                  /line_fill, orientation=45

                cc = get_kbrd(1)
                
             endif

          endif
             
       endfor 
    
; append the two structures and return

       indices = create_struct('indices',[colornames,licknames],colors,indices)
;      indices = create_struct(colors,indices)
       
    endelse 

return, indices
end
