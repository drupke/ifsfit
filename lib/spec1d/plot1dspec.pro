;+
; NAME:
;	PLOT1DSPEC
;
; PURPOSE:
;	Plot one-dimensional extracted spectra.
;
; CALLING SEQUENCE:
;
; INPUTS:
;	specname - name of the multi-extension spectrum
;
; OPTIONAL INPUTS:
;	datapath - path to the spectrum (default PWD)
;	outpath  - 
;	nsmooth  - boxcar smooth with NSMOOTH width
;	psname   - name of the postscript plot (default to OBJECT) 
;	extra    - keywords for SPLOT
;
; KEYWORD PARAMETERS:
;	postscript - generate postscript output
;	all        - plot all FITS files in the current directory, or
;                    in DATAPATH
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;
;
; EXAMPLE:
;
;
; PROCEDURES USED:
;	SPLOT, SOPLOT, SMOOTH, MRDFITS(), TEXTOIDL, PLOTFAVES, CWD() 
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 October 19, U of A
;-

pro plotinfo, specname, header, xtitle=xtitle, ytitle=ytitle, $
  legendbox=legendbox, psname=psname, _extra=extra, all=all

    galaxy = sxpar(header,'GALAXY',count=count)
    if count eq 0L then title = '' else title = strupcase(galaxy)
    object = ''
;   object = ' ('+strupcase(sxpar(header,'OBJECT'))+' )'
    extract = sxpar(header,'EXTRACT')
    scanlen = sxpar(header,'SCANLEN')
    if scanlen eq float(0) then scanlen = 2.5
    aperture = strn(extract,format='(G0.0)')+'" x '+strn(scanlen,format='(G0.0)')+'" Aperture'

    legendbox = [title+object,aperture]

;   fluxcor = sxpar(header,'FLUXCOR')
;   if fluxcor eq 0L then ytitle = 'Counts' else $
;     ytitle = textoidl('f_{\lambda} (erg s^{-1} cm^{-2} '+angstrom()+'^{-1})')

    psname = strmid(specname,0,strpos(specname,'.ms.fits'))

return
end

pro plot1dspec, specname, datapath=datapath, nsmooth=nsmooth, postscript=postscript, $
  outpath=outpath, psname=psname, header=header, overplot=overplot, _extra=extra, $
  left=left, right=right, all=all, normalize=normalize

    nspec = n_elements(specname)
    if (nspec eq 0L) and (n_elements(all) eq 0L) then begin
       print, 'Syntax - plot1dspec, specname, [datapath=], [nsmooth=], _extra=extra'
       return
    endif

    if not keyword_set(datapath) then datapath = cwd()
    if not keyword_set(outpath) then outpath = datapath

; read all the FITS files in the current directory
    
    if n_elements(all) ne 0L then begin
       pushd, datapath
       specname = findfile('*.ms.fits',count=nspec)
       popd
    endif
    
    if keyword_set(overplot) and xregistered('splot') then begin

       scube = rd1dspec(specname[0],datapath=datapath)
       soplot, scube.wave, scube.spec, ps=10, color='green'
       icleanup, scube, /main
       
    endif

    skywaves = [5577.339,5889.950,6300.32] ; sky wavelengths

    if n_elements(normalize) eq 1L then scale = 1.0 else scale = 1E17
    
    for i = 0L, nspec-1L do begin
    
       scube = rd1dspec(specname[i],/silent,normalize=normalize,datapath=datapath)
       header = scube.header

; interpolate over "bad" sky pixels

       skymask = scube.mask*0B
       get_element, scube.wave, skywaves, skypix
       skymask[skypix] = 1B
       skymask = smooth(float(skymask),5) gt 0B
       
       scube.spec = djs_maskinterp(scube.spec,skymask,scube.wave)
       
; smooth the spectrum

       if keyword_set(nsmooth) then begin
          scube.spec = smooth(scube.spec,nsmooth)
          scube.sigspec = smooth(scube.sigspec,nsmooth)
          scube.sky = smooth(scube.sky,nsmooth)
       endif

; grab plot information

       plotinfo, specname[i], header, xtitle=xtitle, ytitle=ytitle, $
         legendbox=legendbox, psname=psname

; S/N statistics

       snr = scube.spec/scube.sigspec
       snrmean = mean(snr)
       snrmedian = median(snr)
       snrsig = stddev(snr)

       snrstr = 'S/N = '+string(snrmedian,format='(F5.1)')
;      snrstr = 'S/N = '+string(snrmedian,format='(F5.1)')+' ('+string(snrmean,format='(F5.1)')+' +/- '+$
;        string(snrsig,format='(F5.1)')+')'
       
       legendbox = [legendbox,snrstr]
       
       if keyword_set(postscript) then begin
          ps_open, outpath+psname, /ps_fonts;, /portrait
          device, /inches, /times, /landscape, _extra=extra;, xsize=7, ysize=7
       endif
          
       if n_elements(normalize) eq 1L then $
         ytitle = textoidl('Normalized f_{\lambda}') else $
         ytitle = textoidl('f_{\lambda} (10^{-17} '+flam_units()+')')
;      ytitle = 'Relative Flux'
       xtitle = 'Wavelength ('+angstrom()+')'

       goodpix = where((scube.wave lt 5500.0) or (scube.wave gt 5640.0)) ; avoid bad sky pixels
       yrange = scale*[min(scube.spec[goodpix]),1.1*max(scube.spec[goodpix])]
       
       if (nspec eq 1L) and (not keyword_set(postscript)) then begin ; SPLOT cannot be blocked
          splot, scube.wave, scale*scube.spec, ps=10, xtitle=xtitle, $
            ytitle=ytitle, charsize=1.8, charthick=2.0, xsty=3, ysty=3, $
            yrange=yrange, xmargin=[12,3], _extra=extra
;         soplot, scube.wave, scube.sigspec, ps=10, line=2
       endif else begin
          if (i eq 0L) and (not keyword_set(postscript)) then window, 0, xs=850, ys=600
          djs_plot, scube.wave, scale*scube.spec, ps=10, xtitle=xtitle, $
            ytitle=ytitle, charsize=1.8, charthick=2.0, xsty=3, ysty=3, $
            yrange=yrange, xmargin=[12,3], _extra=extra
;         djs_oplot, scube.wave, scube.sigspec, ps=10, line=2
       endelse
          
;      legend, legendbox, left=left, right=right, /top, box=0, $
;        charthick=2.0, charsize=2.0

       if keyword_set(postscript) then ps_close

       if (nspec gt 1L) and (i le nspec-2L) then begin
          prompt:

          print, 'Object '+strn(i,length=3)+' '+legendbox[0]+' [Options: b,g,o,r,s,q]'
          cc = strupcase(get_kbrd(1))
          case strlowcase(strcompress(cc,/remove)) of
             'b': i = i-2L ; back
             'g': begin    ; goto 
                number = ''
                read, number, prompt='Goto spectrum number (0-'+strn(nspec-1L)+'): '
                number = 0 > long(number-1L) < (nspec-2L)
                i = number
             end
             'r': i = i-1L ; redraw
             'o': begin    ; overplot another spectrum
                number = ''
                read, number, prompt='Overplot spectrum number (0-'+strn(nspec-1L)+'): '
                number = 0 > long(number) < (nspec-1L)
                j = number                
                print, 'Overplotting '+specname[j]+'.'
                ospec = rd1dspec(specname[j],/silent,normalize=normalize,datapath=datapath)
;               norm = interpol(scube.spec,scube.wave,5500.0)
;               djs_oplot, ospec.wave, norm*ospec.spec/interpol(ospec.spec,ospec.wave,5500.0), $
;                 ps=10, color='green', thick=0.5
                djs_oplot, ospec.wave, scale*ospec.spec, ps=10, color='green', thick=0.5
                goto, prompt
             end
             's': begin ; overplot the sky spectrum
                norm = interpol(scube.spec,scube.wave,5500.0)
                djs_oplot, scube.wave, scale*norm*scube.sky/interpol(scube.sky,scube.wave,5500.0), $
                  ps=10, color='red', thick=0.5
                goto, prompt
             end
             'q': return ; quit
             else: 
          endcase
       endif
        
   endfor 
    
; clean up memory

   icleanup, scube
    
return
end
