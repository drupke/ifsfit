;+
; NAME:
;       PSEUDO_CONTINUUM()
;
; PURPOSE:
;       Measure the pseudo-continuum around an emission line. 
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;       wave     - wavelength vector [NPIX]
;       flux     - spectrum flux [NPIX]
;       ferr     - error spectrum [NPIX]
;       linewave - emission line wavelengths [NLINE]
;
; OPTIONAL INPUTS:
;       llimit   - lower continuum window wavelength (default LINEWAVE
;                  minus 10 Angstrom)
;       lwidth   - width of the lower continuum window (default 15
;                  Angstrom) 
;       ulimit   - upper continuum window wavelength (default LINEWAVE
;                  plus 10 Angstrom)
;       uwidth   - width of the upper continuum window (default 15
;                  Angstrom) 
;
; KEYWORD PARAMETERS:
;       doplot   - generate a plot of the continuum fit
;
; OUTPUTS:
;       pseudo   - a [2,NLINE] array with the continuum measurement in
;                  pseudo[0,*] and the error in pseudo[1,*]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       A fixed continuum error of 10% is assigned.  All wavelengths
;       are in Angstroms. 
;
; PROCEDURE:
;       The continuum windows are defined by LLIMIT, ULIMIT, LWIDTH,
;       and UWIDTH.  The mean continuum in each window is computed and
;       a line is fitted between the points.  The continnum is
;       interpolated at the central wavelength of the emission line.
;
; EXAMPLE:
;
; PROCEDURES USED:
;       IM_SYMBOLS(), GET_ELEMENT, DJS_PLOT, DJS_ITERSTAT, 
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 October 31, U of A
;-

function pseudo_continuum, wave, flux, ferr, linewave, llimit=llimit, $
  lwidth=lwidth, ulimit=ulimit, uwidth=uwidth, doplot=doplot

    npix = n_elements(wave)
    if npix eq 0L then begin
       splog, 'No wavelength vector defined!'
       return, -1
    endif

    im_symbols, 'circle', /fill, psize=2

    nline = n_elements(linewave)

    if n_elements(llimit) eq 0L then llimit = linewave-10.0
    if n_elements(ulimit) eq 0L then ulimit = linewave+10.0

    if n_elements(lwidth) eq 0L then lwidth = replicate(15.0,nline)
    if n_elements(uwidth) eq 0L then uwidth = replicate(15.0,nline)

    if nline gt 1L then begin

       for k = 0L, nline-1L do begin
          pseudo1 = pseudo_continuum(wave,flux,ferr,linewave[k],llimit=llimit[k],$
            lwidth=lwidth[k],ulimit=ulimit[k],uwidth=uwidth[k],doplot=doplot)
          if k eq 0L then pseudo = pseudo1 else pseudo = [ [pseudo], [pseudo1] ]
       endfor

       return, pseudo

    endif
    
    if n_elements(cmethod) eq 0L then cmethod = 0L

    npix = n_elements(wave)
    
    get_element, wave, [llimit,ulimit]+[-20,+20], ww ; plotting range
    
    if keyword_set(doplot) then begin

       djs_plot, wave[(ww[0]-20)>0L:(ww[1]+20)<(npix-1)], flux[(ww[0]-20)>0L:(ww[1]+20)<(npix-1)], $
         color='yellow', xsty=3, ysty=3, ps=10, xtitle='Wavelength ('+angstrom()+')', $
         ytitle='Flux', charsize=1.5, charthick=2.0

    endif

    get_element, wave, llimit+[-lwidth,+lwidth]/2.0, aa
    get_element, wave, ulimit+[-uwidth,+uwidth]/2.0, bb

    localwave = wave[aa[0]:bb[1]]
    
    djs_iterstat, flux[aa[0]:aa[1]], sigrej=3.0, mean=lpoint, mask=lmask
    djs_iterstat, flux[bb[0]:bb[1]], sigrej=3.0, mean=rpoint, mask=rmask

    lwave = djs_mean(wave[aa[0]:aa[1]])
    rwave = djs_mean(wave[bb[0]:bb[1]])

    slope = (rpoint-lpoint)/(rwave-lwave)
    coeff = [rpoint-slope*rwave,slope]
;   coeff = linfit([lwave,rwave],[lpoint,rpoint],/double)
    fit = poly(localwave,coeff)

    pseudo = interpol(fit,localwave,linewave)*[1.0,0.10] ; 10% error
    
    if keyword_set(doplot) then begin
       
       polyfill, [wave[aa[0]],wave[aa[1]],wave[aa[1]],wave[aa[0]]], $ ; blue
         [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
         /line_fill, orientation=45, color=djs_icolor('blue')
       
       polyfill, [wave[bb[0]],wave[bb[1]],wave[bb[1]],wave[bb[0]]], $ ; red
         [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
         /line_fill, orientation=45, color=djs_icolor('red')

       plots, lwave, lpoint, ps=8, color=djs_icolor('white')
       plots, rwave, rpoint, ps=8, color=djs_icolor('white')
       djs_oplot, localwave, fit, thick=2.0, color='cyan'
       plots, linewave, pseudo[0], ps=8, color=djs_icolor('white')

       cc = get_kbrd(1)
       
    endif

return, pseudo
end
