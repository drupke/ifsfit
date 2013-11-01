;+
; NAME:
;	IREDSHIFT()
;
; PURPOSE:
;	Derive the redshift of a galaxy by using PCA eigentemplatesand
;	chi-squared fitting.  
;
; CALLING SEQUENCE:
;	zans = iredshift(speclist,[datapath=,zmin=,zmax=,npoly=], $
;		[obsname=],doplot=doplot,update=update,_extra=extra)
;
; INPUTS:
;	speclist - FITS list of spectra
;
; OPTIONAL INPUTS:
;	datapath - path to the data and the write directory
;	zmin     - minimum redshift to consider (can be negative for
;                  blushifts) [default -0.01] 
;	zmax     - maximum redshift to consider [default 10000 km/s]
;	npoly    - number of polynomial "eigenspectra" to append to
;                  the PCA eigentemplates
;	obsname  - name of the observatory when computing the
;                  heliocentric correction (default: KPNO)
;	extra    - extra inputs for ZCOMPUTE()
;	
; KEYWORD PARAMETERS:
;	doplot   - generate a plot of the spectrum and best-fitting
;                  linear combination of eigentemplates
;	write    - write the updated header and spectra to disk
;
; OUTPUTS:
;	zans     - redshift computation output from ZCOMPUTE()
;
; OPTIONAL OUTPUTS:
;	zinfo    - redshift output structure for all input spectra 
;
; COMMENTS:
;	The heliocentric correction to the redshift is placed in the
;	header, but not returned to the user.
;
; EXAMPLE:
;
; PROCEDURES USED:
;	RD1DSPEC(), COMBINE1FIBER, ZCOMPUTE(), DJS_MEAN(), SXPAR(),
;	HMS2DEC(), OBSERVATORY, HELIOCENTRIC(), WRITE1DSPEC,
;	ICLEANUP 
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2002 January 20, U of A
;-

function iredshift, speclist, datapath=datapath, zmin=zmin, zmax=zmax, $
                     npoly=npoly, obsname=obsname, doplot=doplot, update=update, $
                     eigendir=eigendir, eigenfile=eigenfile, columns=columns, _extra=extra

    nspec = n_elements(speclist)
    if nspec eq 0L then begin
       print, 'Syntax - zans = iredshift(speclist,[datapath=,zmin=,zmax=] $'
       print, '   [npoly=,obsname=],doplot=doplot,update=update,_extra=extra'
       return, -1
    endif

    if not keyword_set(datapath) then datapath = cwd()
    if n_elements(npoly) eq 0L then npoly = 2L

; redshift search radius

    light = 2.99792458D5 ; speed of light [km/s]
    czmax = 15000.0D     ; Hubble flow velocity [km/s]

    if n_elements(zmin) eq 0L then zmin = -0.01
    if n_elements(zmax) eq 0L then zmax = +czmax/light
    
    if n_elements(obsname) eq 0L then obsname = 'KPNO' ; observatory structure
    observatory, obsname, obs 

; if multiple spectra have been passed, call this routine recursively

    if nspec gt 1L then begin
       t0 = systime(/seconds)
       for k = 0L, nspec-1L do begin

          zans1 = iredshift(speclist[k],datapath=datapath,zmin=zmin,zmax=zmax,$
                             npoly=npoly,obsname=obsname,doplot=doplot,$
                             update=update,eigendir=eigendir, eigenfile=eigenfile,$
                             columns=columns,_extra=extra)
          if k eq 0L then zans = zans1 else zans = [[zans],[zans1]]
          splog, 'Object #', k, '  Elap time=', systime(/seconds)-t0, $
            ' (sec)  z=', zans1[0].z, ' (pix)'
          
       endfor
       return, zans
    endif

; determine the logarithmic wavelength spacing of the eigentemplates

    if n_elements(eigendir) eq 0L then eigendir = $
      concat_dir(getenv('IDLSPEC2D_DIR'),'templates')
    if n_elements(eigenfile) eq 0L then eigenfile = 'spEigenGal-*'
    allfiles = findfile(djs_filepath(eigenfile,root_dir=eigendir),count=ct)
    if (ct eq 0) then message, 'Unable to find EIGENFILE matching '+eigenfile+'.'
    thisfile = allfiles[(reverse(sort(allfiles)))[0]]

    eigenflux = readfits(thisfile,shdr)
    naxis = sxpar(shdr,'NAXIS1')
    coeff0 = sxpar(shdr,'COEFF0')
    coeff1 = sxpar(shdr,'COEFF1')
    
    if keyword_set(npoly) then eigenflux = [[eigenflux],[poly_array(naxis,npoly)]]
    
    lam = 10D^(coeff0+coeff1*findgen(naxis)) ; logarithmic wavelength vector
    newloglam = alog10(lam)
          
    scube = rd1dspec(speclist,datapath=datapath) ; read in the galaxy spectrum

    objflux = scube.spec        ; object flux [erg/s/cm^2/A]
    objsig = scube.sigspec      ; object sigma
    objivar = 1.0/objsig^2.0    ; object inverse variance
    wave = scube.wave           ; wavelength [A]
    logwave = alog10(wave)

    header = scube.header
    hdr = header
    sxaddpar, hdr, 'COEFF0', coeff0
    sxaddpar, hdr, 'COEFF1', coeff1
       
; resample the galaxy spectrum into log-lamba, at the spacing of the
; eigenflux log-lambda vector
    
    combine1fiber, logwave, objflux, objivar, newloglam=newloglam, $
      newflux=newflux, newivar=newivar

    fluxnorm = double(djs_mean(newflux)) ; normalization constants
    ivarnorm = double(djs_mean(newivar))

    zans = zfind(newflux/fluxnorm,newivar/ivarnorm,hdr=hdr,$
                 eigenfile=eigenfile,eigendir=eigendir,npoly=npoly,$
                 zmin=zmin,zmax=zmax,doplot=doplot,_extra=extra)

; compute the heliocentric correction
       
    ra = 15D*hms2dec(sxpar(header,'RA'))
    dec = hms2dec(sxpar(header,'DEC'))
    epoch = sxpar(header,'EPOCH')
    jd = sxpar(header,'JD')

    vcorr = heliocentric(ra,dec,epoch,jd=jd,longitude=(360-obs.longitude),$
                         latitude=obs.latitude,altitude=obs.altitude)

    z_helio = zans[0].z + vcorr[0]/light

; compute the velocity dispersion

;   vdans = vdispfit(newflux/fluxnorm,newivar/ivarnorm,hdr=hdr,$
;                    newloglam,npoly=npoly,zobj=zans[0,*].z,eigendir=eigendir,$
;                    eigenfile='spEigenElodie.fits',$
;                    columns=lindgen(24),yfit=dispflux,/doplot)

;  save, newloglam, newflux, dispflux, fluxnorm, filename='test.idlsave'
    
; update the header

    newhead = header
    sxaddpar, newhead, 'Z', z_helio, ' heliocentric redshift computed '+im_today(), $
      before='HISTORY', format='(F17.9)'
    sxaddpar, newhead, 'Z_ERR', zans.z_err, ' heliocentric redshift error', $
      before='HISTORY', format='(F17.9)'

    if keyword_set(update) then begin
       print, 'Updating the header for '+speclist+'.'
       djs_modfits, scube.datapath+scube.specname, 0, newhead
    endif

; plot the spectrum and the best fit.  this plot doesn't work if the
; columns keyword is used

;    if keyword_set(doplot) then begin
;
;       window, 2, xs=450, ys=450
;
;       if n_elements(columns) ne 0L then neigen = n_elements(columns) else $
;         neigen = (size(eigenflux,/dimension))[1]
;       bestspec = newflux*0.0
;       for j = 0L, neigen-1L do bestspec = bestspec + zans.theta[j]*eigenflux[*,j]
;       
;       djs_plot, lam/(1+(zans.z+vcorr/light)), newflux/fluxnorm, ps=10, $
;         xsty=3, ysty=3, charsize=1.5, charthick=2.0, xthick=2.0, $
;         ythick=2.0, xtitle='Wavelength', ytitle='Relative Flux', $
;         color='red'
;       djs_oplot, lam, bestspec, ps=10, line=0, color='green'
;
;    endif

; cleanup memory
    
    icleanup, scube
    
return, zans
end
