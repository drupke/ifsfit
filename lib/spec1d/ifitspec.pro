;+
; NAME:
;       IFITSPEC()
;
; PURPOSE:
;	Fit a galaxy spectrum, including the continuum,
;	absorption-line and the emission-line spectrum. 
;
; CALLING SEQUENCE:
;       ifitspec, flux, wave, eigenflux, eigenwave, $
;          [invvar=,zobj=,specres=,eigenres=,snrcut=,vmaxshift=,$
;          vminwave=,vmaxwave=,maskwidth=,linepars=abslines=,$
;          zindex=,windex=,findex=,fvalue=,backfit=,linefit=,$
;          absfit=,indices=,speclinefit=,speclineewfit=,starflux=],$
;          /doplot,/nozupdate,/combine_blends,_extra=extra
;
; INPUTS:
;       flux      - input spectrum [erg/s/cm2/Angstrom] [NPIX]
;       wave      - corresponding wavelength vector [Angstrom] [NPIX] 
;       eigenflux - eigen-template flux vector [NEIGENPIX,NSTAR]
;       eigenwave - corresponding wavelength vector [NEIGENPIX]
;       linepars  - input emission line structure array as read by 
;                   READ_LINEPARS() [NLINE]
;          line   - emission line name (string)
;          wave   - central wavelength (air wavelength) [Angstrom]
;          zindex - see ILINEFIT() documentation
;          windex - see ILINEFIT() documentation
;          findex - see ILINEFIT() documentation
;          fvalue - see ILINEFIT() documentation
;          blend  - lines with the same "blend" string value are
;                   considered blended in the spectrum; see the
;                   PROCEDURE section for what this means in practice 
;
; OPTIONAL INPUTS:
;       invvar    - inverse variance spectrum corresponding to FLUX
;                   [NPIX]; if not provided then the output errors
;                   will not be correct
;       zobj      - approximate galaxy redshift (default 0.0)
;       specres   - FWHM spectral resolution of FLUX [Angstrom];  can
;                   be either a scalar or an [NPIX] vector (default 5) 
;       eigenres  - scalar FWHM spectral resolution of EIGENFLUX
;                   [NEIGENPIX] (default 3 Angstrom)
;       snrcut    - compute upper limits on lines with S/N < SNRCUT
;                   (default 3.0)
;       vmaxshift - do not allow the continuum redshift to change
;                   (through cross-correlation) by more than
;                   +/-VMAXSHIFT [km/s] (default 100.0, see PROCEDURE) 
;       vminwave  - compute the cross-correlation only between
;                   VMINWAVE and VMAXWAVE [default 4200.0>min(WAVE)]
;       vmaxwave  - compute the cross-correlation only between
;                   VMINWAVE and VMAXWAVE [default 5800.0<max(WAVE)]
;       maskwidth - when masking emission lines, mask all pixels
;                   within MASKWIDTH Angstroms of the line center
;                   (default 30.0)
;       zindex    - overwrite LINEPARS.ZINDEX on second iteration 
;       windex    - overwrite LINEPARS.WINDEX on second iteration
;       findex    - overwrite LINEPARS.FINDEX on second iteration
;       fvalue    - overwrite LINEPARS.FVALUE on second iteration
;	extra     - extra parameters for IBACKFIT() and ILINEFIT() 
;	
; KEYWORD PARAMETERS:
;       doplot    - plot the individual fits of each absorption and
;                   emission line and wait for a keystroke (note: the
;                   DOPLOT option in IM_ZTWEAK() and SPECTRAL_INDICES()
;                   have been suppressed in all cases)
;       nozupdate - if set then do not improve the redshift of the
;                   absorption-line spectrum using the best-fitting 
;                   continuum (see PROCEDURE)
;       combine_blends - combine the individual [O II] doublet lines
;                        into a single emission line
;
; OUTPUTS:
;       All outputs are optional.  For example, the fitting can be
;       carried out and examined with DOPLOT.
;
; OPTIONAL OUTPUTS:
;       backfit       - see IBACKFIT()
;       linefit       - see ILINEFIT()
;       absfit        - see IBACKFIT()
;       indices       - see SPECTRAL_INDICES()
;       speclinefit   - fit to the pure emission-line spectrum  
;                       returned by ILINEFIT() [NPIX]
;       speclineewfit - fit to the emission lines in the data spectrum
;                       to determine the emission-line EW
;       starflux      - identical to EIGENFLUX but broadened and
;                       wavelength resampled to match FLUX and WAVE
;                       [NPIX,NSTAR]
;
; PROCEDURE:
;       The essential pieces of this routine are a spectrum, an
;       arbitrary number of template spectra, and a data structure
;       (LINEPARS) specifying the emission lines to fit.  IFITSPEC()
;       matches the wavelength range and spectral resolution of the
;       templates to the spectrum [IMATCH_TEMPLATES()] and finds the
;       best-fitting linear combination of templates for the spectrum
;       [IBACKFIT()] based on an initial guess of the redshift of the
;       galaxy.  A set of absorption lines (ABSLINES) in the
;       best-fitting continuum spectrum can also be fitted
;       [IABSLINEFIT()].
;
;       The main iteration loop proceeds in the following way: [1]
;       Match the templates (in wavelength spacing and resolution) to
;       the data [IMATCH_TEMPLATES()]; [2] Mask pixels near emission
;       lines [EMISSION_MASK()]; [3] Fit the continuum [IBACKFIT()];
;       [4] Unless NOZUPDATE=1, improve the redshift of the continuum
;       by cross-correlating the best-fitting continuum spectrum with
;       the data [IM_ZTWEAK()]; [5] Subtract the continuum spectrum
;       and fit the pure emission-line spectrum [ILINEFIT()]; [6]
;       Steps [1-5] are iterated a second time to improve the fitting;
;       [7] Fit the absorption lines specified in ABSLINES
;       [IABSLINEFIT()]; [8] Measure color and Lick indices in the
;       input spectrum [SPECTRAL_INDICES()]; [9] Measure Balmer- and
;       forbidden-line equivalent widths using the results of the
;       simultaneous Gaussian fitting to constrain the line profiles.
;
;       Upper limits on emission lines are also computed in 
;       IUPPER_LIMITS().  
;
;       In the emission-line fitting we identify three cases:
;       [1] well-measured line, [2] upper limits, and [3] un-measured
;       lines.  The error codes for cases [1] and [2] are detailed
;       below:
;
;       Case [2]: Upper limits are computed on lines that were either
;       dropped from the fit by MPFIT [see ILINEFIT()], or they have
;       S/N < SNRCUT.  They are identified by the following error
;       codes and quantities:
;
;          LINESIGMA         =  either tied to a well-detected line,
;                               or equal to LINERES at line center 
;          LINESIGMA_ERR     = -3.0
;
;          LINEAREA          =  computed in IUPPER_LIMITS()
;          LINEAREA_ERR      = -3.0
;
;          LINEBOX           =  not meaningful 
;          LINEBOX_ERR       =  not meaningful 
;
;          LINECONTLEVEL     =  computed in IUPPER_LIMITS()
;          LINECONTLEVEL_ERR = -3.0
;
;          LINEEW_AREA       =  based on LINEAREA and LINECONTLEVEL 
;          LINEEW_AREA_ERR   = -3.0
;
;          LINEEW_BOX        =  not meaningful
;          LINEEW_BOX_ERR    =  not meaningful
;
;          LINENPIX          =  not meaningful
;          LINEDOF           =  not meaningful
;          LINECHI2          = -3.0
;
;       Case [3]: Un-measured lines are either completely or partially
;       out of the wavelength range of the spectrum.  They have the
;       following (equivalent) error codes:
;
;          LINEZ             =  ZOBJ
;          LINEZ_ERR         = -2.0
;          LINESIGMA         =  0.0
;          LINESIGMA_ERR     = -2.0
;          LINEAREA          =  0.0
;          LINEAREA_ERR      = -2.0
;          LINEBOX           =  0.0
;          LINEBOX_ERR       = -2.0
;          LINECONTLEVEL     =  0.0
;          LINECONTLEVEL_ERR = -2.0
;          LINEEW_AREA       =  0.0
;          LINEEW_AREA_ERR   = -2.0
;          LINEEW_BOX        =  0.0
;          LINEEW_BOX_ERR    = -2.0
;          LINENPIX          =    0
;          LINEDOF           =  0.0
;          LINECHI2          = -2.0
;
; COMMENTS:
;       The dispersion of the eigentemplates [Angstrom/pixel] is
;       assumed to be constant.  
; 
;       The cross-correlation code for improving the radial velocity
;       of the continuum is not general.  Currently we use the
;       wavelength range [4200,5800] Angstrom.
;
;       Whereas SPECRES can be either a scalar or an NPIX vector,
;       EIGENRES must be a scalar.  A variation in the spectral
;       resolution of the templates as a function of wavelength is
;       currently not implemented.
;
;       The Lick indices are uncorrected for emission.
;
; TODO:
;       [1] Optionally allow EIGENRES to be a vector.
;
; EXAMPLES:
;
;
; INTERNAL SUPPORT ROUTINES:
;       BALMERINDX(), IFITSPEC_COMBINE_BLENDS()
;
; PROCEDURES USED:
;       SPLOG, IM_NORMALIZE(), IMATCH_TEMPLATES(), EMISSION_MASK(),
;       IBACKFIT(), IM_ZTWEAK(), STRUCT_ADDTAGS(), ILINEFIT(),
;       IABSLINEFIT(), SPECTRAL_INDICES(), MEASURE_LINEEW(),
;       FILL_EW(), REMOVE, MATCH_STRING()
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002-2003, U of A
;       jm03jul28uofa - added IFITSPEC_COMBINE_BLENDS() sub-routine
;                       and COMBINE_BLENDS keyword to properly treat
;                       the [O II] doublet 
;-

function balmerindx, linewave, blend, forbidden=forbidden
; index non-blended Balmer lines
    
    nline = n_elements(linewave)
    forbidden = lindgen(nline)   ; default is to assume no Balmer lines 
    
; reject non-blended Balmer lines.  find the Balmer lines by matching
; the wavelength LINEWAVE to the nearest tenth of an Angstrom

    balmerwaves = [3889.049,3970.072,4101.734,4340.464,4861.325,6562.80] ; NOT GENERAL!
    nbalmer = n_elements(balmerwaves)
    balmerindx = lonarr(nbalmer)
    
; BALMERINDX indexes all Balmer lines in LINEWAVE.  if no Balmer lines
; were fitted then return
    
    for i = 0L, nbalmer-1L do balmerindx[i] = where(long(10*balmerwaves[i]) eq long(10*linewave))

    pos = where(balmerindx ne -1L,nbalmerindx,comp=neg)
    if nbalmerindx ne 0L then balmerindx = balmerindx[pos] else return, -1L

; now remove from BALMERINDX any blended Balmer lines (e.g., H-alpha) 
    
    rem = lonarr(nbalmerindx)-1L
    for i = 0L, nbalmerindx-1L do $
      if total(strmatch(blend,blend[balmerindx[i]])) gt float(1.0) then $
      rem[i] = balmerindx[i]

; if the only Balmer line in BALMERINDX is a blend then return
    
    match = where(rem ne -1L,nmatch)
    if (nmatch eq 1L) and (nbalmerindx eq 1L) then return, -1L
    if nmatch ne 0L then remove, match, balmerindx

    nbalmerindx = n_elements(balmerindx)
    
; finally, remove all the non-Balmer (forbidden) lines from the fully
; index linelist to generate the FORBIDDEN variable.  if NBALMERINDX
; equals NLINE that means that no forbidden lines were successfully
; fitted 

    if (nbalmerindx lt nline) then $
      remove, forbidden[balmerindx], forbidden else $
      forbidden = -1L

return, balmerindx
end

function ifitspec_combine_blends, linefit
; only [O II] is treated properly in this sub-routine

    allnames = linefit.linename
    indx = where(strmatch(allnames,'*OII_3727*') eq 1B,nindx,comp=keep,ncomp=nkeep)

    if nindx ne 2L then begin
       splog, '[O II] doublet line has not been measured.'
       return, linefit
    endif
    
    s = linefit[indx]
    s1 = s[0]
    s2 = s[1]

; initialize the output structure
    
    result = icreate_linestruct(1)

    errflag = 0.0
    if (s1.linearea_err eq -1.0) or (s2.linearea_err eq -1.0) then errflag = -1.0 ; "dropped"
    if (s1.linearea_err eq -2.0) or (s2.linearea_err eq -2.0) then errflag = -2.0 ; "not measured"

    case errflag of

       0.0: begin ; both lines are well-measured

; linearly add the redshift and add the error in quadrature
          result.linez = djs_mean(s.linez)
          result.linez_err = sqrt(s1.linez_err^2 + s2.linez_err^2)

; add the sigma line widths and errors in quadrature
          result.linesigma = sqrt(total(s.linesigma^2))
          result.linesigma_err = sqrt(total(s.linesigma_err^2))

; linearly add the Gaussian flux and add the error in quadrature
          result.linearea = s1.linearea + s2.linearea
          result.linearea_err = sqrt(s1.linearea_err^2 + s2.linearea_err^2)
          
; linearly add the box flux and add the error in quadrature
          result.linebox = s1.linebox + s2.linebox
          result.linebox_err = sqrt(s1.linebox_err^2 + s2.linebox_err^2)

; add the number of pixels linearly (not strictly correct)
          result.linenpix = total(s.linenpix)

; compute the mean number of degrees of freedom (not strictly correct)
          result.linedof = djs_mean(s.linedof)

; assign the chi2 to be the mean chi2
          result.linechi2 = djs_mean(s.linechi2)

       end

       -1.0: begin ; one line was dropped from the fit by MPFIT

; linearly add the redshift and add the error in quadrature
          result.linez = djs_mean(s.linez)
          result.linez_err = sqrt(s1.linez_err^2 + s2.linez_err^2)

; add the sigma line widths in quadrature
          result.linesigma = sqrt(total(s.linesigma^2))
          result.linesigma_err = -1.0

; linearly add the Gaussian flux and add the error in quadrature
          result.linearea = 0.0
          result.linearea_err = -1.0
          
; linearly add the box flux and add the error in quadrature
          result.linebox = 0.0
          result.linebox_err = -1.0

; add the number of pixels linearly (not strictly correct)
          result.linenpix = total(s.linenpix)

; compute the mean number of degrees of freedom (not strictly correct)
          result.linedof = djs_mean(s.linedof)

; assign the chi2 to be the mean chi2
          result.linechi2 = djs_mean(s.linechi2)

       end
       
       -2.0: begin              ; set the combined line as "not measured"

; simply copy the "unmeasured" line structure into RESULT          
          
          nm = where(s.linearea_err eq -2.0)
          result = s[nm[0]]
          
       end

    endcase

; fill the output structure with mean quantities 

    result.linename = 'OII_3727'
    result.linewave = djs_mean(s.linewave)
    result.line_blend = 'good3727'
    
; concatenate RESULT to the input structure and return
    
    if nkeep ne 0L then $
      sout = struct_append(result,linefit[keep]) else $
      sout = result

return, sout
end

pro ifitspec, flux, wave, eigenflux, eigenwave, invvar=invvar, $
  zobj=zobj, specres=specres, eigenres=eigenres, snrcut=snrcut, $
  vmaxshift=vmaxshift, vminwave=vminwave, vmaxwave=vmaxwave, $
  maskwidth=maskwidth, linepars=linepars, abslines=abslines, $
  zindex=zindex, windex=windex, findex=findex, fvalue=fvalue, $
  backfit=backfit, linefit=linefit, absfit=absfit, indices=indices, $
  speclinefit=speclinefit, speclineewfit=speclineewfit, $
  starflux=starflux, doplot=doplot, nozupdate=nozupdate, $
  combine_blends=combine_blends, _extra=extra

    light = 2.99792458D5        ; speed of light [km/s]

; defaults and error checking

    nwave = n_elements(wave)
    npix = n_elements(flux)
    neigenwave = n_elements(eigenwave)

    if (nwave eq 0L) or (npix eq 0L) or (neigenwave eq 0L) or (n_elements(eigenflux) eq 0L) then begin
       print, 'Syntax - ifitspec, flux, wave, eigenflux, eigenwave, $'
       print, '   [invvar=,zobj=,specres=,eigenres=,snrcut=,vmaxshift=,$'
       print, '   vminwave=,vmaxwave=,maskwidth=,linepars=abslines=,$'
       print, '   zindex=,windex=,findex=,fvalue=,backfit=,linefit=,$'
       print, '   absfit=,indices=,speclinefit=,speclineewfit=,starflux=],$'
       print, '   /doplot,/nozupdate,/combine_blends,_extra=extra'
       return
    endif

    if (nwave ne npix) then begin
       splog, 'WAVE and FLUX do not have the same number of elements.'
       return
    endif

    eigen_ndim = size(eigenflux,/n_dimension)
    if eigen_ndim ne 2L then begin
       splog, 'EIGENFLUX must be a two-dimensional array.'
       return
    endif

    eigendim = size(eigenflux,/dimension)
    nstar = eigendim[1]

    if (neigenwave ne eigendim[0]) then begin
       splog, 'EIGENWAVE and EIGENFLUX do not have the same number of elements.'
       return
    endif

    nline = n_elements(linepars)
    if nline eq 0L then begin
       splog, 'Structure LINEPARS must be passed.'
       return
    endif
    
    if n_elements(zobj) eq 0L then zobj = float(0.0) else zobj = float(zobj)
    
    nivar = n_elements(invvar)
    if nivar eq 0L then begin
       splog, 'WARNING:  INVVAR not passed - errors will be incorrect.'
       invvar = make_array(npix,value=1.0,/float) 
    endif else begin
       if nivar ne npix then begin
          splog, 'FLUX and INVVAR do not have the same number of elements.'
          return
       endif
    endelse

; make specres an NPIX element vector if it is a scalar
    
    if n_elements(specres) eq 0L then specres = replicate(5.0,npix)     ; FWHM [Angstrom]
    if n_elements(specres) eq 1L then specres = replicate(specres,npix)

    if n_elements(specres) ne npix then begin
       splog, 'WARNING: SPECRES does not have the correct number of elements.'
       return
    endif 

    nspecres = n_elements(specres)
    
    neigenres = n_elements(eigenres)
    if neigenres eq 0L then eigenres = 3.0 ; FWHM [Angstrom]

    if n_elements(snrcut) eq 0L then snrcut = 3.0

    if n_elements(vmaxshift) eq 0L then vmaxshift = 100.0 ; [km/s]
    if n_elements(vminwave) eq 0L then vminwave = 4200.0>min(wave) ; [Angstrom]
    if n_elements(vmaxwave) eq 0L then vmaxwave = 5800.0<max(wave) ; [Angstrom]

    if n_elements(maskwidth) eq 0L then maskwidth = 30.0  ; [Angstrom]

; normalize FLUX to its maximum value; only normalize INVVAR if it was 
; passed by the user
    
    objflux = im_normalize(flux,/max,const=normconst) 
    objwave = wave

    var = 1.0/invvar
    if nivar eq 0L then begin
       objivar = invvar 
       objsig = invvar
    endif else begin
       objivar = invvar*normconst^2.0
       objsig = sqrt(var)/normconst
    endelse

; use 105 km/s as an initial guess for the intrinsic Gaussian line
; width for all the emission lines (the result is not very sensitive
; to this guess, but helps speed up convergence) 

    sigguess = replicate(105.0,nline) / light / alog(10.0) ; log-Angstrom units

; convert the Gaussian sigma line width of the emission and absorption
; lines due to instrumental broadening to km/s

    fwhm2sig = 2.0*sqrt(2.0*alog(2.0))

    lineres1 = interpol(specres,objwave,linepars.wave*(1+zobj))       ; [Angstrom]
    elineres = light * lineres1 / (linepars.wave*(1+zobj)) / fwhm2sig ; sigma width [km/s]

    if (n_elements(abslines) ne 0L) then begin
       lineres1 = interpol(specres,objwave/(1+zobj),abslines.wave) ; [Angstrom]
       alineres = light * lineres1 / abslines.wave / fwhm2sig      ; sigma width [km/s]
    endif

; initialize variables for the iterative fitting

    itermax = 1L
    influx = objflux

    ziter = float(zobj)  ; absorption-line redshift guess [updated by IM_ZTWEAK()]
    zguess = float(zobj) ; emission-line redshift [overwritten by ILINEFIT()]

    zabs = float(0.0)      ; default absorption-line redshift and error
    zabs_err = float(-1.0)

; ------------------------------------------------------------
; begin main iteration loop    
; ------------------------------------------------------------
    
    for iter = 0L, itermax do begin 

       restwave = objwave/(1.0+ziter) ; rest wavelength

; prepare the eigentemplates for fitting

       splog, 'Re-sampling and broadening the eigentemplates.'
       starflux = imatch_templates(eigenflux,eigenwave,restwave,$
          starres=eigenres,newres=specres)

; mask emission-line wavelengths       

       mask = emission_mask(restwave,width=maskwidth,/telluric)

; fit the continuum spectrum

       t0 = systime(1)
       backfit = ibackfit(influx,restwave,starflux,invvar=objivar,mask=mask,_extra=extra)
       splog, format='("CPU time for continuum fitting on iteration ",I1,"/",I1,": ",G0," '+$
         'seconds.")', i, itermax, systime(1)-t0
       
       continuum = backfit.continuum
       continuum_sigma = backfit.continuum_sigma

; improve the redshift guess by cross-correlating the best-fitting
; continuum spectrum with the data

       if not keyword_set(nozupdate) then begin

          zupdate = im_ztweak(influx,restwave,continuum,restwave,$
            specivar=objivar,mask=mask,vmin=-vmaxshift,vmax=+vmaxshift,$
            minwave=vminwave,maxwave=vmaxwave,/silent,doplot=0)
          
          if (zupdate.errflag eq 0L) then begin
             splog, 'Cross-correlation converged . . . updating the absorption-line redshift.'
             zabs = float(ziter + zupdate.zshift) ; absorption-line redshift and error
             zabs_err = float(zupdate.zshift_err)
             ziter = float(zabs)
          endif

       endif 

; append the absorption-line redshift to the continuum fit structure
       
       backfit = struct_addtags(backfit,{z_abs: zabs,z_abs_err:zabs_err})
       
; fit the emission-line spectrum; on the second iteration allow the
; emission line constraints to overwrite those specified in LINEPARS
          
       espectrum = objflux-continuum 
       espectrum_err = sqrt(objsig^2.0 + continuum_sigma^2.0) ; quadrature
       espectrum_ivar = 1.0/espectrum_err^2.0

       if iter eq itermax then begin

          if n_elements(zindex) eq nline then zindex1 = zindex
          if n_elements(windex) eq nline then windex1 = windex
          if n_elements(findex) eq nline then findex1 = findex
          if n_elements(fvalue) eq nline then fvalue1 = fvalue

       endif else begin

          zindex1 = linepars.zindex
          windex1 = linepars.windex
          findex1 = linepars.findex
          fvalue1 = linepars.fvalue
          
       endelse

;      struct_print, linepars

       t0 = systime(1)
       splog, 'Fitting the emission-line spectrum.'
       linefit = ilinefit(espectrum,objwave,linepars.wave,elineres,$
         invvar=espectrum_ivar,linename=linepars.line,zindex=zindex1,$
         windex=windex1,findex=findex1,fvalue=fvalue1,zguess=zguess,$
         sigguess=sigguess,specfit=speclinefit,bfit=bfit,_extra=extra)
       splog, format='("CPU time for line fitting on iteration ",I1,"/",I1,": ",'+$
         'G0," seconds.")', i, itermax, systime(1)-t0

; initialize input variables for the second iteration.  use the
; emission-line redshifts and sigma widths as initial guesses for the
; second iteration, but only for well detected lines

       influx = influx - speclinefit

       goodlines = where(linefit.linearea_err gt 0.0,ngoodlines)
       if ngoodlines ne 0L then begin
          sigguess = linefit.linesigma/alog(10.0)/light ; log-Angstrom
          zguess = float(linefit.linez)
       endif

    endfor 

; ------------------------------------------------------------
; end main iteration loop    
; ------------------------------------------------------------
    
; fit the absorption lines.  note that we're using the error in the
; original data as the error in the continuum because the continuum
; error returned by IBACKFIT() is probably an underestimate of the
; true error in the continuum  

    if n_elements(abslines) ne 0L then begin
    
       absfit = iabslinefit(restwave,continuum,objsig,abslines,$
         lineres=alineres,doplot=doplot,_extra=extra)
       
    endif
    
; put back the normalization constant into the spectral fits and
; the absorption and emission-line fluxes

    speclinefit = speclinefit*normconst
    speclineewfit = speclinefit*0.0 ; initialize this variable

    bfit = bfit*normconst
    backfit.continuum = backfit.continuum*normconst
    backfit.continuum_sigma = backfit.continuum_sigma*normconst

; Gaussian fluxes
    
    goodfit = where(linefit.linearea_err gt 0.0,ngoodfit)
    if ngoodfit ne 0L then begin
       linefit[goodfit].linearea = linefit[goodfit].linearea*normconst
       linefit[goodfit].linearea_err = linefit[goodfit].linearea_err*normconst
    endif

; box fluxes

    goodfit = where(linefit.linebox_err gt 0.0,ngoodfit)
    if ngoodfit ne 0L then begin
       linefit[goodfit].linebox = linefit[goodfit].linebox*normconst
       linefit[goodfit].linebox_err = linefit[goodfit].linebox_err*normconst
    endif

;;; for "dropped" lines set the velocity width equal to the spectral
;;; resolution with zero error    
;;    
;;    setsigma = where(linefit.linearea_err eq -1.0,nsetsigma)
;;    if nsetsigma ne 0L then begin
;;       setwave = linefit[setsigma].linewave
;;       linefit[setsigma].linesigma = light * interpol(specres,restwave,setwave) / setwave / fwhm2sig
;;       linefit[setsigma].linesigma_err = 0.0
;;    endif

    if n_elements(absfit) ne 0L then begin
       
       goodfit = where(absfit.linearea_err gt 0.0,ngoodfit)
       if ngoodfit ne 0L then begin
          absfit[goodfit].linecontlevel     = absfit[goodfit].linecontlevel*normconst
          absfit[goodfit].linecontlevel_err = absfit[goodfit].linecontlevel_err*normconst
          absfit[goodfit].linearea          = absfit[goodfit].linearea*normconst
          absfit[goodfit].linearea_err      = absfit[goodfit].linearea_err*normconst
          absfit[goodfit].linebox           = absfit[goodfit].linebox*normconst
          absfit[goodfit].linebox_err       = absfit[goodfit].linebox_err*normconst
       endif

    endif
       
    linefit.line_blend = linepars.blend ; blended lines

; combine the [O II] doublet 

    linefit = ifitspec_combine_blends(linefit)
    
; measure the Lick spectral indices.  should the best-fitting
; emission-line spectrum be subtracted first?  alternatively we could
; use the best-fitting continuum to measure the indices, but I think
; this is less physical in the case of template incompleteness or
; template mismatch.  these indices will need to be corrected for
; emission contamination later

    splog, 'Measuring color and Lick indices.'
    indices = spectral_indices(restwave,backfit.continuum,$
      backfit.continuum_sigma,doplot=doplot)
;   indices = spectral_indices(restwave,flux,sqrt(var),doplot=0)

; measure emission-line equivalent widths.  the Balmer lines have to
; be treated separately because they appear both in absorption and in
; emission.  we use MEASURE_LINEEW() to measure the EW of forbidden
; lines and blended Balmer lines (such as H-alpha).  then we measure
; the Balmer emission-line EW by taking the measurement of the
; continuum from the absorption-line fit and divide it by the line
; flux derived in ILINEFIT().  first flag lines that were not measured
; (no points to fit = -2.0).  we will derive upper limits for lines
; with LINEAREA_ERR = -1.0

    measured = where(linefit.linearea_err ne -2.0,nmeasured,comp=unmeasured,ncomp=nunmeasured)
    
    if (nunmeasured ne 0L) then begin

       splog, 'The following lines were not measured: ['+strjoin(linefit[unmeasured].linename,', ')+'].'

    endif

    if (nmeasured ne 0L) then begin

       balmerindx = balmerindx(linefit[measured].linewave,$
         linefit[measured].line_blend,forbidden=forbidden)

; ---------------------------------------------------------------------------       
; forbidden lines and blended Balmer lines
; ---------------------------------------------------------------------------       

       if forbidden[0] eq -1L then begin

          splog, 'No forbidden emission line(s) or blended Balmer emission line(s) were successfully fitted.' 

       endif else begin

          splog, 'Computing EWs for '+strn(n_elements(forbidden))+$
            ' forbidden emission line(s) and blended Balmer emission line(s).'

          linefit[measured[forbidden]] = measure_lineew(linefit[measured[forbidden]],$
            wave,flux,sqrt(var),lineres=elineres[measured[forbidden]],speclineewfit=$
            speclineewfit,snrcut=snrcut,doplot=doplot)

       endelse 

; ---------------------------------------------------------------------------       
; unblended Balmer lines
; ---------------------------------------------------------------------------       

       if balmerindx[0] eq -1L then begin

          splog, 'No unblended Balmer emission line(s) were successfully fitted.' 

       endif else begin

; these individual Balmer lines were either dropped from the
; multi-Gaussian fit or have low S/N
          
          ulimit = where((linefit[measured[balmerindx]].linearea_err eq -1.0) or $
            (linefit[measured[balmerindx]].linearea/linefit[measured[balmerindx]].linearea_err lt $
            snrcut),nulimit,comp=goodfit,ncomp=ngoodfit)

          if nulimit ne 0L then begin

; compute upper limits on unblended Balmer line fluxes and EW's
             
             splog, 'Computing EWs for '+strn(nulimit)+' unblended Balmer emission line(s).'
             linefit[measured[balmerindx[ulimit]]] = measure_lineew(linefit[measured[balmerindx[ulimit]]],$
               wave,flux,sqrt(var),lineres=elineres[measured[balmerindx[ulimit]]],snrcut=snrcut,$
               speclineewfit=speclineewfit,doplot=doplot)

          endif

; the remaining Balmer lines have good line fluxes measured but we
; can't use MEASURE_LINEEW() because of the contaminating absorption.
; in order to compute an EW we need a measure of the continuum at the
; line center.  if the corresponding absorption line was measured
; (e.g., H-beta) then use the continuum measurement from that line to
; compute the emission-line EW

          if (n_elements(absfit) eq 0L) and (ngoodfit ne 0L) then begin

             splog, 'Unable to compute EWs for the following unblended Balmer emission line(s): ['+$
               strjoin(linefit[measured[balmerindx[goodfit]]].linename)+'] (no absorption line information).'

             linefit[measured[balmerindx[goodfit]]].lineew_area = 0.0
             linefit[measured[balmerindx[goodfit]]].lineew_area_err = -2.0
             linefit[measured[balmerindx[goodfit]]].lineew_box = 0.0
             linefit[measured[balmerindx[goodfit]]].lineew_box_err = -2.0

          endif else begin
          
             if (ngoodfit ne 0L) then begin

                amatch = match_string(linefit[measured[balmerindx[goodfit]]].linename,$
                  absfit.linename,index=aindx,/exact)
                good = where(aindx ne -1L,ngood)
                
                if ngood ne 0L then begin

                   splog, 'Computing EWs for ['+$
                     strjoin(linefit[measured[balmerindx[goodfit[good]]]].linename,', ')+$
                     '] using the absorption line continuum.'

                   linefit[measured[balmerindx[goodfit[good]]]].linecontlevel = absfit[aindx[good]].linecontlevel
                   linefit[measured[balmerindx[goodfit[good]]]].linecontlevel_err = absfit[aindx[good]].linecontlevel_err

                   linefit[measured[balmerindx[goodfit[good]]]] = fill_ew(linefit[measured[balmerindx[goodfit[good]]]])

                   if keyword_set(doplot) then begin

                      for k = 0L, ngood-1L do begin

                         linename = linefit[measured[balmerindx[goodfit[good[k]]]]].linename
                         meanz = linefit[measured[balmerindx[goodfit[good[k]]]]].linez
                         meanwave = linefit[measured[balmerindx[goodfit[good[k]]]]].linewave*(1+meanz)
                         lineres = elineres[measured[balmerindx[goodfit[good[k]]]]]
                         
                         linecont = linefit[measured[balmerindx[goodfit[good[k]]]]].linecontlevel

                         leftbox  = meanwave - 13.0*meanwave*lineres/light ; [Angstrom]
                         rightbox = meanwave + 13.0*meanwave*lineres/light

                         scale = 1E15
                         djs_plot, wave, scale*flux, xsty=11, ysty=3, xrange=[leftbox,rightbox]*[0.99,1.01], $
                           ps=10, charsize=2.0, charthick=2.0, xthick=2.0, ythick=2.0, ymargin=[4,3], $
                           xtitle='Wavelength '+angstrom()+')', ytitle='f_{\lambda} (10^{15} '+flam_units()+')'
                         axis, /xaxis, xrange=[leftbox,rightbox]*[0.99,1.01]/(1+meanz), xthick=2.0, $
                           charsize=2.0, charthick=2.0, xtitle='Rest Wavelength ('+angstrom()+')', xsty=3
                         djs_oplot, wave, scale*(backfit.continuum+speclinefit), color='red', thick=5.0, ps=10
                         djs_oplot, [leftbox,rightbox], scale*linecont*[1,1], line=0, thick=5.0, color='navy'
                         legend, linename, /left, /top, box=0, charsize=2.0, charthick=2.0

                         cc = get_kbrd(1)
                         
                      endfor
                         
                   endif
                   
                endif else begin

                   splog, 'Unable to compute EWs for the following unblended Balmer emission line(s): ['+$
                     strjoin(linefit[measured[balmerindx[goodfit]]])+'] (no absorption line information).'

                   linefit[measured[balmerindx[goodfit]]].lineew_area = 0.0
                   linefit[measured[balmerindx[goodfit]]].lineew_area_err = -2.0
                   linefit[measured[balmerindx[goodfit]]].lineew_box = 0.0
                   linefit[measured[balmerindx[goodfit]]].lineew_box_err = -2.0

                endelse
                
             endif 

          endelse 

       endelse 
          
    endif 
    
;   djs_plot, wave, flux, xsty=3, ysty=3, ps=10
;   djs_oplot, wave, backfit.continuum, ps=10, color='green', thick=3.0
;   djs_oplot, wave, backfit.continuum+speclinefit, ps=10, color='red'

return
end
