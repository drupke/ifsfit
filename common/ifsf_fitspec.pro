; docformat = 'rst'
;
;+
;
; This function is the core routine to fit the continuum and emission
; lines of a spectrum.
;
; The function requires an initialization structure with one required
; and a bunch of optional tags, specified in INITTAGS.txt.
;
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    A structure that contains the fit and much else ...
;
; :Params:
;    lambda: in, required, type=dblarr(npix)
;      Spectrum, observed-frame wavelengths.
;    flux: in, required, type=dblarr(npix)
;      Spectrum, fluxes.
;    err: in, required, type=dblarr(npix)
;      Spectrum, flux errors.
;    zstar: in, required, type=structure
;      Initial guess for stellar redshift
;    linelist: in, required, type=hash(lines)
;      Emission line rest frame wavelengths.
;    linelistz: in, required, type=hash(lines\,ncomp)
;      Emission line observed frame wavelengths.
;    ncomp: in, required, type=hash(lines)
;      Number of components fit to each line.
;    initdat: in, required, type=structure
;      Structure of initialization parameters, with tags specified in
;      INITTAGS.txt.
;
; :Keywords:
;    maskwidths: in, optional, type=hash(lines\,maxncomp)
;      Widths, in km/s, of regions to mask from continuum fit. If not
;      set, routine defaults to +/- 500 km/s. Can also be set in INITDAT. 
;      Routine prioritizes the keyword definition.
;    peakinit: in, optional, type=hash(lines\,maxncomp)
;      Initial guesses for peak emission-line flux densities. If not
;      set, routine guesses from spectrum. Can also be set in INITDAT.
;      Routine prioritizes the keyword definition.
;    siginit_gas: in, optional, type=hash(lines\,maxncomp)
;      Initial guess for emission line widths for fitting.
;    tweakcntfit: in, optional, type=dblarr(3\,nregions)
;      Parameters for tweaking continuum fit with localized polynomials. For 
;      each of nregions regions, array contains lower limit, upper limit, and 
;      polynomial degree.
;    quiet: in, optional, type=byte
;      Use to prevent detailed output to screen. Default is to print
;      detailed output.
; 
; :Author:
;    David S. N. Rupke::
;      Rhodes College
;      Department of Physics
;      2000 N. Parkway
;      Memphis, TN 38104
;      drupke@gmail.com
;
; :History:
;    ChangeHistory::
;      2009, DSNR, copied base code from Harus Jabran Zahid
;      2009may, DSNR, tweaked for LRIS data
;      2009jun/jul, DSNR, rewritten
;      2010jan28, DSNR, fitting now done in observed frame, not rest frame
;      2010mar18, DSNR, added ct_coeff output to continuum fit
;      2013sep, DSNR, complete re-write
;      2013nov13, DSNR, renamed, added license and copyright 
;      2013nov25, DSNR, changed structure tags of output spectra for clarity
;      2013dec09, DSNR, removed stellar z and sig optimization;
;                       added PPXF option
;      2013dec10, DSNR, removed docs of initdat tags, since it's
;                       repeated in INITTAGS.txt; removed linelabel
;                       parameter, since it's in initdat; changed
;                       'initstr' parameter to 'initdat', for
;                       consistency with IFSF; testing and bug fixes
;      2013dec11, DSNR, added MASK_HALFWIDTH variable; changed value
;                       from 500 to 1000 km/s
;      2013dec12, DSNR, added SIGINIT_GAS_DEFAULT variable
;      2013dec17, DSNR, started propagation of hashes through code and 
;                       implementation of new calling sequence rubric
;      2014jan13, DSNR, propagated use of hashes
;      2014jan16, DSNR, updated treatment of redshifts; bugfixes
;      2014jan17, DSNR, bugfixes; implemented SIGINIT_GAS, TWEAKCNTFIT keywords
;      2014feb17, DSNR, removed code that added "treated" templates
;                       prior to running a generic continuum fitting
;                       routine (rebinning, adding polynomials, etc.);
;                       i.e., generic continuum fitting routine is now
;                       completely generic
;      2014feb26, DSNR, replaced ordered hashes with hashes
;      2014apr23, DSNR, changed MAXITER from 1000 to 100 in call to MPFIT
;      2016jan06, DSNR, allow no emission line fit with initdat.noemlinfit
;      2016feb02, DSNR, handle cases with QSO+stellar PPXF continuum fits
;         
; :Copyright:
;    Copyright (C) 2013--2016 David S. N. Rupke
;
;    This program is free software: you can redistribute it and/or
;    modify it under the terms of the GNU General Public License as
;    published by the Free Software Foundation, either version 3 of
;    the License or any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;    General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program.  If not, see
;    http://www.gnu.org/licenses/.
;
;-
function ifsf_fitspec,lambda,flux,err,zstar,linelist,linelistz,$
                      ncomp,initdat,maskwidths=maskwidths,$
                      peakinit=peakinit,quiet=quiet,siginit_gas=siginit_gas,$
                      tweakcntfit=tweakcntfit

  c = 299792.458d        ; speed of light, km/s
  maskwidths_def = 1000d ; default half-width in km/s for emission line masking
  siginit_gas_def = 100d ; default sigma for initial guess 
                         ; for emission line widths
  nlines = n_elements(initdat.lines)

  if keyword_set(quiet) then quiet=1b else quiet=0b
  if tag_exist(initdat,'fcnlinefit') then fcnlinefit=initdat.fcnlinefit $
  else fcnlinefit='ifsf_manygauss'
  if tag_exist(initdat,'argslinefit') then argslinefit=initdat.argslinefit
  if tag_exist(initdat,'nomaskran') then nomaskran=initdat.nomaskran $
  else nomaskran=0b
  if tag_exist(initdat,'startempfile') then istemp = 1b else istemp=0b
  if tag_exist(initdat,'loglam') then loglam=1b else loglam=0b
  if tag_exist(initdat,'vacuum') then vacuum=1b else vacuum=0b
  if tag_exist(initdat,'dored') then redinit=1d else redinit=[]

  if istemp then begin
;    Get stellar templates
     restore,initdat.startempfile
;    Redshift stellar templates
     templatelambdaz = $
        reform(template.lambda,n_elements(template.lambda)) * (1d + zstar)
     if vacuum then airtovac,templatelambdaz
  endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Pick out regions to fit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  gd_indx = indgen(n_elements(flux))

; Find where spectrum extends blueward or redward of template
  if istemp then begin
     bt_indx = where(lambda lt min(templatelambdaz) OR $
                     lambda gt max(templatelambdaz),ctbt)
     if ctbt gt 0 then gd_indx = cmset_op(gd_indx,'AND',/NOT2,bt_indx)
  endif
; Find where spectrum is outside fit range
  if tag_exist(initdat,'fitran') then begin
     fitran = initdat.fitran
     or_indx = where(lambda lt initdat.fitran[0] OR $
                     lambda gt initdat.fitran[1],ctor)
     if ctor gt 0 then $
        gd_indx = cmset_op(gd_indx,'AND',/NOT2,or_indx)
  endif else begin
     fitran = 0
  endelse

; Limit data to "good" regions
  npix     = n_elements(gd_indx)
  gdflux   = flux[gd_indx]
  gdlambda = lambda[gd_indx]
  gderr    = err[gd_indx]

  if fitran[0] eq 0 then fitran = [gdlambda[0],gdlambda[npix-1]]

; Find where flux is negative and inverse variance is undefined;
; otherwise MPFIT chokes.
  neg_indx = where(gdflux lt 0,ct)
  if ct gt 0 then begin
     gdflux[neg_indx]=-1d*gdflux[neg_indx]
     gderr[neg_indx]=max(gderr)
     if not quiet then $
        print,"Setting ",ct," points from neg. flux to pos. and max(err).",$
              format='(A,I0,A)'
  endif
  zer_indx = where(gdflux eq 0 OR gderr le 0,ct)
  if ct gt 0 then begin
     gdflux[zer_indx]=median(gdflux)
     gderr[zer_indx]=max(gderr)
     if not quiet then $
        print,"Setting ",ct,$
              " points from zero flux or error to med(flux) and max(err).",$
              format='(A,I0,A)'
  endif
  inf_indx = where(finite(gderr,/infinity),ct)
  if ct gt 0 then begin
     gdflux[inf_indx]=median(gdflux)
     gderr[inf_indx]=max(gderr)
     if ~ quiet then $
        print,"Setting ",ct,$
              " points from inf. error to med(flux) and max(err).",$
              format='(A,I0,A)'
  endif

; Weight
  gdweight = 1d/gderr^2

; timer
  fit_time0 = systime(1)
      

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Fit continuum
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if tag_exist(initdat,'fcncontfit') then begin
     
;    Mask emission lines
     if ~ tag_exist(initdat,'noemlinfit') then begin
        if not keyword_set(maskwidths) then $
           if tag_exist(initdat,'maskwidths') then $
              maskwidths = initdat.maskwidths $
           else begin
              maskwidths = hash(initdat.lines)
              foreach line,initdat.lines do $
                 maskwidths[line] = dblarr(initdat.maxncomp)+maskwidths_def
           endelse
        ct_indx  = ifsf_masklin(gdlambda, linelistz, maskwidths, $
                                nomaskran=nomaskran)
     endif else ct_indx = indgen(n_elements(gdlambda))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Option 1: Input function
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

     if initdat.fcncontfit ne 'ppxf' then begin

        if tag_exist(initdat,'argscontfit') then continuum = $
           call_function(initdat.fcncontfit,gdlambda,gdflux,$
                         gdweight,new_temp,ct_indx,ct_coeff,$
                         quiet=quiet,_extra=initdat.argscontfit) $
        else continuum = $
           call_function(initdat.fcncontfit,gdlambda,gdflux,$
                         gdweight,new_temp,ct_indx,ct_coeff,quiet=quiet)
        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Option 2: PPXF
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
     endif else if (istemp AND $
                    tag_exist(initdat,'siginit_stars')) then begin

;       Log rebin galaxy spectrum
        log_rebin,[gdlambda[0],gdlambda[n_elements(gdlambda)-1]],gdflux,$
                  gdflux_log,gdlambda_log,velscale=velscale
        log_rebin,[gdlambda[0],gdlambda[n_elements(gdlambda)-1]],gderr^2d,$
                  gderrsq_log
        gderr_log = sqrt(gderrsq_log)
        
;       Interpolate template to same grid as data
        temp = ifsf_interptemp(gdlambda,templatelambdaz,template.flux)
        temp_log = ifsf_interptemp(gdlambda_log,alog(templatelambdaz),$
                                   template.flux)

;       Mask emission lines in log space
        ct_indx_log = $
           ifsf_masklin(exp(gdlambda_log), linelistz, maskwidths, $
                        nomaskran=nomaskran)

;       Check polynomial degree
        polyterms = 4
        if tag_exist(initdat,'ppxf_maxdeg_addpoly') then $
           polyterms = initdat.ppxf_maxdeg_addpoly

;       This ensures PPXF doesn't look for lambda if no reddening is done
        if n_elements(redinit) eq 0 then redlambda = [] else redlambda=gdlambda

; Attempt to add QSO template as sky spectrum ... didn't work.
        if tag_exist(initdat,'qsotempfile') then begin
           restore,initdat.qsotempfile
           log_rebin,[gdlambda[0],gdlambda[n_elements(gdlambda)-1]],$
                     struct.cont_fit[gd_indx],$
                     gdqsotemp_log,gdlambda_log_tmp
           sky=gdqsotemp_log       
        endif else sky=0b

        ppxf,temp_log,gdflux_log,gderr_log,velscale,$
             [0,initdat.siginit_stars],sol,$
             goodpixels=ct_indx_log,bestfit=continuum_log,moments=2,$
             degree=polyterms,polyweights=polyweights,quiet=quiet,$
             weights=ct_coeff,reddening=redinit,lambda=redlambda,sky=sky

;       Resample the best fit into linear space
        continuum = interpol(continuum_log,gdlambda_log,ALOG(gdlambda))

; IFSF_CMPCONTPPXF is not yet functional ...    
;        continuum = ifsf_cmpcontppxf(gdlambda,gdlambda_log,temp,ct_coeff,$
;                                     polyterms,polyweights)

;       Adjust stellar redshift based on fit
        zstar += sol[0]/c

     endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Option to tweak cont. fit with local polynomial fits
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
     if keyword_set(tweakcntfit) then begin
;       Arrays holding emission-line-masked data
        ct_lambda=gdlambda[ct_indx]
        ct_flux=gdflux[ct_indx]
        ct_err=gderr[ct_indx]
        ct_cont = continuum[ct_indx]
        for i=0,n_elements(tweakcntfit[0,*])-1 do begin
;          Indices into full data
           tmp_ind = where(gdlambda ge tweakcntfit[0,i] AND $
                           gdlambda le tweakcntfit[1,i],ct_ind)
;          Indices into masked data
           tmp_ctind = where(ct_lambda ge tweakcntfit[0,i] AND $
                             ct_lambda le tweakcntfit[1,i],ct_ctind)
           if ct_ind gt 0 AND ct_ctind gt 0 then begin
              parinfo = replicate({value:0d},tweakcntfit[2,i]+1)
              tmp_pars = mpfitfun('poly',ct_lambda[tmp_ctind],$
                                  ct_flux[tmp_ctind] - ct_cont[tmp_ctind],$
                                  ct_err[tmp_ctind],parinfo=parinfo,/quiet)
              continuum[tmp_ind] += poly(gdlambda[tmp_ind],tmp_pars)
           endif
        endfor
     endif
     
     if tag_exist(initdat,'dividecont') then begin
        gdflux_nocnt = gdflux / continuum - 1
        gdweight_nocnt = gdweight * continuum^2
        gderr_nocnt = gderr / continuum
        method   = 'CONTINUUM DIVIDED'
     endif else begin
        gdflux_nocnt = gdflux - continuum
        gdweight_nocnt = gdweight
        gderr_nocnt = gderr
        method   = 'CONTINUUM SUBTRACTED'
     endelse

  endif else begin
  
     gdflux_nocnt = gdflux
     gderr_nocnt = gderr
     method   = 'NO CONTINUUM FIT'
     ct_coeff = 0d
     ct_indx = 0d

  endelse
  
  fit_time1 = systime(1)
  if not quiet then print,'IFSF_FITSPEC: Continuum fit took ',$
                          fit_time1-fit_time0,$
                          ' s.',format='(A,D0.1,A)'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Fit emission lines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if ~ tag_exist(initdat,'noemlinfit') then begin

; Initial guesses for emission line peak fluxes (above continuum)
; If initial guess is negative, set to 0 to prevent MPFITFUN from choking 
;   (since we limit peak to be >= 0).
  if not keyword_set(peakinit) then $
     if tag_exist(initdat,'peakinit') then $
        peakinit = initdat.peakinit $
     else begin
        peakinit = hash(initdat.lines)
        foreach line,initdat.lines do begin
           peakinit[line] = interpol(gdflux_nocnt, gdlambda, linelistz[line])
           neg = where(peakinit[line] lt 0, ct)
           if ct gt 0 then peakinit[line,neg] = 0
        endforeach
     endelse
; Initial guesses for emission line widths
  if not keyword_set(siginit_gas) then $
     if not tag_exist(initdat,'siginit_gas') then begin
        siginit_gas = hash(initdat.lines)
        foreach line,initdat.lines do $
           siginit_gas[line] = dblarr(initdat.maxncomp)+siginit_gas_def
     endif else siginit_gas = initdat.siginit_gas

; Fill out parameter structure with initial guesses and constraints
  if tag_exist(initdat,'argsinitpar') then parinit = $
     call_function(initdat.fcninitpar,linelist,linelistz,$
                   initdat.linetie,peakinit,siginit_gas,initdat.maxncomp,ncomp,$
                   _extra=initdat.argsinitpar) $
  else parinit = $
     call_function(initdat.fcninitpar,linelist,linelistz,$
                   initdat.linetie,peakinit,siginit_gas,initdat.maxncomp,ncomp)

  testsize = size(parinit)
  if testsize[0] eq 0 then begin
     print,'IFSF_FITSPEC: Bad initial parameter guesses. Aborting.'
     outstr = 0
     goto,finish
  endif

  specfit = DBLARR(npix)
  param = Mpfitfun(fcnlinefit,gdlambda,gdflux_nocnt,gderr_nocnt,$
                   parinfo=parinit,perror=perror,maxiter=100,$
                   bestnorm=chisq,covar=covar,yfit=specfit,dof=dof,$
                   nfev=nfev,niter=niter,status=status,quiet=quiet,$
                   npegged=npegged,ftol=1D-6,functargs=argslinefit,$
                   errmsg=errmsg)
  if status eq 0 OR status eq -16 then begin
     print,'IFSF_FITSPEC: Error in MPFIT. Aborting.'
     outstr = 0
     goto,finish
  endif

  cont_dat = gdflux - specfit

  endif else begin
    cont_dat = gdflux
    specfit = 0
    chisq = 0
    dof = 1
    niter = 0
    status = 0
    linelist = 0
    parinit = 0
    param = 0
    perror = 0
    covar = 0
  endelse

  ; This sets the output reddening to a numerical 0 instead of NULL
  if n_elements(redinit) eq 0 then redinit=0d
  
  fit_time2 = systime(1)
  if not quiet then print,'IFSF_FITSPEC: Line fit took ',$
                          fit_time2-fit_time1,' s.',$
                          format='(A,D0.1,A)'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output structure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  outstr = {$
           fitran: fitran, $
;          Continuum fit parameters
           ct_method: method, $
           ct_coeff: ct_coeff, $
           ct_ebv: redinit, $
           zstar: zstar, $
;          Spectrum in various forms
           wave: gdlambda, $
           spec: gdflux, $      ; data
           spec_err: gderr, $
           cont_dat: cont_dat, $ ; cont. data (all data - em. line fit)
           cont_fit: continuum, $     ; cont. fit
           emlin_dat: gdflux_nocnt, $ ; em. line data (all data - cont. fit)
           emlin_fit: specfit, $      ; em. line fit
           ct_indx: ct_indx, $        ; where emission is not masked
           gd_indx: gd_indx, $        ; cuts on various criteria
;          Line fit parameters
           redchisq: chisq/dof, $
           niter: niter, $
           fitstatus: status, $
           linelist: linelist, $
           linelabel: initdat.lines, $
           parinfo: parinit, $
           param: param, $
           perror: perror, $
           covar: covar $
           }

finish:

  return, outstr
  
end
