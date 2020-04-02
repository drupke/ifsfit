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
;    siglim_gas: in, optional, type=dblarr(2)
;      Sigma limits for line fitting.
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
;      2016feb12, DSNR, changed treatment of sigma limits for emission lines
;                       so that they can be specified on a pixel-by-pixel basis
;      2016aug31, DSNR, added option to mask continuum range(s) by hand with
;                       INITDAT tag MASKCTRAN
;      2016sep13, DSNR, added internal logic to check if emission-line fit present
;      2016sep16, DSNR, allowed MASKWIDTHS_DEF to come in through INITDAT
;      2016sep22, DSNR, tweaked continuum function call to allow new continuum
;                       fitting capabilities; moved logging of things earlier
;                       instead of ensconcing in PPXF loop, for use of PPXF 
;                       elsewhere; new output tag CONT_FIT_PRETWEAK
;      2016oct03, DSNR, multiply PERROR by reduced chi-squared, per prescription
;                       in MPFIT documentation
;      2016oct11, DSNR, added calculation of fit residual
;      2016nov17, DSNR, changed FTOL in MPFITFUN call from 1d-6 to 
;                       default (1d-10)
;      2018mar05, DSNR, added option to convolve template with spectral resolution
;                       profile
;      2018may30, DSNR, added option to adjust XTOL and FTOL for line fitting
;      2018jun25, DSNR, added NOEMLINMASK switch, distinct from NOEMLINFIT
;         
; :Copyright:
;    Copyright (C) 2013--2018 David S. N. Rupke
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
function ifsf_fitspec,lambda,flux,err,dq,zstar,linelist,linelistz,$
                      ncomp,initdat,maskwidths=maskwidths,$
                      peakinit=peakinit,quiet=quiet,siginit_gas=siginit_gas,$
                      siglim_gas=siglim_gas,tweakcntfit=tweakcntfit,$
                      col=col,row=row

  flux_out = flux
  err_out = err

  c = 299792.458d        ; speed of light, km/s
  siginit_gas_def = 100d ; default sigma for initial guess 
                         ; for emission line widths
  if tag_exist(initdat,'lines') then begin
     nlines = n_elements(initdat.lines)
     linelabel = initdat.lines
  endif else begin
     linelabel = 0b
  endelse

  if keyword_set(quiet) then quiet=1b else quiet=0b
  if keyword_set(siglim_gas) then siglim_gas=siglim_gas else siglim_gas=0b
  if tag_exist(initdat,'fcnlinefit') then fcnlinefit=initdat.fcnlinefit $
  else fcnlinefit='ifsf_manygauss'
  if tag_exist(initdat,'argslinefit') then argslinefit=initdat.argslinefit
  if tag_exist(initdat,'nomaskran') then nomaskran=initdat.nomaskran $
  else nomaskran=0b
  if tag_exist(initdat,'startempfile') then istemp = 1b else istemp=0b
  if tag_exist(initdat,'loglam') then loglam=1b else loglam=0b
  if tag_exist(initdat,'vacuum') then vacuum=1b else vacuum=0b
  if tag_exist(initdat,'ebv_star') then ebv_star=initdat.ebv_star else ebv_star=[]
  if tag_exist(initdat,'maskwidths_def') then $
     maskwidths_def = initdat.maskwidths_def $
  else maskwidths_def = 1000d ; default half-width in km/s for emission line masking
  if tag_exist(initdat,'mpfit_xtol') then mpfit_xtol=initdat.mpfit_xtol $
  else mpfit_xtol=1d-10
  if tag_exist(initdat,'mpfit_ftol') then mpfit_ftol=initdat.mpfit_ftol $
  else mpfit_ftol=1d-10

  noemlinfit = 0b
  if tag_exist(initdat,'noemlinfit') then ct_comp_emlist = 0 $
  else nocomp_emlist = ncomp.where(0,complement=comp_emlist,ncomp=ct_comp_emlist)
  if ct_comp_emlist eq 0 then noemlinfit=1b

  noemlinmask = 0b
  if noemlinfit AND ~ tag_exist(initdat,'doemlinmask') then noemlinmask = 1b

  if istemp then begin
;    Get stellar templates
     restore,initdat.startempfile
;    Redshift stellar templates
     templatelambdaz = reform(template.lambda,n_elements(template.lambda))
     if ~ tag_exist(initdat,'keepstarz') then $
        templatelambdaz *= 1d + zstar
     if vacuum then airtovac,templatelambdaz
     if tag_exist(initdat,'waveunit') then $
        templatelambdaz *= initdat.waveunit
     if tag_exist(initdat,'fcnconvtemp') then begin
        if tag_exist(initdat,'argsconvtemp') then $
           newtemplate = $
              call_function(initdat.fcnconvtemp,templatelambdaz,$
                            template,_extra=initdat.argsconvtemp) $
        else $
           newtemplate = $
              call_function(initdat.fcnconvtemp,templatelambdaz,$
                            template)
     endif
  endif else begin
     templatelambdaz = lambda
  endelse

; Set up error in zstar
  zstar_err = 0d

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Pick out regions to fit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  flux_raw = flux
  err_raw = err

  if tag_exist(initdat,'fitran') then fitran_tmp = initdat.fitran $
  else fitran_tmp = [lambda[0],lambda[n_elements(lambda)-1]]
; indices locating good data and data within fit range
  gd_indx_full = where(flux ne 0 AND err gt 0 $
                       AND ~ finite(flux,/nan) $
                       AND ~ finite(flux,/infinity) $
                       AND ~ finite(err,/nan) $
                       AND ~ finite(err,/infinity) AND $
                       dq eq 0 AND $
                       lambda ge min(templatelambdaz) AND $
                       lambda le max(templatelambdaz) AND $
                       lambda ge fitran_tmp[0] AND $
                       lambda le fitran_tmp[1],ctgd_full)

  fitran = [min(lambda[gd_indx_full]),max(lambda[gd_indx_full])]

; Find where flux is <= 0 or error is <= 0 or infinite or NaN
; (Otherwise MPFIT chokes.)
  neg_indx = where(flux lt 0,ctneg)
  zerinf_indx = where(flux eq 0 OR err le 0 $
                      OR finite(flux,/infinity) $
                      OR finite(flux,/nan) $
                      OR finite(err,/infinity) $
                      OR finite(err,/nan),ctzerinf)
  maxerr = max(err[gd_indx_full])
;  if ctneg gt 0 then begin
;     flux[neg_indx]=-1d*flux[neg_indx]
;     err[neg_indx]=maxerr*100d
;     if not quiet then print,'Setting ',ctneg,' points from neg. flux to pos. '+$
;        'and max(err)x100.',format='(A,I0,A)'
;  endif
  if ctzerinf gt 0 then begin
     flux[zerinf_indx]=median(flux[gd_indx_full])
     err[zerinf_indx]=maxerr*100d
     if not quiet then print,'Setting ',ctzerinf,' points from zero/inf./NaN flux or '+$
        'neg./zero/inf./NaN error to med(flux) and max(err)x100.',format='(A,I0,A)'
  endif

; indices locating data within actual fit range
  fitran_indx = where(lambda ge fitran[0] AND lambda le fitran[1],ctfitran)

; indices locating good regions within lambda[fitran_indx]
  gd_indx_full_rezero = gd_indx_full - fitran_indx[0]
  max_gd_indx_full_rezero = max(fitran_indx) - fitran_indx[0]
  i_gd_indx_full_rezero = where(gd_indx_full_rezero ge 0 AND $
                                gd_indx_full_rezero le max_gd_indx_full_rezero,$
                                ctgd)
  gd_indx = gd_indx_full_rezero[i_gd_indx_full_rezero]

; Limit data to fitrange
  npix     = n_elements(fitran_indx)
  gdflux   = flux[fitran_indx]
  gdlambda = lambda[fitran_indx]
  gderr    = err[fitran_indx]

; Weight
  gdweight = 1d/gderr^2

; Log rebin galaxy spectrum for finding bad regions in log space
  log_rebin,fitran,flux_raw[fitran_indx],gdflux_log
  log_rebin,fitran,err_raw[fitran_indx]^2d,gderrsq_log
  gderr_log = sqrt(gderrsq_log)
;  neg_indx_log = where(gdflux_log lt 0,ctneg_log)
  zerinf_indx_log = where(gdflux_log eq 0 OR gderr_log le 0 OR $
                          finite(gderr_log,/infinity),ctzerinf_log)
  gd_indx_log = indgen(ctfitran)
;  if ctneg_log gt 0 then $
;     gd_indx_log = cgsetdifference(gd_indx_log,neg_indx_log)
  if ctzerinf_log gt 0 then $
     gd_indx_log = cgsetdifference(gd_indx_log,zerinf_indx_log)
  

; Log rebin galaxy spectrum for use with PPXF, this time with 
; errors corrected before rebinning
  log_rebin,fitran,gdflux,gdflux_log,gdlambda_log,velscale=velscale
  log_rebin,fitran,gderr^2d,gderrsq_log
  gderr_log = sqrt(gderrsq_log)

; timer
  fit_time0 = systime(1)
      

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Fit continuum
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if tag_exist(initdat,'fcncontfit') then begin

;    Mask emission lines
     if ~ noemlinmask then begin
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
;       Mask emission lines in log space
        ct_indx_log = ifsf_masklin(exp(gdlambda_log), linelistz, $
                                   maskwidths, nomaskran=nomaskran)
     endif else begin
        ct_indx = indgen(n_elements(gdlambda))
        ct_indx_log = indgen(n_elements(gdlambda_log))
     endelse

     ct_indx = cgsetintersection(ct_indx,gd_indx)
     ct_indx_log = cgsetintersection(ct_indx_log,gd_indx_log)
           
;;    Mask other regions
;;    Now doing this in IFSF_FITLOOP with CUTRANGE tag
;     if tag_exist(initdat,'maskctran') then begin
;        mrsize = size(initdat.maskctran)
;        nreg = 1
;        if mrsize[0] gt 1 then nreg = mrsize[2]
;        for k=0,nreg-1 do begin
;            indx_mask = where(gdlambda ge initdat.maskctran[0,k] AND $
;                              gdlambda le initdat.maskctran[1,k],ct)
;            indx_mask_log = where(exp(gdlambda) ge initdat.maskctran[0,k] AND $
;                                  exp(gdlambda) le initdat.maskctran[1,k],ct)
;            if ct gt 0 then begin
;               ct_indx = cgsetdifference(ct_indx,indx_mask)
;               ct_indx_log = cgsetdifference(ct_indx_log,indx_mask_log)
;            endif
;        endfor
;     endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Option 1: Input function
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

     if initdat.fcncontfit ne 'ppxf' then begin

        if istemp then begin
           templatelambdaz_tmp = templatelambdaz
           templateflux_tmp = template.flux
        endif else begin
           templatelambdaz_tmp = 0b
           templateflux_tmp = 0b
        endelse

        if tag_exist(initdat,'argscontfit') then begin
           argscontfit_use = initdat.argscontfit
           if initdat.fcncontfit eq 'ifsf_fitqsohost' then $
              argscontfit_use = create_struct(argscontfit_use,'fitran',fitran)
           if tag_exist(initdat.argscontfit,'uselog') then $
              argscontfit_use = $
                 create_struct(argscontfit_use,'index_log',ct_indx_log)
           if tag_exist(initdat.argscontfit,'usecolrow') AND $
              keyword_set(col) AND keyword_set(row) then $  
              argscontfit_use = $
                 create_struct(argscontfit_use,'colrow',[col,row])
           continuum = $
              call_function(initdat.fcncontfit,gdlambda,gdflux,$
                            gdweight,templatelambdaz_tmp,templateflux_tmp,$
                            ct_indx,ct_coeff,zstar,$
                            quiet=quiet,_extra=argscontfit_use)
           ppxf_sigma=0d
           if initdat.fcncontfit eq 'ifsf_fitqsohost' AND $
              tag_exist(initdat.argscontfit,'refit') then $ 
              ppxf_sigma=ct_coeff.ppxf_sigma
        endif else begin
           continuum = $
              call_function(initdat.fcncontfit,gdlambda,gdflux,$
                            gdweight,templatelambdaz_tmp,templateflux_tmp,$
                            ct_indx,ct_coeff,zstar,quiet=quiet)
           ppxf_sigma=0d
        endelse
        add_poly_weights=0d
        ct_rchisq=0d
        ppxf_sigma_err=0d

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Option 2: PPXF
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
     endif else if (istemp AND $
                    tag_exist(initdat,'siginit_stars')) then begin
        
;       Interpolate template to same grid as data
        temp_log = ifsf_interptemp(gdlambda_log,alog(templatelambdaz),$
                                   template.flux)

;       Check polynomial degree
        add_poly_degree = 4
        if tag_exist(initdat,'argscontfit') then $
           if tag_exist(initdat.argscontfit,'add_poly_degree') then $
              add_poly_degree = initdat.argscontfit.add_poly_degree

;       This ensures PPXF doesn't look for lambda if no reddening is done
        if n_elements(ebv_star) eq 0 then redlambda = [] else redlambda=exp(gdlambda_log)

;       Add QSO template as sky spectrum so that it doesn't get convolved with
;       anything.
        if tag_exist(initdat,'qsotempfile') then begin
           restore,initdat.qsotempfile
           log_rebin,[gdlambda[0],gdlambda[n_elements(gdlambda)-1]],$
                     struct.cont_fit,$
                     gdqsotemp_log,gdlambda_log_tmp
           sky=gdqsotemp_log       
        endif else sky=0b

        ppxf,temp_log,gdflux_log,gderr_log,velscale,$
             [0,initdat.siginit_stars],sol,$
             goodpixels=ct_indx_log,bestfit=continuum_log,moments=2,$
             degree=add_poly_degree,polyweights=add_poly_weights,quiet=quiet,$
             weights=ct_coeff,reddening=ebv_star,lambda=redlambda,sky=sky,$
             error=solerr

;       Resample the best fit into linear space
        continuum = interpol(continuum_log,gdlambda_log,ALOG(gdlambda))

;       Adjust stellar redshift based on fit
        zstar += sol[0]/c
        ppxf_sigma=sol[1]

;     From PPXF docs:
;     - These errors are meaningless unless Chi^2/DOF~1 (see parameter SOL below).
;       However if one *assume* that the fit is good, a corrected estimate of the
;       errors is: errorCorr = error*sqrt(chi^2/DOF) = error*sqrt(sol[6]).
        ct_rchisq = sol[6]
        solerr *= sqrt(sol[6])
        zstar_err = sqrt(zstar_err^2d + (solerr[0]/c)^2d)
        ppxf_sigma_err=solerr[1]

     endif else begin
        add_poly_weights=0d
        ct_rchisq=0d
        ppxf_sigma=0d
        ppxf_sigma_err=0d
     endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Option to tweak cont. fit with local polynomial fits
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
     if keyword_set(tweakcntfit) then begin
        continuum_pretweak=continuum
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
     endif else continuum_pretweak=continuum
     
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

  if ~ noemlinfit then begin

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
  if not keyword_set(siginit_gas) then begin
     siginit_gas = hash(initdat.lines)
     foreach line,initdat.lines do $
        siginit_gas[line] = dblarr(initdat.maxncomp)+siginit_gas_def
  endif

;; Normalize data so it's near 1. Use 95th percentile of flux. If it's far from
;; 1, results are different, and probably less correct b/c of issues of numerical
;; precision in calculating line properties.
;  ifsort = sort(gdflux_nocnt)
;  fsort = gdflux_nocnt[ifsort]
;  i95 = fix(n_elements(gdflux_nocnt)*0.95d)
;  fnorm = fsort[i95]
;  gdflux_nocnt /= fnorm
;  gderr_nocnt /= fnorm
;  foreach line,initdat.lines do peakinit[line] /= fnorm

; Fill out parameter structure with initial guesses and constraints
  if tag_exist(initdat,'argsinitpar') then parinit = $
     call_function(initdat.fcninitpar,linelist,linelistz,$
                   initdat.linetie,peakinit,siginit_gas,initdat.maxncomp,$
                   ncomp,siglim=siglim_gas,$
                   _extra=initdat.argsinitpar) $
  else parinit = $
     call_function(initdat.fcninitpar,linelist,linelistz,$
                   initdat.linetie,peakinit,siginit_gas,initdat.maxncomp,$
                   ncomp,siglim=siglim_gas)

  testsize = size(parinit)
  if testsize[0] eq 0 then begin
     message,'Bad initial parameter guesses.'
   ;     outstr = 0
   ;     goto,finish
  endif

  specfit = DBLARR(npix)
  param = Mpfitfun(fcnlinefit,gdlambda,gdflux_nocnt,gderr_nocnt,$
                   parinfo=parinit,perror=perror,maxiter=1000,$
                   bestnorm=chisq,covar=covar,yfit=specfit,dof=dof,$
                   nfev=nfev,niter=niter,status=status,quiet=quiet,$
                   npegged=npegged,functargs=argslinefit,$
                   errmsg=errmsg,xtol=mpfit_xtol,ftol=mpfit_ftol)

;; Un-normalize fit.
;  specfit *= fnorm
;  gdflux_nocnt *= fnorm
;  gderr_nocnt *= fnorm
;  foreach line,linelist.keys() do begin
;     iline = where(parinit.line eq line)
;     ifluxpk = cgsetintersection(iline,where(parinit.parname eq 'flux_peak'),$
;                                 count=ctfluxpk)
;     param[ifluxpk] *= fnorm
;     perror[ifluxpk] *= fnorm
;  endforeach

  if status eq 0 OR status eq -16 then begin
     message,'MPFIT: '+errmsg
;     outstr = 0
;     goto,finish
  endif
  if status eq 5 then message,'MPFIT: Max. iterations reached.',/cont

; Errors from covariance matrix ...
  perror *=  sqrt(chisq/dof)
; ... and from fit residual.
  resid=gdflux-continuum-specfit
  perror_resid = perror
  sigrange = 20d
  foreach line,linelist.keys() do begin
     iline = where(parinit.line eq line)
     ifluxpk = cgsetintersection(iline,where(parinit.parname eq 'flux_peak'),$
                                 count=ctfluxpk)
     isigma = cgsetintersection(iline,where(parinit.parname eq 'sigma'))
     iwave = cgsetintersection(iline,where(parinit.parname eq 'wavelength'))
     for i=0,ctfluxpk-1 do begin
        waverange = sigrange*$
                    sqrt((param[isigma[i]]/c*param[iwave[i]])^2d + param[2]^2d)
        wlo = value_locate(gdlambda,param[iwave[i]]-waverange/2d)
        whi = value_locate(gdlambda,param[iwave[i]]+waverange/2d)
        if gdlambda[wlo] lt gdlambda[0] OR wlo eq -1 then wlo=0
        if gdlambda[whi] gt gdlambda[n_elements(gdlambda)-1] OR whi eq -1 then $
           whi=n_elements(gdlambda)-1
        if param[ifluxpk[i]] gt 0 then $
           perror_resid[ifluxpk[i]] = sqrt(mean(resid[wlo:whi]^2d))
      endfor
  endforeach
   
  outlinelist = linelist ; this bit of logic prevents overwriting of linelist
  cont_dat = gdflux - specfit

  endif else begin
    cont_dat = gdflux
    specfit = 0
    chisq = 0
    dof = 1
    niter = 0
    status = 0
    outlinelist = 0
    parinit = 0
    param = 0
    perror = 0
    perror_resid = 0
    covar = 0
  endelse

  ; This sets the output reddening to a numerical 0 instead of NULL
  if n_elements(ebv_star) eq 0 then ebv_star=0d
  
  fit_time2 = systime(1)
  if not quiet then print,'IFSF_FITSPEC: Line fit took ',$
                          fit_time2-fit_time1,' s.',$
                          format='(A,D0.1,A)'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output structure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; restore initial values
  flux = flux_out
  err = err_out

  outstr = {$
           fitran: fitran, $
;          Continuum fit parameters
           ct_method: method, $
           ct_coeff: ct_coeff, $
           ct_ebv: ebv_star, $
           zstar: zstar, $
           zstar_err: zstar_err,$
           ct_add_poly_weights: add_poly_weights,$
           ct_ppxf_sigma: ppxf_sigma,$
           ct_ppxf_sigma_err: ppxf_sigma_err,$
           ct_rchisq: ct_rchisq,$
;          Spectrum in various forms
           wave: gdlambda, $
           spec: gdflux, $      ; data
           spec_err: gderr, $
           cont_dat: cont_dat, $ ; cont. data (all data - em. line fit)
           cont_fit: continuum, $     ; cont. fit
           cont_fit_pretweak: continuum_pretweak, $ ; cont. fit before tweaking
           emlin_dat: gdflux_nocnt, $ ; em. line data (all data - cont. fit)
           emlin_fit: specfit, $      ; em. line fit
;          gd_indx is applied, and then ct_indx
           gd_indx: gd_indx, $        ; cuts on various criteria
           fitran_indx: fitran_indx, $; cuts on various criteria
           ct_indx: ct_indx, $        ; where emission is not masked
;          Line fit parameters
           noemlinfit: noemlinfit,$   ; was emission line fit done?
           noemlinmask: noemlinmask,$ ; were emission lines masked?
           redchisq: chisq/dof, $
           niter: niter, $
           fitstatus: status, $
           linelist: outlinelist, $
           linelabel: linelabel, $
           parinfo: parinit, $
           param: param, $
           perror: perror, $
           perror_resid: perror_resid, $ ; error from fit residual
           covar: covar, $
           siglim: siglim_gas $
           }

finish:

  return, outstr
  
end
