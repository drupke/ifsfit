; docformat = 'rst'
;
;+
;
; This function is the core routine to fit the continuum and emission
; lines of a spectrum.
;
; The function requires an initialization structure with one required
; and a bunch of optional tags. They are as follows::
;
;      fcninitpar: in, required, type=string
;        Name of function for initializing continuum.
;      argsaddpoly2temp: in, optional, type=structure
;        Arguments for IFSF_ADDPOLY2TEMP call.
;      argscontfit: in, optional, type=structure
;        Arguments for continuum fit routine.
;      argsinitpar: in, optional, type=structure
;        Arguments for parameter initialization routine.
;      argslinefit: in, optional, type=structure
;        Arguments for line fitting routine
;      argsoptstelz: in, optional, type=structure
;        Arguments for stellar redshift optimization.
;      fcncontfit: in, optional, type=string
;        Name of continuum fitting function. If not specified, no
;        continuum is fit. If this equals 'ppxf,' then PPXF is used to
;        fit the stellar velocity, velocity dispersion, and template
;        weights. PPXF requires a template to be chosen and an initial
;        guess for sigma (siginit_stars) to be specified.
;      fcnlinefit: in, optional, type=string
;        Name of line fitting function. Default: IFSF_MANYGAUSS
;      fcnoptstelsig: in, optional, type=string
;        Name of routine to optimize stellar dispersion.
;      fcnoptstelz: in, optional, type=string
;        Name of routine to optimize stellar redshift.
;      fitran: in, optional, type=dblarr(2)
;        Range of fitting, in observed frame. If not set, default is
;        entire range of data / template intersection.
;      loglam: in, optional, type=byte
;        Set if data has constant log(lambda) dispersion.
;      maskwidths: in, optional, type=dblarr(nlines*ncomp)
;        Width, in km/s, of regions to mask from continuum fit. If not
;        set, routine defaults to +/- 500 km/s. If parameter has one
;        value, then this half-width is applied to all emission
;        lines. If it has multiple values, it should have exactly the
;        same number of elements as lines that are being fit.
;      nomaskran: in, optional, type=dblarr(2)
;        Wavelength region *not* to mask.
;      peakinit: in, optional, type=dblarr(nlines,ncomp)
;        Initial peak flux guesses.
;      siginit_gas: in, optional, type=dblarr(nlines,ncomp)
;        Initial line width guesses, in sigma and km/s.
;      siginit_stars: in, optional, type=double
;        Initial sigma value, in km/s, for a Gaussian kernel for
;        convolving with stellar template. Convolution only performed
;        if this param is set.
;      sigfitvals: in, optional, type=dblarr
;        If this param is set, routine cross-correlates data with
;        continua convolved with each sigma value in this array, and
;        chooses the sigma with the highest correlation coeff.
;      startempfile: in, optional, type=structure
;        File containing IDL save file (usually ending in .xdr) of
;        stellar templates. Tags are lambda [type=dblarr(nwave)] and
;        flux [type=dblarr(nwave,ntemplates)].
;      dividecont: in, optional, type=byte
;        Set this param to divide the data by the continuum
;        fit. Default is to subtract.
;      vacuum: in, optional, type=byte
;        Set this param to shift stellar templates from air to
;        vacuum wavelengths.
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
;    z: in, required, type=structure
;      Structure of initial guesses for redshifts. Only required tag
;      is STAR (in, type=double; redshift used to shift template to
;      observed frame).
;    linelabel: in, required, type=strarr(nlines)
;      Emission line labels.
;    linewave: in, required, type=dblarr(nlines)
;      Emission line rest frame wavelengths.
;    linewavez: in, required, type=dblarr(nlines\,ncomp)
;      Emission line observed frame wavelengths.
;    linetie: in, required, type=hash(linelabels\,tielabels)
;      Name of emission line to which each emission line is tied
;      (in redshift and linewidth).
;    ncomp: in, required, type=???arr(nlines)
;      Number of components fit to each line.
;    initstr: in, required, type=structure
;      Structure of initialization parameters.;
;
; :Keywords:
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
;    
; :Copyright:
;    Copyright (C) 2013 David S. N. Rupke
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
function ifsf_fitspec,lambda,flux,err,z,linelabel,linewave,linewavez,$
                      linetie,ncomp,initstr,quiet=quiet

  c = 299792.458d               ; speed of light, km/s

  if keyword_set(quiet) then quiet=1b else quiet=0b
  if tag_exist(initstr,'fcnlinefit') then fcnlinefit=initstr.fcnlinefit $
  else fcnlinefit='ifsf_manygauss'
  if tag_exist(initstr,'argslinefit') then argslinefit=initstr.argslinefit
  if tag_exist(initstr,'nomaskran') then nomaskran=1b else nomaskran=0b
  if tag_exist(initstr,'startempfile') then istemp = 1b else istemp=0b
  if tag_exist(initstr,'loglam') then loglam=1b else loglam=0b
  if tag_exist(initstr,'vacuum') then vacuum=1b else vacuum=0b

  if istemp then begin
;    Get stellar templates
     restore,initstr.startempfile
;    Redshift stellar templates
     templatelambdaz = template.lambda * (1d + z.star)
     if vacuum then airtovac,templatelambdaz
  endif

;  Initialize emission line list
  nlines = n_elements(linelabel)

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
  if tag_exist(initstr,'fitran') then begin
     fitran = initstr.fitran
     or_indx = where(lambda lt initstr.fitran[0] OR $
                     lambda gt initstr.fitran[1],ctor)
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

  if tag_exist(initstr,'fcncontfit') then begin
     
;    Mask emission lines
     masklines = reform(linewavez,max(ncomp)*nlines)
;    Estimated sigma for masking emission lines and initiating fit
     if not tag_exist(initstr,'maskwidths') then $
        maskwidths = replicate(500,max(ncomp)*nlines) $
     else if n_elements(initstr.maskwidths) eq 1 then $
        maskwidths = replicate(initstr.maskwidths,max(ncomp)*nlines) $
     else maskwidths = initstr.maskwidths
     ct_indx  = ifsf_masklin(exp(gdlambda), masklines, maskwidths, $
                             nomaskran=nomaskran)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Option 1: Input function
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

     if initstr.fcncontfit ne 'ppxf' then begin
        
;       Prepare templates
        if istemp then begin
           new_temp = template.flux
;          Interpolate template to same grid as data
           new_temp = ifsf_interptemp(gdlambda,templatelambdaz,new_temp)
;          If requested, convolve template with Gaussian
           if tag_exist(initstr,'siginit_stars') then begin
              new_temp_undisp = new_temp
              new_temp = ifsf_disptemp(new_temp,gdlambda,$
                                       initstr.siginit_stars,loglam=loglam)
           endif
;          Add polynomials to templates
           if tag_exist(initstr,'argsaddpoly2temp') then new_temp = $
              call_function('ifsf_addpoly2temp',new_temp,$
                            _extra=initstr.argsaddpoly2temp) $
           else new_temp = call_function('ifsf_addpoly2temp',new_temp)
        endif else begin
           new_temp = 0
        endelse

        if tag_exist(initstr,'argscontfit') then continuum = $
           call_function(initstr.fcncontfit,gdlambda,gdflux,$
                         gdweight,new_temp,ct_indx,ct_coeff,$
                         quiet=quiet,_extra=initstr.argscontfit) $
        else continuum = $
           call_function(initstr.fcncontfit,gdlambda,gdflux,$
                         gdweight,new_temp,ct_indx,ct_coeff,quiet=quiet)
        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Option 2: PPXF
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
     endif else if (istemp AND $
                    tag_exist(initstr,'siginit_stars')) then begin

;       Log rebin galaxy spectrum
        log_rebin,[gdlambda[0],gdlambda[1]],gdflux,gdflux_log,gdlambda_log,$
                  velscale=velscale
        log_rebin,[gdlambda[0],gdlambda[1]],gderr^2d,gderrsq_log
        gderr_log = sqrt(gderrsq_log)
        
;       Interpolate template to same grid as data
        temp = ifsf_interptemp(gdlambda,templatelambdaz,template.flux)
        temp_log = ifsf_interptemp(gdlambda_log,templatelambdaz,template.flux)

;       Mask emission lines in log space
        ct_indx_log  = ifsf_masklin(exp(gdlambda_log), masklines, maskwidths, $
                                    nomaskran=nomaskran)

;       Check polynomial degree
        if tag_exist(initstr.argsaddpoly2temp,'nterms') then $
           polyterms = initstr.argsaddpoly2temp.nterms $
        else polyterms = 4

        ppxf,temp_log,gdflux_log,gderr_log,velscale,[0,siginit_stars],$
             sol=sol,goodpixels=ct_indx,bestfit=continuum_log,moments=2,$
             degree=polyterms,polyweights=polyweights,quiet=quiet,$
             weights=ct_coeff
     
        continuum = ifsf_cmpcontppxf(gdlambda,gdlambda_log,temp,ct_coeff,$
                                     polyterms,polyweights)

     endif
        
     if tag_exist(initstr,'dividecont') then begin
        gdflux_nocnt = gdflux / continuum - 1
        gdweight_nocnt = gdweight * continuum^2
        gderr_nocnt = gderr / continuum
        method   = 'CONTINUUM DIVIDED'
        cont_sig = stdev(gdflux_nocnt[ct_indx])
     endif else begin
        gdflux_nocnt = gdflux - continuum
        gdweight_nocnt = gdweight
        gderr_nocnt = gderr
        method   = 'CONTINUUM SUBTRACTED'
        cont_sig = stdev(gdflux_nocnt[ct_indx])
     endelse

  endif else begin
  
     gdflux_nocnt = gdflux
     gderr_nocnt = gderr
     cont_sig = 0d
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

; Initial guesses for emission line peak fluxes (above continuum)
; If initial guess is negative, set to 0 to prevent MPFITFUN from choking 
;   (since we limit peak to be >= 0).
  if not tag_exist(initstr,'peakinit') then begin
     peakinit = dblarr(nlines,max(ncomp))
     for i=0,max(ncomp)-1 do $
        peakinit[*,i] = interpol(gdflux_nocnt, gdlambda, linewavez[*,i])
     neg = where(peakinit lt 0, ct)
     if ct gt 0 then peakinit[neg] = 0
  endif else peakinit = initstr.peakinit
; Initial guesses for emission line widths
  if not tag_exist(initstr,'siginit_gas') then $
     siginit_gas = dblarr(nlines,max(ncomp))+100d $
  else siginit_gas = initstr.siginit_gas

; Fill out parameter structure with initial guesses and constraints
  if tag_exist(initstr,'argsinitpar') then parinit = $
     call_function(initstr.fcninitpar,linelabel,linewave,linewavez,$
                   linetie,peakinit,siginit_gas,ncomp,$
                   _extra=initstr.argsinitpar) $
  else parinit = $
     call_function(initstr.fcninitpar,linelabel,linewave,linewavez,$
                   linetie,peakinit,siginit_gas,ncomp)

  testsize = size(parinit)
  if testsize[0] eq 0 then begin
     print,'IFSF_FITSPEC: Bad initial parameter guesses. Aborting.'
     outstr = 0
     goto,finish
  endif

  specfit = DBLARR(npix)
  param = Mpfitfun(fcnlinefit,gdlambda,gdflux_nocnt,gderr_nocnt,$
                   parinfo=parinit,perror=perror,maxiter=1000,$
                   bestnorm=chisq,covar=covar,yfit=specfit,dof=dof,$
                   nfev=nfev,niter=niter,status=status,quiet=quiet,$
                   npegged=npegged,ftol=1D-6,functargs=argslinefit,$
                   errmsg=errmsg)
  if status eq 0 OR status eq -16 then begin
     print,'IFSF_FITSPEC: Error in MPFIT. Aborting.'
     outstr = 0
     goto,finish
  endif
  
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
           method: method, $
           sig_stars: bestsig_stars, $
           ct_sig: cont_sig, $
           ct_coeff: ct_coeff, $
;          Spectrum in various forms
           wave: gdlambda, $
           spec: gdflux, $      ; data
           spec_err: gderr, $
           cont_dat: continuum-specfit, $ ; cont. data (all data - em. line fit)
           cont_fit: continuum, $         ; cont. fit
           emlin_dat: gdflux_nocnt, $ ; em. line data (all data - cont. fit)
           emlin_fit: specfit, $      ; em. line fit
           ct_indx: ct_indx, $        ; where emission is not masked
           gd_indx: gd_indx, $        ; cuts on various criteria
;          Line fit parameters
           redchisq: chisq/dof, $
           niter: niter, $
           fitstatus: status, $
           linewave: linewave, $
           linelabel: linelabel, $
           param: param, $
           perror: perror, $
           covar: covar $
           }

finish:

  return, outstr
  
end
