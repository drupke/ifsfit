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
; Criteria for calling a parameter as an independent parameter in the function
; call rather than calling a parameter as part of initdat structure: Parameter 
; must be absent from allowable initdat tags defined in initialization procedure 
; (i.e., no adding new structure tags in IFSF if they are not definable in the 
; initialization procedure). If they are definable in the initialization 
; procedure but have not been, they can be added in to initdat in IFSF. 
; (This ultimately would be better handled by objects or hashes, probably.) If 
; they are defined in the initialization procedure but need to be redefined 
; with a different dimensionality (e.g., we want a different # of components 
; in each spaxel, but IFSF_FITSPEC doesn't recognize spaxels), then they 
; should be separate parameters. 
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
;    linewave: in, required, type=dblarr(nlines)
;      Emission line rest frame wavelengths.
;    linewavez: in, required, type=dblarr(nlines\,ncomp)
;      Emission line observed frame wavelengths.
;    linetie: in, required, type=strarr(nlines)
;      Name of emission line to which each emission line is tied
;      (in redshift and linewidth).
;    ncomp: in, required, type=???arr(nlines)
;      Number of components fit to each line.
;    initdat: in, required, type=structure
;      Structure of initialization parameters, with tags specified in
;      INITTAGS.txt.
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
function ifsf_fitspec,lambda,flux,err,z,linewave,linewavez,$
                      ncomp,initdat,quiet=quiet

  c = 299792.458d        ; speed of light, km/s
  mask_halfwidth = 1000d ; default half-width in km/s for emission line masking
  siginit_gas_default = 100d ; default sigma for initial guess for emission line widths

  if keyword_set(quiet) then quiet=1b else quiet=0b
  if tag_exist(initdat,'fcnlinefit') then fcnlinefit=initdat.fcnlinefit $
  else fcnlinefit='ifsf_manygauss'
  if tag_exist(initdat,'argslinefit') then argslinefit=initdat.argslinefit
  if tag_exist(initdat,'nomaskran') then nomaskran=initdat.nomaskran else nomaskran=0b
  if tag_exist(initdat,'startempfile') then istemp = 1b else istemp=0b
  if tag_exist(initdat,'loglam') then loglam=1b else loglam=0b
  if tag_exist(initdat,'vacuum') then vacuum=1b else vacuum=0b

  if istemp then begin
;    Get stellar templates
     restore,initdat.startempfile
;    Redshift stellar templates
     templatelambdaz = $
        reform(template.lambda,n_elements(template.lambda)) * (1d + z.star)
     if vacuum then airtovac,templatelambdaz
  endif

;  Initialize emission line list
  nlines = n_elements(initdat.lines)

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
     masklines = reform(linewavez,initdat.maxncomp*nlines)
;    Estimated sigma for masking emission lines and initiating fit
     if not tag_exist(initdat,'maskwidths') then $
        maskwidths = replicate(mask_halfwidth,initdat.maxncomp*nlines) $
     else if n_elements(initdat.maskwidths) eq 1 then $
        maskwidths = replicate(initdat.maskwidths,initdat.maxncomp*nlines) $
     else maskwidths = initdat.maskwidths
     ct_indx  = ifsf_masklin(gdlambda, masklines, maskwidths, $
                             nomaskran=nomaskran)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Option 1: Input function
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

     if initdat.fcncontfit ne 'ppxf' then begin
        
;       Prepare templates
        if istemp then begin
           new_temp = template.flux
;          Interpolate template to same grid as data
           new_temp = ifsf_interptemp(gdlambda,templatelambdaz,new_temp)
;          If requested, convolve template with Gaussian
           if tag_exist(initdat,'siginit_stars') then begin
              new_temp_undisp = new_temp
              new_temp = ifsf_disptemp(new_temp,gdlambda,$
                                       initdat.siginit_stars,loglam=loglam)
           endif
;          Add polynomials to templates
           if tag_exist(initdat,'argsaddpoly2temp') then new_temp = $
              call_function('ifsf_addpoly2temp',new_temp,$
                            _extra=initdat.argsaddpoly2temp) $
           else new_temp = call_function('ifsf_addpoly2temp',new_temp)
        endif else begin
           new_temp = 0
        endelse

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
        log_rebin,[gdlambda[0],gdlambda[n_elements(gdlambda)-1]],gderr^2d,gderrsq_log
        gderr_log = sqrt(gderrsq_log)
        
;       Interpolate template to same grid as data
        temp = ifsf_interptemp(gdlambda,templatelambdaz,template.flux)
        temp_log = ifsf_interptemp(gdlambda_log,alog(templatelambdaz),$
                                   template.flux)

;       Mask emission lines in log space
        ct_indx_log = $
           ifsf_masklin(exp(gdlambda_log), masklines, maskwidths, $
                        nomaskran=nomaskran)

;       Check polynomial degree
        polyterms = 4
        if tag_exist(initdat,'argsaddpoly2temp') then $
           if tag_exist(initdat.argsaddpoly2temp,'nterms') then $
              polyterms = initdat.argsaddpoly2temp.nterms

        ppxf,temp_log,gdflux_log,gderr_log,velscale,$
             [0,initdat.siginit_stars],sol,$
             goodpixels=ct_indx_log,bestfit=continuum_log,moments=2,$
             degree=polydeg,polyweights=polyweights,quiet=quiet,$
             weights=ct_coeff
     
        continuum = ifsf_cmpcontppxf(gdlambda,gdlambda_log,temp,ct_coeff,$
                                     polydeg,polyweights)

     endif
        
     if tag_exist(initdat,'dividecont') then begin
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
  if not tag_exist(initdat,'peakinit') then begin
     peakinit = dblarr(nlines,initdat.maxncomp)
     for i=0,initdat.maxncomp-1 do $
        peakinit[*,i] = interpol(gdflux_nocnt, gdlambda, linewavez[*,i])
     neg = where(peakinit lt 0, ct)
     if ct gt 0 then peakinit[neg] = 0
  endif else peakinit = initdat.peakinit
; Initial guesses for emission line widths
  if not tag_exist(initdat,'siginit_gas') then $
     siginit_gas = dblarr(nlines,initdat.maxncomp)+siginit_gas_default $
  else siginit_gas = initdat.siginit_gas

; Fill out parameter structure with initial guesses and constraints
  if tag_exist(initdat,'argsinitpar') then parinit = $
     call_function(initdat.fcninitpar,initdat.lines,linewave,linewavez,$
                   linetie,peakinit,siginit_gas,initdat.maxncomp,ncomp,$
                   _extra=initdat.argsinitpar) $
  else parinit = $
     call_function(initdat.fcninitpar,initdat.lines,linewave,linewavez,$
                   linetie,peakinit,siginit_gas,initdat.maxncomp,ncomp)

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
           linelabel: initdat.lines, $
           param: param, $
           perror: perror, $
           covar: covar $
           }

finish:

  return, outstr
  
end
