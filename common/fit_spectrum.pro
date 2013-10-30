;+
; NAME:
;     fit_spectrum
;
; PURPOSE:
;     Fit optical spectra with a stellar continuum and emission line model.
;
; CALLING SEQUENCE:
;     result = fit_spectra(lambda, flux, err, [/disperse_temp, /subtract, $
;                        /no_fit_cont, /dust, obj_id = obj_id, $
;                        vdisp = vdisp, nterms = nterms, cfite$
;                        linelist = linelist, linename = linename] )
;
; INPUTS
;   startempfile - File of stellar templates, saved to variable
;                  TEMPLATE in XDR format.
;
; OPTIONAL/KEYWORDS
;   disperse_temp - convolve template spectra with gaussian of sigma
;                   equal to vdisp
;   vdisp         - the dispersion velocity, default = 150 km which is used
;                   for masking out lines for the continuum fit and also to
;                   "smear" template spectra if the disperse_temp
;   subtract      - subtract continuum instead of dividing by continuum
;   linelist      - 
;   obj_id        - object identifier
;   no_fit_cont   - set keyword to not fit continuum
;   nterms        - sets the number of polynomial terms to be added to stellar
;                   templates. 0=constant, 1=linear, 2=quadratic etc. set to
;                   -1 to add no polynomial terms. (default = 2)
;   dust          - set if you want to fit stellar continuum for dust extinction
;
; OUTPUTS:
;
; REVISION HISTORY:
;    Written  by HJZ oct-22-2008  incorporating code from J. Moustakas
;    modified by HJZ mar-19-2009  added equivalent width to structure
;    modified by HJZ mar-26-2009  removed flux_calib keyword, need to redo 
;                                 at some point
;    modified by HJZ mar-26-2009  added ivar keyword
;    modified by HJZ apr-23-2009  changed method for calculating EW, added
;                                 subtract keyword, added input for linelist,
;                                 added method tag to output structure, took
;                                 out vdisp (dispersion velocity) from output
;                                 structure, changed spec_sig in output 
;                                 structure to cont_sig
;    modified by HJZ apr-30-2009  added object id number tag to structure
;    modified by HJZ may-7-2009   added no continuum fit keyword, added
;                                 template as a common block, changed method
;                                 tag in output structure to a string
;    modified by HJZ may-8-2009   changed keyword of no_dispersion to 
;                                 disperse_temp, default is now to not do 
;                                 the continuum template smearing, added
;                                 nterms keyword for adding polynomial terms
;                                 to stellar templates, added dust keyword
;    modified by HJZ may-10-2009  added gd_indx to output structure
;    09may      DSNR  tweaked for LRIS data
;    09jun/jul  DSNR  rewritten
;    10jan28    DSNR  fitting now done in observed frame, not rest frame
;    10mar18    DSNR  added ct_coeff output to continuum fit
;-

function fit_spectrum,lambda,flux,err,startempfile,z,ncomp,$
                      disperse_temp=disperse_temp,$
                      subtract=subtract,linelist=linelist,nterms=nterms,$
                      obj_id=obj_id,no_fit_cont=no_fit_cont,dust=dust,$
                      time=time,quiet=quiet,fitran=fitran,$
                      peakguess=peakguess,sigguess=sigguess,vdisp=vdisp,$
                      maskwidths=maskwidths,masklines=masklines,$
                      fcninitpar=fcninitpar,argsinitpar=argsinitpar,$
                      fcncontfit=fcncontfit,argscontfit=argscontfit,$
                      fcnlinefit=fcnlinefit,argslinefit=argslinefit,$
                      loglam=loglam,fcnz=fcnz,argsz=argsz,$
                      argspoly=argspoly,addmwstar=addmwstar,$
                      vacuum=vacuum

c = 2.99792458d5 ; speed of light, km/s
nsig = 5d ; sigma range for masking emission lines

if   keyword_set(dust)       then no_dust = 0 else no_dust = 1
if ~ keyword_set(vdisp)      then vdisp = 150 ;km/s
if ~ keyword_set(obj_id)     then obj_id = -1
if ~ keyword_set(nterms)     then nterms = 2
if ~ keyword_set(fcncontfit) then fcncontfit = 'fit_continuum'
if ~ keyword_set(fcnlinefit) then fcnlinefit = 'manygauss'
if ~ keyword_set(quiet)      then quiet = 0 else quiet=1
if ~ keyword_set(loglam)     then loglam = 0d
if ~ keyword_set(vacuum)     then vacuum = 0
if   ncomp eq 0              then noem = 1 else noem = 0
if startempfile ne ''        then istemp = 1 else istemp = 0

if istemp then begin
; Get stellar templates
   restore,startempfile
; Redshift stellar templates
   if keyword_set(argsz) then begin 
      templatelambdaz = call_function(fcnz,template.lambda,z,_extra=argsz)
      argsz_old = argsz
   endif else begin
      templatelambdaz = call_function(fcnz,template.lambda,z)
      argsz_old = 0
   endelse
   if keyword_set(vacuum) then airtovac,templatelambdaz
endif

;  Initialize emission line list
linewave = linelist.wave
linelabel = linelist.label
nlines = n_elements(linewave)

if ~ noem then begin

;  Redshift emission line list
   linewavez = dblarr(nlines,ncomp)
   if keyword_set(argsz) then begin
      argsz = jjadd_tag(argsz,'icomp',0)
      argsz = jjadd_tag(argsz,'gas',1)
   endif else $
      argsz = {icomp:0,gas:1}
   for i=0,ncomp-1 do begin
      argsz.icomp = i
      linewavez[*,i] = call_function(fcnz,linewave,z,_extra=argsz)
   endfor

endif

gd_indx = indgen(n_elements(flux))
if istemp then begin
; Find where spectrum extends blueward or redward of template
   bt_indx = where(lambda lt min(templatelambdaz) OR $
                   lambda gt max(templatelambdaz),ctbt)
   if ctbt gt 0 then gd_indx = cmset_op(gd_indx,'AND',/NOT2,bt_indx)
endif
; Find where spectrum is outside fit range
if keyword_set(fitran) then begin
   or_indx = where(lambda lt fitran[0] OR $
                   lambda gt fitran[1],ctor)
   if ctor gt 0 then $
      gd_indx = cmset_op(gd_indx,'AND',/NOT2,or_indx)
endif

; Limit data to "good" regions
npix     = n_elements(gd_indx)
gdflux   = flux[gd_indx]
gdlambda = lambda[gd_indx]
gderr    = err[gd_indx]

; Find where flux is negative and inverse variance is undefined;
; otherwise MPFIT chokes.
neg_indx = where(gdflux lt 0,ct)
if ct gt 0 then begin
   gdflux[neg_indx]=-1d*gdflux[neg_indx]
   gderr[neg_indx]=max(gderr)
   if ~ quiet then $
      print,"Setting ",ct," points from neg. flux to pos. and max(err).",$
            format='(A,I0,A)'
endif
zer_indx = where(gdflux eq 0 OR gderr le 0,ct)
if ct gt 0 then begin
   gdflux[zer_indx]=median(gdflux)
   gderr[zer_indx]=max(gderr)
   if ~ quiet then $
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

; Mean sigma
meansig = vdisp/c * mean(gdlambda)

if istemp then begin
   new_temp = template.flux
; Interpolate template to same grid as data
   if keyword_set(addmwstar) then begin
      if keyword_set(argsz_old) then begin
         ztmp = z
         ztmp.star = addmwstar
         temp_lam_rest = call_function(fcnz,template.lambda,ztmp,$
                                       _extra=argsz_old)
      endif else $
         temp_lam_rest = call_function(fcnz,template.lambda,ztmp)
   endif else temp_lam_rest = 0
   new_temp = interpol_template(gdlambda,templatelambdaz,new_temp,$
                                temp_lam_rest=temp_lam_rest)
; If requested, convolve template with Gaussian
   if keyword_set(disperse_temp) then begin
      dispsig = meansig
      if keyword_set(loglam) then dispsig = vdisp
      new_temp = disperse_template(new_temp, gdlambda, dispsig, loglam=loglam)
   endif
; Add polynomials to templates
   if keyword_set(argspoly) $
   then new_temp = call_function('add_poly2template',new_temp,nterms,$
                                 _extra=argspoly) $
   else new_temp = call_function('add_poly2template',new_temp,nterms)
endif else begin
   new_temp = 0
endelse

; Mask emission lines
if ~ noem then begin
   if ~ keyword_set(masklines) then $
      masklines = reform(linewavez,ncomp*nlines)
   if ~ keyword_set(maskwidths) then $
      maskwidths = nsig*replicate(meansig,nlines*ncomp)
   ct_indx  = mask_emission(gdlambda, masklines, maskwidths)
endif else $
   ct_indx = indgen(n_elements(gd_indx))

; Fit stellar continuum
if not keyword_set(no_fit_cont) then begin

   if keyword_set(argscontfit) then $
      continuum = call_function(fcncontfit,gdlambda,gdflux,gdweight,$
                                new_temp,ct_indx,ct_coeff,$
                                no_dust=no_dust,quiet=quiet,$
                                _extra=argscontfit) $
   else continuum = call_function(fcncontfit,gdlambda,gdflux,gdweight,$
                                  new_temp,ct_indx,ct_coeff,$
                                  no_dust=no_dust,quiet=quiet)

   if keyword_set(subtract) then begin
      gdflux_nocnt = gdflux - continuum
      gdweight_nocnt = gdweight
      gderr_nocnt = gderr
      method   = 'CONTINUUM SUBTRACTED'
      cont_sig = stdev(gdflux_nocnt[ct_indx])
   endif else begin
      gdflux_nocnt = gdflux / continuum - 1
      gdweight_nocnt = gdweight * continuum^2
      gderr_nocnt = gderr / continuum
      method   = 'CONTINUUM DIVIDED'
      cont_sig = stdev(gdflux_nocnt[ct_indx])
   endelse
   
endif else begin
  
   gdflux_nocnt = gdflux
   gderr_nocnt = gderr
   cont_sig = 0d
   method   = 'NO CONTINUUM FIT'
   ct_coeff = 0d

endelse


fit_time1 = systime(1)
if keyword_set(time) then print,'    Continuum fit took ',fit_time1-fit_time0,$
                                ' s.',format='(A,D0.1,A)'

if ~ noem then begin

; Initial guesses for Gaussian peak
; If initial guess is negative, set to 0 to prevent MPFITFUN from choking 
;   (since we limit peak to be >= 0).
   if ~ keyword_set(peakguess) then begin
      peakguess = dblarr(nlines,ncomp)
      for i=0,ncomp-1 do $
         peakguess[*,i] = interpol(gdflux_nocnt, gdlambda, linewavez[*,i])
      neg = where(peakguess lt 0, ct)
      if (ct gt 0) then peakguess[neg] = 0
   endif
; Initial guesses for Gaussian widths.  We assume here that if
; loglinear dispersion is set, that sigmas are expressed in velocity
; space.
   if ~ keyword_set(sigguess) then begin
      if loglam then sigguess = dblarr(nlines,ncomp)+vdisp $
      else sigguess = vdisp*linewavez/c
   endif
;
; Fit emission lines
;
; Fill out parameter structure with initial guesses and constraints
   if keyword_set(argsinitpar) then begin
      parinit = $
         call_function(fcninitpar,linelabel,linewave,linewavez,$
                       peakguess,sigguess,z,_extra=argsinitpar)
   endif else parinit = $
      call_function(fcninitpar,linelabel,linewave,$
                    linewavez,peakguess,sigguess,z)
   specfit = DBLARR(npix)
; Check parinit initial values vs. limits
   badpar = where((parinit.limited[0] AND $
                   parinit.value lt parinit.limits[0]) OR $
                  (parinit.limited[1] AND $
                   parinit.value gt parinit.limits[1]),ct)
   if ct gt 0 then begin
      print,'FIT_SPECTRUM: Init. Values are outside limits for PARINIT.'
      print,'FIT_SPECTRUM: Offending parameters:'
      print,'Index','Value','Lower limit','Upper limit',$
            format='(A7,3A15)'
      for i=0,ct-1 do begin
         j = badpar[i]
         print,j,parinit[j].value,parinit[j].limits[0],$
               parinit[j].limits[1],format='(I7,3E15.6)'
      endfor
   endif

   param = Mpfitfun(fcnlinefit,gdlambda,gdflux_nocnt,gderr_nocnt,$
                    parinfo=parinit,perror=perror,maxiter=1000,$
                    bestnorm=chisq,covar=covar,yfit=specfit,dof=dof,$
                    nfev=nfev,niter=niter,status=status,quiet=quiet,$
                    npegged=npegged,ftol=1D-6,functargs=argslinefit,$
                    errmsg=errmsg)
   if status eq -16 then begin
      print,'FIT_SPECTRUM: Function overflow.  Skipping this spectrum.'
      redchisq = 0d
      niter = 0d
      specfit = dblarr(n_elements(gd_indx))
      param = 0d
      perror = 0d
      covar = 0d
   endif

; Reduced chi-squared
   redchisq = chisq/dof

   fit_time2 = systime(1)
   if keyword_set(time) then $
      print,'    Line fit took ',fit_time2-fit_time1,' s.',$
            format='(A,D0.1,A)'

endif else begin

   redchisq = 0d
   niter = 0d
   specfit = dblarr(n_elements(gd_indx))
   param = 0d
   perror = 0d
   covar = 0d
   status = 0d

endelse

; Output structure
outstr = {obj_id       : obj_id, $
;         Fit parameters
          redchisq     : redchisq, $
          niter        : niter, $
          fitstatus    : status, $
;         Continuum fit parameters
          nterms       : nterms, $
          method       : method, $
          ct_sig       : cont_sig, $
          ct_coeff     : ct_coeff, $
;         Spectrum in various forms
          wave         : gdlambda, $
          spec         : gdflux, $
          spec_err     : gderr, $
          spec_nocnt   : gdflux_nocnt, $
          specfit      : specfit, $
          ct_indx      : ct_indx, $ ; where emission is not masked
          gd_indx      : gd_indx, $ ; cuts on various criteria
;         Line parameters
          linewave     : linewave, $
          linelabel    : linelabel, $
          param        : param, $
          perror       : perror, $
          covar        : covar }

return, outstr
  
end
