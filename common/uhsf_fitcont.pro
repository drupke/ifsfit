; docformat = 'rst'
;
;+
;
; Function to fit continuum to spectrum. Options to fit are:
;   1) fit continuum w/ templates, w/ poly. renormalization of both
;      continuum and data
;   2) fit continuum w/ templates, w/o poly renormalization (select
;      /nopoly)
;   3) fit continuum w/ polynomial (select /nobvls)
;   4) create template-based continuum using ct_coeff and renormalize
;      w/ poly. (select /nobvls and /ctinput)
;   5) create template-based continuum using ct_coeff and do not
;      renormalize w/ poly. (select /nobvls, /ctinput, and /nopoly)
;
; :Categories:
;    UHSPECFIT
;
; :Returns:
;    The best fit continuum spectrum (over all wavelengths, not just
;    those fit).
;
; :Params:
;    lambda: in, required, type=dblarr(N)
;    flux: in, required, type=dblarr(N)
;    weight: in, required, type=dblarr(N)
;    template_flux: in, required, type=dblarr(N\,M)
;      M is # of templates
;    index: in, required, type=intarr
;      Contains indices of continuum regions to fit
;    ct_coeff: out, required, type=dblarr(M)
;      Coefficients used to combine stellar templates. If /nobvls is
;      set, then the input value is used to compute the continuum. If
;      not, then the best fit coefficients are output to this
;      variable. Set to 0 if /nobvls is set and /nopoly is not.
;
; :Keywords:
;    ctinput: in, optional, type=byte
;      Use ct_coeff as an input to reproduce the specified
;      template-based continuum.
;    dust: in, optional, type=byte
;      Turn on fitting dust extinction of stellar cont.
;    fitord: in, optional, type=integer
;      Specifies order of polynomial renormalization (default = 3).
;    nopoly: in, optional, type=byte
;      Turns off polynomial renormalization.
;    nobvls: in, optional, type=byte
;      Turns off stellar template fitting.
;    quiet: in, optional, type=byte
;    refit: in, optional, type=structure
;      If set, contains structure with array of continuum regions to
;      re-fit [tag RAN, type=dblarr(2, X), where X is the number of
;      regions to refit] and array of polynomial orders to fit [tag
;      ORD, type=intarr(X)].
;
; :Author:
;    David Rupke
;
; :History:
;    Change History::
;      2009aug14, DSNR, created
;      2009dec11, DSNR, fixed bug in call to ibackfit: invar->invvar
;      2010mar18, DSNR, added ct_coeff output
;      2013oct, DSNR, complete re-write of software
;-
function uhsf_fitcont,lambda,flux,weight,template_flux,index,$
                      ct_coeff,ctinput=ctinput,dust=dust,$
                      fitord=fitord,nobvls=nobvls,$
                      nopoly=nopoly,quiet=quiet,refit=refit,$
                      addnorm

  if keyword_set(nopoly) AND $
     keyword_set(nobvls) AND $
     ~ keyword_set(ctinput) then begin
     print, 'UHSF_FITCONT: ERROR: Cannot select /nobvls and '+$
            '/nopoly w/o /ctinput.'
     print, 'UHSF_FITCONT: Returning no continuum.'
     ct_coeff = 0
     return,dblarr(n_elements(lambda))
  endif

  ilambda=lambda[index]
  iflux=flux[index]
  iweight=weight[index]

  if ~ keyword_set(nobvls) AND ~ keyword_set(ctinput) then begin
     if keyword_set(dust) then nodust=0 else nodust=1
     backfit = ibackfit(iflux,ilambda,template_flux[index, *], $
                        invvar=iweight, nodust=nodust, quiet=quiet)
     continuum = backfit.starcoeff##template_flux
     ct_coeff = backfit.starcoeff
  endif else if keyword_set(ctinput) then begin
     continuum = ct_coeff##template_flux
  endif

; Convert weights to errors
  ierr = 1d/sqrt(iweight)

  if ~ keyword_set(nopoly) then begin
;    Fit line to continuum and data for normalization
     if ~ keyword_set(fitord) then fitord=3
     fluxfit = poly_fit(ilambda,iflux,fitord,measure=ierr)
     if ~ keyword_set(nobvls) then begin
;       Normalize the fitted continuum to the flux.
        continuum *= poly(lambda,fluxfit)
        cntfit = poly_fit(ilambda,continuum[index],fitord,measure=ierr)
        continuum /= poly(lambda,cntfit)
     endif else begin
        continuum = poly(lambda,fluxfit)
        ct_coeff = 0
     endelse
  endif

  icontinuum = continuum[index]
  
  if keyword_set(refit) then begin
     for i = 0,n_elements(refit.ord)-1 do begin
        tmp_ind = where(lambda ge refit.ran[0,i] AND $
                        lambda le refit.ran[1,i])
        tmp_iind = where(ilambda ge refit.ran[0,i] AND $
                         ilambda le refit.ran[1,i])
        parinfo = replicate({value:0d},refit.ord[i]+1)
        tmp_pars = mpfitfun('poly',ilambda[tmp_iind],$
                            iflux[tmp_iind]-icontinuum[tmp_iind],$
                            ierr[tmp_iind],parinfo=parinfo,/quiet)
        continuum[tmp_ind] += poly(lambda[tmp_ind],tmp_pars)
     endfor
  endif

  return,continuum

end
