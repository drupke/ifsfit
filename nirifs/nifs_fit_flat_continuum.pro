;
; History
;  13apr30  DSNR  copies from gmos_fit_flat_continuum
;
function nifs_fit_flat_continuum,lambda,flux,weight,template_flux,index,$
                                 ct_coeff,cfitran=cfitran,fitord=fitord,$
                                 no_dust=no_dust,quiet=quiet,z=z,$
                                 fitblr=fitblr,norefit=norefit

  if ~keyword_set(fitord) then fitord=3
  if ~keyword_set(norefit) then norefit=0 else norefit=1

  ilambda=lambda[index]
  iflux=flux[index]
  iweight=weight[index]
  err = 1d/sqrt(iweight)

  if keyword_set(cfitran) then begin
     ifit = where(ilambda ge cfitran[0] AND ilambda le cfitran[1])
     fluxfit = poly_fit(ilambda[ifit],iflux[ifit],fitord,$
                        measure=err[ifit],yfit=yfit)
  endif else begin
     fluxfit = poly_fit(ilambda,iflux,fitord,measure=err,yfit=yfit)
  endelse
  continuum = poly(lambda,fluxfit)
  
; Re-fit around emission lines if ncomp != 0. Note that the fit ranges
; below are tailored specifically to F08572+3915:NW.

  if z.gas[0] ne -1 AND ~ norefit then begin

     if keyword_set(fitblr) then begin

        el_fitord = [1,1,1,1]
        el_ran_rest_lo = [19600d,$ ; around Paa and HeI 1.87; limits 
                                ; avoid telluric (H2O) feature and
                                ; funny cnt shape at edge of spectrum
                          21350d,$ ; around H2 S(2)
                          22250d,$ ; around H2 S(1)
                          22750d]  ; around Brg
        
        el_ran_rest_hi = [20400d,$
                          21750d,$
                          22700d,$
                          23150d]
        
     endif else begin

        el_fitord = [2,1,1,1,1]
        el_ran_rest_lo = [19650d,$ ; around Paa and HeI 1.87; limits 
                                   ; avoid telluric (H2O) feature and
                                   ; funny cnt shape at edge of spectrum
                          20400d,$ ; around Brd and H2 S(3)
                          21350d,$ ; around H2 S(2)
                          22250d,$ ; around H2 S(1)
                          22750d]  ; around Brg

        el_ran_rest_hi = [20000d,$
                          21000d,$
                          21750d,$
                          22700d,$
                          23150d]
    
     endelse

     for i=0,n_elements(el_fitord)-1 do begin
        el_ran = [el_ran_rest_lo[i],el_ran_rest_hi[i]]
        i_el = where(lambda ge el_ran[0] AND lambda le el_ran[1])
        ii_el = where(ilambda ge el_ran[0] AND ilambda le el_ran[1])
        if i eq 0 AND keyword_set(fitblr) then begin
           el_fit = mpfitpeak(ilambda[ii_el],$
                              iflux[ii_el]-yfit[ii_el],el_pars,$
                              error=1d/sqrt(iweight[ii_el]),nterms=3,$
                              /positive)
           continuum[i_el] += gaussian(lambda[i_el],el_pars)
           print,'Gaussian center: ',el_pars[1]
           print,'Gaussian FWHM: ',el_pars[2]*2.35/el_pars[1]*299792
        endif else begin
           parinfo = replicate({value:0d},el_fitord[i]+1)
           el_pars = mpfitfun('poly',ilambda[ii_el],$
                              iflux[ii_el]-yfit[ii_el],$
                              1d/sqrt(iweight[ii_el]),$
                              parinfo=parinfo,/quiet) 
           continuum[i_el] += poly(lambda[i_el],el_pars)
        endelse
     endfor

  endif

  ct_coeff=0d

  return,continuum

end
