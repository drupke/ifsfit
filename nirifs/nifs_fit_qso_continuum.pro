;
; History
;  13mar05  DSNR  created
;


function nifs_fit_qso_continuum,lambda,flux,weight,template_flux,index,$
                                ct_coeff,fitord=fitord,qsotmp=qsotmp,$
                                qsoord=qsoord,z=z,expterms=expterms,$
                                no_dust=no_dust,quiet=quiet,qsoonly=qsoonly

  if ~keyword_set(fitord) then fitord=3
  if ~keyword_set(expterms) then expterms=0
  if ~keyword_set(qsoonly) then qsoonly=0 else fitord=0

  ncomp = n_elements(z.gas)

  ilambda=lambda[index]
  iflux=flux[index]
  iweight=weight[index]

  qsotmpspec = readspec(qsotmp)
  qsotmpwave = qsotmpspec[*,0]
  qsotmpflux = qsotmpspec[*,1]

  iqsotmpflux = interpol(qsotmpflux,qsotmpwave,lambda)
  iiqsotmpflux = interpol(qsotmpflux,qsotmpwave,ilambda)

  npar = 2+(fitord+qsoord+2*expterms)*2
  param = dblarr(npar)

  fcn = 'nifs_qso_cnt_fcn'
  fcnargs = {fitord:fitord,$
             qsoflux:iiqsotmpflux,$
             qsoord:qsoord,$
             qsoonly:qsoonly,$
             expterms:expterms}
  parinfo = REPLICATE({$
            limited:[1b,0b],$
            limits:[0d,0d],$
            fixed:0b,$
            mpprint:0b $
                      },npar)
  
  fluxfit = mpcurvefit(ilambda,iflux,iweight,param,function_name=fcn,$
                       functargs=fcnargs,/noderiv,parinfo=parinfo,$
                       /quiet)

  nifs_qso_cnt_fcn,lambda,param,continuum,expterms=expterms,$
                   fitord=fitord,qsoflux=iqsotmpflux,qsoord=qsoord,$
                   qsoonly=qsoonly

; Re-fit around emission lines if ncomp != 0

  if z.gas[0] ne -1 then begin
  
     if ncomp le 2 then begin
        el_fitord = [2,2,2,2,2]
        ;; el_fitord = [2,1,1,1,2]
        el_ran_rest_lo = [18756d -150d,$
                          19576d -150d,$
                          20338d -150d,$
                          21218d -150d,$
                          21661d -150d]
        el_ran_rest_hi = [18756d +150d,$
                          19576d +150d,$
                          20338d +150d,$
                          21218d +150d,$
                          21661d +150d]
     endif else if ncomp eq 3 then begin
        el_fitord = [2,1,1,1,2]
        el_ran_rest_lo = [18756d -150d,$
                          19576d -200d,$
                          20338d -200d,$
                          21218d -200d,$
                          21661d -150d]
        el_ran_rest_hi = [18756d +150d,$
                          19576d +200d,$
                          20338d +200d,$
                          21218d +200d,$
                          21661d +150d] 
     endif
    
     for i=0,n_elements(el_fitord)-1 do begin
        el_ran = nifs_redshift_spec([el_ran_rest_lo[i],el_ran_rest_hi[i]],z)
        i_el = where(lambda ge el_ran[0] AND lambda le el_ran[1])
        ii_el = where(ilambda ge el_ran[0] AND ilambda le el_ran[1])
        parinfo = replicate({value:0d},el_fitord[i]+1)
        el_pars = mpfitfun('poly',ilambda[ii_el],$
                           iflux[ii_el]-fluxfit[ii_el],$
                           1d/sqrt(iweight[ii_el]),$
                           parinfo=parinfo,/quiet)
        
        continuum[i_el] += poly(lambda[i_el],el_pars)
     endfor

  endif

  ct_coeff=param
  return,continuum

end
