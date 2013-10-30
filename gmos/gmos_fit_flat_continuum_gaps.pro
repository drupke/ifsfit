;
; History
;  11jul20  DSNR  copied from gmos_fit_qso_cnt.pro
;
function gmos_fit_flat_continuum_gaps,lambda,flux,weight,template_flux,$
                                      index,ct_coeff,cfitran=cfitran,$
                                      fitord=fitord,no_dust=no_dust,$
                                      quiet=quiet,expterms=expterms,$
                                      gaps=gaps,z=z

  if ~keyword_set(fitord) then fitord=3
  if ~keyword_set(expterms) then expterms=0
  if ~keyword_set(gaps) then $
     print,'GMOS_FIT_FLAT_CONTINUUM: ERROR.  Gap wavelengths not specified.'

  gapcent=[mean(gaps[*,0]),mean(gaps[*,1])]

  ilambda=lambda[index]
  iflux=flux[index]
  iweight=weight[index]
  err = 1d/sqrt(iweight)

  npar = 2+(fitord+expterms)*2
  param = dblarr(npar)

  fcn = 'gmos_flat_continuum_gaps'
  fcnargs = {fitord:fitord,$
             expterms:expterms,$
             gapcent:gapcent}
  parinfo = REPLICATE({$
            limited:[1b,0b],$
            limits:[0d,0d],$
            fixed:0b,$
            mpprint:0b $
                      },npar)

  ichipgap1 = where(ilambda gt gaps[0,0] AND ilambda lt gaps[1,0])
  ichipgap2 = where(ilambda gt gaps[0,1] AND ilambda lt gaps[1,1])
  ; if max(iweight(ichipgap1)) gt min(iweight) then $
  ;    parinfo[npar-2].fixed = 1b
  ; if max(iweight(ichipgap2)) gt min(iweight) then $
  ;    parinfo[npar-1].fixed = 1b
  parinfo[npar-2].limited[0]=0b
  parinfo[npar-1].limited[0]=0b

  fluxfit = mpcurvefit(ilambda,iflux,iweight,param,function_name=fcn,$
                       functargs=fcnargs,/noderiv,parinfo=parinfo,$
                       /quiet)

  gmos_flat_continuum_gaps,lambda,param,continuum,expterms=expterms,$
                           fitord=fitord,gapcent=gapcent

  o1_fitord = 2
  o1_ran_rest = [6275,6400]
  o1_ran = gmos_redshift_spec(o1_ran_rest,z)
  i_o1 = where(lambda ge o1_ran[0] AND lambda le o1_ran[1])
  ii_o1 = where(ilambda ge o1_ran[0] AND ilambda le o1_ran[1])
  parinfo = replicate({value:0d},o1_fitord+1)
  o1_pars = mpfitfun('poly',ilambda[ii_o1],$
                     iflux[ii_o1]-fluxfit[ii_o1],$
                     1d/sqrt(iweight[ii_o1]),$
                     parinfo=parinfo,/quiet)
  continuum[i_o1] += poly(lambda[i_o1],o1_pars)

  ct_coeff=param

  return,continuum

end
