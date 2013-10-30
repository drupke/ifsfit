;
; History
;  10jul17  DSNR  created
;


function gmos_fit_qso_continuum,lambda,flux,weight,template_flux,index,$
                                ct_coeff,fitord=fitord,qsotmp=qsotmp,$
                                qsoord=qsoord,z=z,expterms=expterms,$
                                no_dust=no_dust,quiet=quiet

  if ~keyword_set(fitord) then fitord=3
  if ~keyword_set(expterms) then expterms=0

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

  fcn = 'gmos_qso_cnt_fcn'
  fcnargs = {fitord:fitord,$
             qsoflux:iiqsotmpflux,$
             qsoord:qsoord,$
             expterms:expterms}
  parinfo = REPLICATE({$
            limited:[1b,0b],$
            limits:[0d,0d],$
            fixed:0b,$
            mpprint:0b $
                      },npar)

  ichipgap1 = where(ilambda gt 6020 AND ilambda lt 6040)
  ichipgap2 = where(ilambda gt 6520 AND ilambda lt 6540)
  if max(iweight(ichipgap1)) gt min(iweight) then $
     parinfo[npar-2].fixed = 1b
  if max(iweight(ichipgap2)) gt min(iweight) then $
     parinfo[npar-1].fixed = 1b
  parinfo[npar-2].limited[0]=0b
  parinfo[npar-1].limited[0]=0b

  fluxfit = mpcurvefit(ilambda,iflux,iweight,param,function_name=fcn,$
                       functargs=fcnargs,/noderiv,parinfo=parinfo,$
                       /quiet)


  gmos_qso_cnt_fcn,lambda,param,continuum,expterms=expterms,$
                   fitord=fitord,qsoflux=iqsotmpflux,qsoord=qsoord

; Re-fit [OI] and Ha/[NII] regions separately
  
; Fitting a Gaussian to Halpha absorption doesn't work well because
; there's just not enough spectra to do it robustly.  Ideally, should
; fit at the same time one fits emission.

  ; n2ha_ran_rest = [6460,6660]
  ; n2ha_ran = gmos_redshift_spec(n2ha_ran_rest,z) 
  ; i_n2ha = where(lambda ge n2ha_ran[0] AND lambda le n2ha_ran[1])
  ; ii_n2ha = where(ilambda ge n2ha_ran[0] AND ilambda le n2ha_ran[1])
  ; parinfo = replicate({value:0d,fixed:0d},n2ha_fitord+1)
  ; n2ha_pars = mpfitfun('poly',ilambda[i_n2ha],$
  ;                      iflux[i_n2ha]-fluxfit[i_n2ha],$
  ;                      1d/sqrt(iweight[i_n2ha]),$
  ;                      parinfo=parinfo,/quiet)
  ; continuum[ii_n2ha] += poly(lambda[ii_n2ha],n2ha_pars)

  ; ha_rest = 6563.4d
  ; z_sys = {star:0.0422d,gas:0.0422d}
  ; ha = gmos_redshift_spec(ha_rest,z_sys)
  ; gmos_qso_cnt_fcn,ilambda,param,icont,expterms=expterms,$
  ;                  fitord=fitord,qsoflux=iiqsotmpflux,qsoord=qsoord
  ; gmos_qso_cnt_fcn,ilambda,param,icont_qso,expterms=expterms,$
  ;                  fitord=fitord,qsoflux=iiqsotmpflux,qsoord=qsoord,$
  ;                  /qsoonly
  ; istel_n2ha = iflux[ii_n2ha]-icont_qso[ii_n2ha]
  ; istel_n2ha_norm = istel_n2ha / max(icont[ii_n2ha]-icont_qso[ii_n2ha])
  ; parinfo = replicate({value:0d,fixed:0b,limited:bytarr(2),$
  ;                      limits:dblarr(2)},3)
  ; parinfo[0].limited = [1b,1b]
  ; parinfo[0].limits = [-1d,0d]
  ; parinfo[1].fixed=1b
  ; parinfo[1].value=ha
  ; parinfo[2].fixed=1b
  ; parinfo[2].value=6d
  ; ha_pars = mpfitfun('gaussian',ilambda[ii_n2ha],$
  ;                    istel_n2ha_norm,$
  ;                    1d/sqrt(iweight[ii_n2ha]),$
  ;                    parinfo=parinfo,/quiet)
  ; ha_pars[0]*= max(istel_n2ha)
  ; continuum[i_n2ha] += gaussian(lambda[i_n2ha],ha_pars)
  
  o1_fitord = 1
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
