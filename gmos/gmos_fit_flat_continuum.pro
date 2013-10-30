function gmos_fit_flat_continuum,lambda,flux,weight,template_flux,index,$
                                 ct_coeff,cfitran=cfitran,fitord=fitord,$
                                 no_dust=no_dust,quiet=quiet,$
                                 refito1=refito1,z=z

; History
;  09aug14  DSNR  created
;  10may27  DSNR  cleaned up

  if ~keyword_set(fitord) then fitord=3

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
    
  continuum = polycomp(lambda,fluxfit)

  if keyword_set(refito1) AND keyword_set(z) then begin
     o1_fitord = 2
     o1_ran_rest = [6275,6400]
     o1_ran = gmos_redshift_spec(o1_ran_rest,z)
     i_o1 = where(lambda ge o1_ran[0] AND lambda le o1_ran[1])
     ii_o1 = where(ilambda ge o1_ran[0] AND ilambda le o1_ran[1])
     parinfo = replicate({value:0d},o1_fitord+1)
     o1_pars = mpfitfun('poly',ilambda[ii_o1],$
                        iflux[ii_o1]-yfit[ii_o1],$
                        1d/sqrt(iweight[ii_o1]),$
                        parinfo=parinfo,/quiet)
     continuum[i_o1] += poly(lambda[i_o1],o1_pars)
  endif

  ct_coeff=0d

  return,continuum

end
