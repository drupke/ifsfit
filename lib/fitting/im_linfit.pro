function linefunc, x, p
; linear model
return, p[0] + p[1] * x
end    

function im_linfit, x, y, invvar=invvar, fixed=fixed, yfit=yfit
; jm02jul18uofa
; fit a line to the data but with the option of constraining either
; the slope or the intercept

; fixed = 0B for free and 1B for fixed

    if n_elements(invvar) eq 0L then invvar = y*0.0+1.0
    if n_elements(fixed) eq 0L then fixed = [0B,0B]
    
    parinfo = {value: 0.0D, $
               fixed: 0L}
;               value: 0.0D,    $
;               fixed: 0L,      $
;               tied: '',       $
;               limited: [0,0], $
;               limits: [0.0D,0.0D]}
    parinfo = replicate(parinfo,2)
    parinfo.fixed = fixed

    coeff = mpfitfun('linefunc',x,y,weights=invvar,parinfo=parinfo,functargs=functargs,$
                     covar=covar,perror=perror,yfit=yfit,nfev=nfev,niter=niter,/quiet)

return, coeff
end    
