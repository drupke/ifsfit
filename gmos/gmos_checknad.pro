;
; History
;  10jul22  DSNR  created
;
; When fitting NaD line with 2 components, check that the fit is a
; "good" one.

function gmos_checknad,specfile,parfile,fitran

  readcol300,specfile,specwave,specflux,specerr,/silent,/skip,$
             format='(D,D,D)'
  gmos_readnadpars,parfile,abspars,empars,opars

  cntwave = where((specwave ge fitran[0]-20d AND $
                   specwave le fitran[0]) OR $
                  (specwave ge fitran[1] AND $
                   specwave le fitran[1]+20d))
  sigcnt = 1d*stddev(specflux[cntwave])
;;   cntwave = where(specwave ge fitran[0] AND $
;;                   specwave le fitran[1])
;;   sigcnt = mean(specerr[cntwave])

; Bad component =  FWHM too big
;               OR Peak flux doesn't fall below continuum RMS
  badcomp=0
  if n_elements(opars.nabs gt 1) then begin
     for i=0,opars.nabs-1 do $
        if (abspars[3,i]*2d*sqrt(alog(2)) gt 1250d OR $
            1d -(1d -abspars[0,i]*(1d -exp(-abspars[1,i]))) le sigcnt) $
        then badcomp=i+1
  endif
     
  return,badcomp

end
