;
; History
;  13mar05  DSNR  created
;

pro nifs_qso_cnt_fcn,x,p,ymod,fitord=fitord,$
                     qsoflux=qsoflux,qsoord=qsoord,$
                     qsoonly=qsoonly,expterms=expterms

  npar = n_elements(p)
  nx = n_elements(x)
  dx = max(x) - min(x)
  xc = x - min(x)
  xc_hi = max(x) - x

  ymod = dblarr(nx)
  qsoscl = dblarr(nx)

  if keyword_set(expterms) then nexp=expterms else nexp=0
  itot = fitord+qsoord+2*nexp

; Stellar continuum
  if ~ keyword_set(qsoonly) then begin
     for i=0,fitord-1 do $
        ymod += p[i]        * (xc-p[itot+i])^double(i)
     if nexp ge 1 then $
        ymod += p[fitord]   * exp(-(xc-p[itot+fitord])/dx)
     if nexp ge 2 then $
        ymod += p[fitord+2] * (1d - exp(-(xc-p[itot+fitord+2])/dx))
     if nexp ge 3 then $
        ymod += p[fitord+1] * exp(-(xc_hi-p[itot+fitord+1])/dx)
     if nexp eq 4 then $
        ymod += p[fitord+3] * (1d - exp(-(xc_hi-p[itot+fitord+3])/dx))
  endif

; QSO continuum
  for i=0,qsoord-1 do $
     qsoscl += p[fitord+nexp+i] * $
                (xc-p[itot+fitord+nexp+i])^double(i)
  if nexp ge 1 then $
     qsoscl += p[fitord+nexp+qsoord] * $
               exp(-(xc-p[itot+fitord+nexp+qsoord])/dx)
  if nexp ge 2 then $
     qsoscl += p[fitord+nexp+qsoord+2] * $
               (1d - exp(-(xc-p[itot+fitord+nexp+qsoord+2])/dx)) 
  if nexp ge 3 then $
     qsoscl += p[fitord+nexp+qsoord+1] * $
               exp(-(xc_hi-p[itot+fitord+nexp+qsoord+1])/dx)
  if nexp eq 4 then $
     qsoscl += p[fitord+nexp+qsoord+3] * $
               (1d - exp(-(xc_hi-p[itot+fitord+nexp+qsoord+3])/dx))
  ymod += qsoscl*qsoflux

end
