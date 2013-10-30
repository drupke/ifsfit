;
; History
;  11jul20  DSNR  copied from gmos_qso_cnt_fcn.pro
;
pro gmos_flat_continuum_gaps,x,p,ymod,fitord=fitord,expterms=expterms,$
                             gapcent=gapcent

  npar = n_elements(p)
  nx = n_elements(x)
  dx = max(x) - min(x)
  xc = x - min(x)
  xc_hi = max(x) - x
  ymod = dblarr(nx)

  if keyword_set(expterms) then nexp=expterms else nexp=0
  itot = fitord+nexp

; Stellar continuum
  for i=0,fitord-1 do $
     ymod += p[i]        * (xc-p[itot+i])^double(i)
  if nexp ge 1 then $
     ymod += p[fitord]   * exp(-(xc-p[itot+fitord])/dx)
  if nexp ge 2 then $
     ymod += p[fitord+1] * (1d - exp(-(xc-p[itot+fitord+1])/dx))
  if nexp ge 3 then $
     ymod += p[fitord+2] * exp(-(xc_hi-p[itot+fitord+2])/dx)
  if nexp eq 4 then $
     ymod += p[fitord+3] * (1d - exp(-(xc_hi-p[itot+fitord+3])/dx))

; Offset terms
  ilo = where(x le gapcent[0])
  ihi = where(x ge gapcent[1])
  ymod[ilo] += p[npar-2]
  ymod[ihi] += p[npar-1]

end
