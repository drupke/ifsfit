;
; History
;  10jan19  DSNR  created
;
; Refine initial guess for redshift by finding emission-line peak
;
function uhsf_refinez,zinit,lambda,flux,lref=lref,searchwidth=searchwidth

  if ~ keyword_set(searchwidth) then searchwidth=10d
  if ~ keyword_set(lref) then lref = 6562.8d
  
  searchwindow = [lref-searchwidth,lref+searchwidth]
  iref0 = value_locate(lambda,searchwindow[0])
  iref1 = value_locate(lambda,searchwindow[1])
  lambdawindow = lambda[iref0:iref1]
  fluxwindow = flux[iref0:iref1]
  fmax = max(fluxwindow,imax)
  znew = zinit + (lambdawindow[imax] / lref) - 1d
  
  return,znew
  
end
