;
; History
;  11jul28  DSNR  created
;

function gmos_compweq,datfile,parfile,z

  readcol300,datfile,specwave,specflux,specerr,/silent,/skip,$
             format='(D,D,D)'

  gmos_readnadpars,parfile,abspars,empars,opars

  modflux = dblarr(n_elements(specflux))+1d
  modabs = dblarr(n_elements(specflux))+1d
  modem = dblarr(n_elements(specflux))+1d
  for i=0,opars.nabs-1 do modabs *= nad(specwave,abspars[*,i])
  for i=0,opars.nem-1 do begin
     arg = ((specwave-empars[1,i])/(empars[1,i]*empars[2,i]/299792d))^2d
     modem += empars[0,i]*exp(-arg)
  endfor
  modflux = modabs * modem

; Compute weq
  dlam = specwave - shift(specwave,1)
  weq = total((1-modabs)*dlam)

  return,weq

end
