;
; History
;  10jul21  DSNR  created
;

pro uhsf_pltnadfit,datfile,parfile,outfile,z,oplot=oplot,xran=xran,$
                   yran=yran,weq=weq

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

; PLOT

  nad1_rest = 5895.92d
  nad2_rest = 5889.95d
  he_rest = 5875.661d
  nad1 = gmos_redshift_spec(nad1_rest,z)
  nad2 = gmos_redshift_spec(nad2_rest,z)
  he = gmos_redshift_spec(he_rest,z)
  
  if ~ keyword_set(oplot) then begin
     set_plot,'Z'
     device,decomposed=0,set_resolution=[1280,960],set_pixel_depth=24
     !P.charsize=1
     !P.charthick=1
  endif

  if ~ keyword_set(xran) then xran = [opars.llo-20d,opars.lhi+20d]
  if ~ keyword_set(yran) then yran = [0d,1.5d]

  loadct,0,/silent
  plot,specwave,specflux,xran=xran,yran=yran,/xsty,/ysty
  loadct,13,/silent
  oplot,specwave,modabs,color=75,thick=4
  oplot,specwave,modflux,color=255,thick=4
  if ~ keyword_set(oplot) then oplot,specwave,specerr,color=75
  if ~ keyword_set(oplot) then restcolor=125 else restcolor=0
  oplot,[nad1,nad1],yran,color=restcolor,linesty=2
  oplot,[nad2,nad2],yran,color=restcolor,linesty=2
  oplot,[he,he],yran,color=restcolor,linesty=2

  if ~ keyword_set(oplot) then begin
     tmpfile = outfile
     img = cgsnapshot(filename=tmpfile,/jpeg,/nodialog,quality=100)
  endif

end
