;
; History
;  10jul19  DSNR  created
;

pro gmos_plotnad,instr,outfile,wavenorm,fluxnorm,parnorm,$
                 ps=ps,zbuf=zbuf,qsocntargs=qsocntargs

  cleanplot,/silent

  dops=0
  dozbuf=0
  if (keyword_set(ps)) then dops=1
  if (keyword_set(zbuf)) then dozbuf=1

  if (dops) then begin
     set_plot,'ps',/copy,/interpolate
     device,filename=outfile+'.eps',/encapsulated,xsize=10,ysize=7.5,$
            bits_per_pixel=8,/color,/inches
     !P.charsize=1
     !P.charthick=1
  endif else if (dozbuf) then begin
     set_plot,'Z'
     device,decomposed=0,set_resolution=[1280,960],set_pixel_depth=24
     !P.charsize=1
     !P.charthick=1
     erase
  endif else begin
     set_plot,'X'
     device,decomposed=0
     window,xsize=1280,ysize=960,xpos=0,ypos=0,retain=2
     !P.charsize=1
     !P.charthick=1
  endelse

  wave = instr.wave
  spectot = instr.spec
  specstars = instr.spec - instr.specfit
  speclines = instr.spec_nocnt
  modtot = instr.specfit + (instr.spec - instr.spec_nocnt)
  modstars = instr.spec - instr.spec_nocnt

  if keyword_set(qsocntargs) then begin
     qsotmpspec = readspec(qsocntargs.qsotmp)
     qsotmpwave = qsotmpspec[*,0]
     qsotmpflux = qsotmpspec[*,1]
     iqsotmpflux = interpol(qsotmpflux,qsotmpwave,wave)
     gmos_qso_cnt_fcn,wave,instr.ct_coeff,modqso,$
                      fitord=qsocntargs.fitord,qsoflux=iqsotmpflux,$
                      qsoord=qsocntargs.qsoord,$
                      expterms=qsocntargs.expterms,/qsoonly
     modresid = modstars - modqso
     specresid = specstars - modqso
  endif else begin
     modresid = modstars
     specresid = specstars
  endelse
  
  norm = max(modstars)
  spectot /= norm
  specstars /= norm
  speclines /= norm
  specresid /= norm
  modtot /= norm
  modstars /= norm
  modresid /= norm

  xran_rest = [5793,5993]
  nad1_rest = 5895.92d
  nad2_rest = 5889.95d
  he_rest = 5875.661d
; Mrk 231
  ; cfit1ran_rest = [5800d,5865d]
  ; cfit2ran_rest = [5900d,5980d]
  cfit1ran_rest = [5800d,5865d]
  cfit2ran_rest = [5905d,5980d]
  cfit1ran = gmos_redshift_spec(cfit1ran_rest,instr.z)
  cfit2ran = gmos_redshift_spec(cfit2ran_rest,instr.z)
  xran1 = gmos_redshift_spec(xran_rest,instr.z)
  nad1 = gmos_redshift_spec(nad1_rest,instr.z)
  nad2 = gmos_redshift_spec(nad2_rest,instr.z)
  he = gmos_redshift_spec(he_rest,instr.z)
  i1 = where(wave gt xran1[0] AND wave lt xran1[1],ct1)
  
  maxthresh=0.2
  ntop = 20
  ntop++
  nbottom = 20
  nbottom--

  ydat = specresid
  ymod = modresid
  yran = [min([ydat[i1],ymod[i1]]),max([ydat[i1],ymod[i1]])]
  ydi = ydat[i1]
  ymodi = ymod[i1]
  y = [ydi-ymodi]
  ny = n_elements(y)
  iysort = sort(y)
  ysort = y[iysort]
  ymodisort = ymodi[iysort]
  if ysort[ny-ntop] lt ysort[ny-1]*maxthresh then $
     yran[1] = max(ysort[0:ny-ntop]+ymodisort[0:ny-ntop])
  if ysort[nbottom] gt ysort[0]*maxthresh then $
     yran[0] = min(ysort[nbottom:ny-1]+ymodisort[nbottom:ny-1])
  if (yran[0] lt 0) then yran[0]=0
  yran[0] = 0d

  multiplot,[1,2],gap=0
  loadct,0,/silent
  plot,wave,specresid,xran=xran1,yran=yran,/xsty,/ysty
  loadct,13,/silent
  oplot,wave,poly(wave,parnorm)/norm,color=255
  oplot,[nad1,nad1],yran,color=125,linesty=2
  oplot,[nad2,nad2],yran,color=125,linesty=2
  oplot,[he,he],yran,color=125,linesty=2
  ytmp = yran[0]+(yran[1]-yran[0])*0.05
  oplot,cfit1ran,[ytmp,ytmp],color=75,linesty=2
  oplot,cfit2ran,[ytmp,ytmp],color=75,linesty=2
  
  multiplot

  yran = [0,1.3]
  loadct,0,/silent
  plot,wavenorm,fluxnorm,xran=xran1,yran=yran,/xsty,/ysty
  loadct,13,/silent
  oplot,xran1,[1,1],color=255
  oplot,[nad1,nad1],yran,color=125,linesty=2
  oplot,[nad2,nad2],yran,color=125,linesty=2
  oplot,[he,he],yran,color=125,linesty=2
  multiplot,/reset

  tmpfile = outfile
  if (dops) then device,/close_file $
  else img = cgsnapshot(filename=tmpfile,/jpeg,/nodialog,quality=100)
 
end
