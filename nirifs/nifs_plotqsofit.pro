;
; History
;   13mar05  DSNR  created
;

pro nifs_plotqsofit,instr,outfile,ps=ps,zbuf=zbuf,$
                    qsocntargs=qsocntargs

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

  qsotmpspec = readspec(qsocntargs.qsotmp)
  qsotmpwave = qsotmpspec[*,0]
  qsotmpflux = qsotmpspec[*,1]
  iqsotmpflux = interpol(qsotmpflux,qsotmpwave,wave)
  nifs_qso_cnt_fcn,wave,instr.ct_coeff,modqso,$
                   fitord=qsocntargs.fitord,qsoflux=iqsotmpflux,$
                   qsoord=qsocntargs.qsoord,expterms=qsocntargs.expterms,$
                   qsoonly=qsocntargs.qsoonly
  modresid = modstars - modqso
  
  norm = max(modstars)
  spectot /= norm
  specstars /= norm
  speclines /= norm
  modtot /= norm
  modstars /= norm
  modqso /= norm
  modresid /= norm

  xran = instr.fitrange
  dxran = xran[1] - xran[0]
  xran1 = [xran[0],xran[0]+dxran/3d]
  xran2 = [xran[0]+dxran/3d,xran[0]+2d*dxran/3d]
  xran3 = [xran[0]+2d*dxran/3d,xran[1]]
  i1 = where(wave gt xran1[0] AND wave lt xran1[1],ct1)
  i2 = where(wave gt xran2[0] AND wave lt xran2[1],ct2)
  i3 = where(wave gt xran3[0] AND wave lt xran3[1],ct3)
  
  maxthresh=0.2
  ntop = 20
  ntop++
  nbottom = 20
  nbottom--
  
  loadct,0,/silent
  multiplot,[1,3],/doyaxis,/doxaxis,ygap=0.02
  ydat = specstars
  ymod = modstars
  yran = [min([ydat[i1],ymod[i1],modresid[i1]]),$
          max([ydat[i1],ymod[i1],modresid[i1]])]
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
  plot,wave,ydat,xran=xran1,yran=yran,/xsty,/ysty
  loadct,13,/silent
  oplot,wave,ymod,color=255
;  oplot,wave,modqso,color=75,linesty=2
  oplot,wave,modresid,color=125,linesty=2
  loadct,0,/silent
  multiplot,/doyaxis,/doxaxis
  ydat = specstars
  ymod = modstars
  yran = [min([ydat[i2],ymod[i2],modresid[i2]]),$
          max([ydat[i2],ymod[i2],modresid[i2]])]
  ydi = ydat[i2]
  ymodi = ymod[i2]
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
  plot,wave,ydat,xran=xran2,yran=yran,/xsty,/ysty
  loadct,13,/silent
  oplot,wave,ymod,color=255
;  oplot,wave,modqso,color=75,linesty=2
  oplot,wave,modresid,color=125,linesty=2
  loadct,0,/silent
  if ct3 gt 0 then begin
     multiplot,/doyaxis,/doxaxis
     ydat = specstars
     ymod = modstars
     yran = [min([ydat[i3],ymod[i3],modresid[i3]]),$
             max([ydat[i3],ymod[i3],modresid[i3]])]
     ydi = ydat[i3]
     ymodi = ymod[i3]
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
     plot,wave,ydat,xran=xran3,yran=yran,/xsty,/ysty
     loadct,13,/silent
     oplot,wave,ymod,color=255
;    oplot,wave,modqso,color=75,linesty=2
     oplot,wave,modresid,color=125,linesty=2
     loadct,0,/silent
  endif else multiplot
  multiplot,/reset

  tmpfile = outfile
  if (dops) then device,/close_file $
  else img = cgsnapshot(filename=tmpfile,/jpeg,/nodialog,quality=100)
 
end
