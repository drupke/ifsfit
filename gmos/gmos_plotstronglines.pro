pro gmos_plotstronglines,instr,outfile,ps=ps,zbuf=zbuf

dops=0
dozbuf=0
if (keyword_set(ps)) then dops=1
if (keyword_set(zbuf)) then dozbuf=1

if (dops) then begin
  set_plot,'ps',/copy,/interpolate
  device,filename=outfile+'.ps',/encapsulated,xsize=20,ysize=10,bits_per_pixel=8,/color
  !P.charsize=0.75
  !P.charthick=1
endif else if (dozbuf) then begin
  set_plot,'Z'
  device,decomposed=0,set_resolution=[1280,640],set_pixel_depth=24
  !P.charsize=1
  !P.charthick=1
  erase
endif else begin
  set_plot,'X'
  device,decomposed=0
  window,xsize=1280,ysize=640,xpos=0,ypos=0,retain=2
  !P.charsize=1
  !P.charthick=1
endelse

wave = instr.wave
spectot = instr.spec
specstars = instr.spec - instr.specfit
speclines = instr.spec_nocnt
modtot = instr.specfit + (instr.spec - instr.spec_nocnt)
modstars = instr.spec - instr.spec_nocnt
modlines = instr.specfit

norm = max(modstars)
spectot /= norm
specstars /= norm
speclines /= norm
modtot /= norm
modstars /= norm
modlines /= norm

; GMOS
;; xran1 = [6270.01,6379.99]
;; xran2 = [6540.01,6589.99]
; KPNO
xran1 = [6540.01,6589.99]
xran2 = [6690.01,6749.99]
i1 = where(wave gt xran1[0] AND wave lt xran1[1],ct1)
i2 = where(wave gt xran2[0] AND wave lt xran2[1],ct2)

loadct,0,/silent
multiplot,[2,3],/rowmajor,xgap=0.02,/doyaxis
if ct1 gt 0 then begin
   ydat = spectot
   ymod = modtot
   yran = [min([ydat[i1],ymod[i1]]),max([ydat[i1],ymod[i1]])]
   plot,wave,ydat,xran=xran1,yran=yran,/xsty,/ysty
   loadct,13,/silent
   oplot,wave,ymod,color=255
   loadct,0,/silent
   multiplot,/doyaxis
   ydat = specstars
   ymod = modstars
   yran = [min([ydat[i1],ymod[i1]]),max([ydat[i1],ymod[i1]])]
   plot,wave,ydat,xran=xran1,yran=yran,/xsty,/ysty
   loadct,13,/silent
   oplot,wave,ymod,color=255
   loadct,0,/silent
   multiplot,/doyaxis
   ydat = speclines
   ymod = modlines
   yran = [min([ydat[i1],ymod[i1]]),max([ydat[i1],ymod[i1]])]
   plot,wave,ydat,xran=xran1,yran=yran,/xsty,/ysty
   loadct,13,/silent
   oplot,wave,ymod,color=255
endif else begin
   multiplot
   multiplot
endelse

loadct,0,/silent
multiplot,/doyaxis
ydat = spectot
ymod = modtot
yran = [min([ydat[i2],ymod[i2]]),max([ydat[i2],ymod[i2]])]
plot,wave,ydat,xran=xran2,yran=yran,/xsty,/ysty
loadct,13,/silent
oplot,wave,ymod,color=255
loadct,0,/silent
multiplot,/doyaxis
ydat = specstars
ymod = modstars
yran = [min([ydat[i2],ymod[i2]]),max([ydat[i2],ymod[i2]])]
plot,wave,ydat,xran=xran2,yran=yran,/xsty,/ysty
loadct,13,/silent
oplot,wave,ymod,color=255
loadct,0,/silent
multiplot,/doyaxis
ydat = speclines
ymod = modlines
yran = [min([ydat[i2],ymod[i2]]),max([ydat[i2],ymod[i2]])]
plot,wave,ydat,xran=xran2,yran=yran,/xsty,/ysty
loadct,13,/silent
oplot,wave,ymod,color=255

loadct,0,/silent
multiplot,/reset

tmpfile = outfile
if (dops) then device,/close_file $
else img = tvread(filename=tmpfile,/jpeg,/nodialog,quality=100)
 
end
