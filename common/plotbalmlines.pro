;
; History
;  09xxxYY  DSNR  created
;

pro plotbalmlines,instr,outfile,ps=ps,zbuf=zbuf

dops=0
dozbuf=0
if keyword_set(ps) then dops=1
if keyword_set(zbuf) then dozbuf=1

if (dops) then begin
   set_plot,'ps',/copy,/interpolate
   device,filename=outfile+'.eps',/encapsulated,/inches,$
          xsize=10,ysize=5,bits_per_pixel=8,/color
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

bx5 = 6563d * (1d + instr.z.star)
bx4 = 4861d * (1d + instr.z.star)
bx3 = 4340d * (1d + instr.z.star)
bx2 = 4104d * (1d + instr.z.star)
bx1 = 3970d * (1d + instr.z.star)
bx0 = 3889d * (1d + instr.z.star)

halfwidth = 10d * (1d + instr.z.star)

xran = [bx0-halfwidth,bx0+halfwidth]
ix = where(wave gt xran[0] AND wave lt xran[1],ctix)
loadct,0,/silent
multiplot,[6,3],/rowmajor
if (ctix gt 0) then begin
  ydat = spectot
  ymod = modtot
  yran = [min([ydat[ix],ymod[ix]]),max([ydat[ix],ymod[ix]])]
  plot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty
  loadct,13,/silent
  oplot,wave,ymod,color=255
  xyouts,xran[0]+(xran[1]-xran[0])*0.9,yran[0]+(yran[1]-yran[0])*0.9,$
         textoidl('H\zeta'),color=255
  loadct,0,/silent
  multiplot
  ydat = specstars
  ymod = modstars
  yran = [min([ydat[ix],ymod[ix]]),max([ydat[ix],ymod[ix]])]
  plot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty
  loadct,13,/silent
  oplot,wave,ymod,color=255
  loadct,0,/silent
  multiplot
  ydat = speclines
  ymod = modlines
  yran = [min([ydat[ix],ymod[ix]]),max([ydat[ix],ymod[ix]])]
  plot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty
  loadct,13,/silent
  oplot,wave,ymod,color=255
endif else begin
  multiplot
  multiplot
endelse

xran = [bx1-halfwidth,bx1+halfwidth]
ix = where(wave gt xran[0] AND wave lt xran[1])
loadct,0,/silent
multiplot
ydat = spectot
ymod = modtot
yran = [min([ydat[ix],ymod[ix]]),max([ydat[ix],ymod[ix]])]
plot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty
loadct,13,/silent
oplot,wave,ymod,color=255
xyouts,xran[0]+(xran[1]-xran[0])*0.9,yran[0]+(yran[1]-yran[0])*0.9,$
       textoidl('H\epsilon'),color=255
loadct,0,/silent
multiplot
ydat = specstars
ymod = modstars
yran = [min([ydat[ix],ymod[ix]]),max([ydat[ix],ymod[ix]])]
plot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty
loadct,13,/silent
oplot,wave,ymod,color=255
loadct,0,/silent
multiplot
ydat = speclines
ymod = modlines
yran = [min([ydat[ix],ymod[ix]]),max([ydat[ix],ymod[ix]])]
plot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty
loadct,13,/silent
oplot,wave,ymod,color=255

xran = [bx2-halfwidth,bx2+halfwidth]
ix = where(wave gt xran[0] AND wave lt xran[1])
loadct,0,/silent
multiplot
ydat = spectot
ymod = modtot
yran = [min([ydat[ix],ymod[ix]]),max([ydat[ix],ymod[ix]])]
plot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty
loadct,13,/silent
oplot,wave,ymod,color=255
xyouts,xran[0]+(xran[1]-xran[0])*0.9,yran[0]+(yran[1]-yran[0])*0.9,$
       textoidl('H\delta'),color=255
loadct,0,/silent
multiplot
ydat = specstars
ymod = modstars
yran = [min([ydat[ix],ymod[ix]]),max([ydat[ix],ymod[ix]])]
plot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty
loadct,13,/silent
oplot,wave,ymod,color=255
loadct,0,/silent
multiplot
ydat = speclines
ymod = modlines
yran = [min([ydat[ix],ymod[ix]]),max([ydat[ix],ymod[ix]])]
plot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty
loadct,13,/silent
oplot,wave,ymod,color=255

xran = [bx3-halfwidth,bx3+halfwidth]
ix = where(wave gt xran[0] AND wave lt xran[1])
loadct,0,/silent
multiplot
ydat = spectot
ymod = modtot
yran = [min([ydat[ix],ymod[ix]]),max([ydat[ix],ymod[ix]])]
plot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty
loadct,13,/silent
oplot,wave,ymod,color=255
xyouts,xran[0]+(xran[1]-xran[0])*0.9,yran[0]+(yran[1]-yran[0])*0.9,$
       textoidl('H\gamma'),color=255
loadct,0,/silent
multiplot
ydat = specstars
ymod = modstars
yran = [min([ydat[ix],ymod[ix]]),max([ydat[ix],ymod[ix]])]
plot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty
loadct,13,/silent
oplot,wave,ymod,color=255
loadct,0,/silent
multiplot
ydat = speclines
ymod = modlines
yran = [min([ydat[ix],ymod[ix]]),max([ydat[ix],ymod[ix]])]
plot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty
loadct,13,/silent
oplot,wave,ymod,color=255

xran = [bx4-halfwidth,bx4+halfwidth]
ix = where(wave gt xran[0] AND wave lt xran[1])
loadct,0,/silent
multiplot
ydat = spectot
ymod = modtot
yran = [min([ydat[ix],ymod[ix]]),max([ydat[ix],ymod[ix]])]
plot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty
loadct,13,/silent
oplot,wave,ymod,color=255
xyouts,xran[0]+(xran[1]-xran[0])*0.9,yran[0]+(yran[1]-yran[0])*0.9,$
       textoidl('H\beta'),color=255
loadct,0,/silent
multiplot
ydat = specstars
ymod = modstars
yran = [min([ydat[ix],ymod[ix]]),max([ydat[ix],ymod[ix]])]
plot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty
loadct,13,/silent
oplot,wave,ymod,color=255
loadct,0,/silent
multiplot
ydat = speclines
ymod = modlines
yran = [min([ydat[ix],ymod[ix]]),max([ydat[ix],ymod[ix]])]
plot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty
loadct,13,/silent
oplot,wave,ymod,color=255

xran = [bx5-halfwidth,bx5+halfwidth]
ix = where(wave gt xran[0] AND wave lt xran[1],ctix)
loadct,0,/silent
multiplot
if ctix gt 0 then begin
  ydat = spectot
  ymod = modtot
  yran = [min([ydat[ix],ymod[ix]]),max([ydat[ix],ymod[ix]])]
  plot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty
  loadct,13,/silent
  oplot,wave,ymod,color=255
  xyouts,xran[0]+(xran[1]-xran[0])*0.9,yran[0]+(yran[1]-yran[0])*0.9,$
         textoidl('H\alpha'),color=255
  loadct,0,/silent
  multiplot
  ydat = specstars
  ymod = modstars
  yran = [min([ydat[ix],ymod[ix]]),max([ydat[ix],ymod[ix]])]
  plot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty
  loadct,13,/silent
  oplot,wave,ymod,color=255
  loadct,0,/silent
  multiplot
  ydat = speclines
  ymod = modlines
  yran = [min([ydat[ix],ymod[ix]]),max([ydat[ix],ymod[ix]])]
  plot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty
  loadct,13,/silent
  oplot,wave,ymod,color=255
endif else begin
  multiplot
  multiplot
endelse

loadct,0,/silent
multiplot,/reset

tmpfile = outfile
if (dops) then device,/close_file $
else img = tvread(filename=tmpfile,/jpeg,/nodialog,quality=100)
 
end
