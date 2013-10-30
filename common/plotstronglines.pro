;
; History
;  10jan20  DSNR  modified to plot Hb/[OIII] separately and
;                 to show [OI] and [SII] lines
;
; Plot strong emission line fits.
;

function componeline,instr,line,comp,velsig=velsig
  
  iline = where(instr.linelabel eq line,ct)
  ppoff = instr.param[0]
  ncomp = instr.param[1]
  ppoff0 = ppoff - (ncomp-1)

  nline = n_elements(instr.linelabel)
  indices = ppoff+(comp-1)*nline*3+iline*3
  indices = indices[0] + indgen(3)
  gausspar = instr.param[indices]
  if keyword_set(velsig) then gausspar[2] *= gausspar[1]/299792d
;  print,comp,indices[0],instr.param[indices]
  flux = gaussian(instr.wave,gausspar,/double)

  return,flux

end

pro plotstronglines,instr,outfile,ps=ps,zbuf=zbuf,comp=comp,$
                    plotilines=plotilines,velsig=velsig

  dops=0
  dozbuf=0
  if keyword_set(ps) then dops=1
  if keyword_set(zbuf) then dozbuf=1
  if ~ keyword_set(comp) then comp=1
  if ~ keyword_set(velsig) then velsig=0

  if (dops) then begin
     set_plot,'ps',/copy,/interpolate
     device,filename=outfile+'.eps',/encapsulated,/inches,$
            xsize=10,ysize=7.5,bits_per_pixel=8,/color
     !P.charsize=1
     !P.charthick=2
     !P.thick=2
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

  defaultXtickint=!X.tickinterval
  defaultXminor=!X.minor
  !X.tickinterval=20
  !X.minor=10

  if keyword_set(plotilines) then ncomp = instr.param[1]

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

  if keyword_set(comp) then zbase = instr.z.gas[comp-1] $
  else zbase = instr.z.star

  lab = textoidl(['[OII]','H\beta','[OIII]','[OI]','H\alpha/[NII]','[SII]'])
  off = [-15d,15d]
  xran1 = (3727.4d + off) * (1d + zbase)
  xran2 = (4861.3d + off) * (1d + zbase)
  xran3 = ([4958.9d,5006.8d]+off) * (1d + zbase)
  xran4 = ([6300.3d,6363.7d]+off) * (1d + zbase)
  xran5 = ([6548.0d,6583.4d]+off) * (1d + zbase)
  xran6 = ([6716.4d,6730.8d]+off) * (1d + zbase)
  i1 = where(wave gt xran1[0] AND wave lt xran1[1],ct1)
  i2 = where(wave gt xran2[0] AND wave lt xran2[1],ct2)
  i3 = where(wave gt xran3[0] AND wave lt xran3[1],ct3)
  i4 = where(wave gt xran4[0] AND wave lt xran4[1],ct4)
  i5 = where(wave gt xran5[0] AND wave lt xran5[1],ct5)
  i6 = where(wave gt xran6[0] AND wave lt xran6[1],ct6)

  loadct,0,/silent
  !X.tickinterval=10
  !X.minor=5
  if ct1 gt 0 then begin
     !P.position = [0.1,0.75,0.35,0.95]
     ydat = spectot
     ymod = modtot
     yran = [min([ydat[i1],ymod[i1]]),max([ydat[i1],ymod[i1]])]
     plot,wave,ydat,xran=xran1,yran=yran,/xsty,/ysty,pos=!P.position,$
          xtickn=replicate(' ',5),ytit='Fit'
     loadct,13,/silent
     oplot,wave,ymod,color=255
     if keyword_set(plotilines) then begin
        for i=1,ncomp do begin
           flux = componeline(instr,'[OII]3726',i,velsig=velsig)
           oplot,wave,yran[0]+flux/norm,color=75,linesty=2
           flux = componeline(instr,'[OII]3729',i,velsig=velsig)
           oplot,wave,yran[0]+flux/norm,color=75,linesty=2
        endfor
     endif
     loadct,0,/silent
     xyouts,xran1[0]+(xran1[1]-xran1[0])*0.05d,$
            yran[0]+(yran[1]-yran[0])*0.85d,$
            lab[0],charsize=1.5,charthick=2
;     multiplot,/doyaxis,/doxaxis
     !P.position = [0.1,0.55,0.35,0.75]
     ydat = specstars
     ymod = modstars
     yran = [min([ydat[i1],ymod[i1]]),max([ydat[i1],ymod[i1]])]
     plot,wave,ydat,xran=xran1,yran=yran,/xsty,/ysty,$
          pos=!P.position,/noerase,ytit='Residual'
     loadct,13,/silent
     oplot,wave,ymod,color=255
  endif
  !X.minor=10
  !X.tickinterval=20
;;   endif else begin
;;      multiplot
;;   endelse

  
  loadct,0,/silent
  !X.tickinterval=10
  !X.minor=5
  if ct2 gt 0 then begin
;     multiplot,/doyaxis
     !P.position = [0.4,0.75,0.65,0.95]
     ydat = spectot
     ymod = modtot
     yran = [min([ydat[i2],ymod[i2]]),max([ydat[i2],ymod[i2]])]
     plot,wave,ydat,xran=xran2,yran=yran,/xsty,/ysty,pos=!P.position,$
          xtickn=replicate(' ',5),/noerase
     loadct,13,/silent
     oplot,wave,ymod,color=255
     if keyword_set(plotilines) then begin
        for i=1,ncomp do begin
           flux = componeline(instr,'Hbeta',i,velsig=velsig)
           oplot,wave,yran[0]+flux/norm,color=75,linesty=2
        endfor
     endif
     loadct,0,/silent
     xyouts,xran2[0]+(xran2[1]-xran2[0])*0.05d,$
            yran[0]+(yran[1]-yran[0])*0.85d,$
            lab[1],charsize=1.5,charthick=2
;     multiplot,/doyaxis,/doxaxis
     !P.position = [0.4,0.55,0.65,0.75]
     ydat = specstars
     ymod = modstars
     yran = [min([ydat[i2],ymod[i2]]),max([ydat[i2],ymod[i2]])]
     plot,wave,ydat,xran=xran2,yran=yran,/xsty,/ysty,pos=!P.position,$
          /noerase
     loadct,13,/silent
     oplot,wave,ymod,color=255
  endif
  !X.tickinterval=20
  !X.minor=10
;;   endif else begin
;;      multiplot
;;      multiplot
;;   endelse


  loadct,0,/silent
  if ct3 gt 0 then begin
;     multiplot,/doyaxis
     !P.position = [0.7,0.75,0.95,0.95]
     ydat = spectot
     ymod = modtot
     yran = [min([ydat[i3],ymod[i3]]),max([ydat[i3],ymod[i3]])]
     plot,wave,ydat,xran=xran3,yran=yran,/xsty,/ysty,pos=!P.position,$
          xtickn=replicate(' ',5),/noerase
     loadct,13,/silent
     oplot,wave,ymod,color=255
     if keyword_set(plotilines) then begin
        for i=1,ncomp do begin
           flux = componeline(instr,'[OIII]4959',i,velsig=velsig)
           oplot,wave,yran[0]+flux/norm,color=75,linesty=2
           flux = componeline(instr,'[OIII]5007',i,velsig=velsig)
           oplot,wave,yran[0]+flux/norm,color=75,linesty=2
        endfor
     endif
     loadct,0,/silent
     xyouts,xran3[0]+(xran3[1]-xran3[0])*0.05d,$
            yran[0]+(yran[1]-yran[0])*0.85d,$
            lab[2],charsize=1.5,charthick=2
;     multiplot,/doyaxis,/doxaxis
     !P.position = [0.7,0.55,0.95,0.75]
     ydat = specstars
     ymod = modstars
     yran = [min([ydat[i3],ymod[i3]]),max([ydat[i3],ymod[i3]])]
     plot,wave,ydat,xran=xran3,yran=yran,/xsty,/ysty,pos=!P.position,$
          /noerase
     loadct,13,/silent
     oplot,wave,ymod,color=255
  endif
;;   endif else begin
;;      multiplot
;;      multiplot
;;   endelse

  loadct,0,/silent
  if ct4 gt 0 then begin
;     multiplot,/doyaxis
     !P.position = [0.1,0.3,0.35,0.5]
     ydat = spectot
     ymod = modtot
     yran = [min([ydat[i4],ymod[i4]]),max([ydat[i4],ymod[i4]])]
     plot,wave,ydat,xran=xran4,yran=yran,/xsty,/ysty,pos=!P.position,$
          xtickn=replicate(' ',5),/noerase,ytit='Fit'
     loadct,13,/silent
     oplot,wave,ymod,color=255
     if keyword_set(plotilines) then begin
        for i=1,ncomp do begin
           flux = componeline(instr,'[OI]6300',i,velsig=velsig)
           oplot,wave,yran[0]+flux/norm,color=75,linesty=2
           flux = componeline(instr,'[OI]6364',i,velsig=velsig)
           oplot,wave,yran[0]+flux/norm,color=75,linesty=2
        endfor
     endif
     loadct,0,/silent
     xyouts,xran4[1]-(xran4[1]-xran4[0])*0.05d,$
            yran[0]+(yran[1]-yran[0])*0.85d,$
            lab[3],charsize=1.5,charthick=2,align=1
;     multiplot,/doyaxis,/doxaxis
     !P.position = [0.1,0.1,0.35,0.3]
     ydat = specstars
     ymod = modstars
     yran = [min([ydat[i4],ymod[i4]]),max([ydat[i4],ymod[i4]])]
     plot,wave,ydat,xran=xran4,yran=yran,/xsty,/ysty,$
          pos=!P.position,/noerase,ytit='Residual'
     loadct,13,/silent
     oplot,wave,ymod,color=255
  endif
;;   endif else begin
;;      multiplot
;;      multiplot
;;   endelse

  loadct,0,/silent
  if ct5 gt 0 then begin
     !P.position = [0.4,0.3,0.65,0.5]
;     multiplot,/doyaxis
     ydat = spectot
     ymod = modtot
     yran = [min([ydat[i5],ymod[i5]]),max([ydat[i5],ymod[i5]])]
     plot,wave,ydat,xran=xran5,yran=yran,/xsty,/ysty,pos=!P.position,$
          xtickn=replicate(' ',5),/noerase
     loadct,13,/silent
     oplot,wave,ymod,color=255
     if keyword_set(plotilines) then begin
        for i=1,ncomp do begin
           flux = componeline(instr,'[NII]6548',i,velsig=velsig)
           oplot,wave,yran[0]+flux/norm,color=75,linesty=2 
           flux = componeline(instr,'Halpha',i,velsig=velsig)
           oplot,wave,yran[0]+flux/norm,color=75,linesty=2 
           flux = componeline(instr,'[NII]6583',i,velsig=velsig)
           oplot,wave,yran[0]+flux/norm,color=75,linesty=2
        endfor
     endif
     loadct,0,/silent
     xyouts,xran5[0]+(xran5[1]-xran5[0])*0.05d,$
            yran[0]+(yran[1]-yran[0])*0.85d,$
            lab[4],charsize=1.5,charthick=2
;     multiplot,/doyaxis,/doxaxis
     !P.position = [0.4,0.1,0.65,0.3]
     ydat = specstars
     ymod = modstars
     yran = [min([ydat[i5],ymod[i5]]),max([ydat[i5],ymod[i5]])]
     plot,wave,ydat,xran=xran5,yran=yran,/xsty,/ysty,$
          pos=!P.position,/noerase
     loadct,13,/silent
     oplot,wave,ymod,color=255
  endif
;;   endif else begin
;;      multiplot
;;      multiplot
;;   endelse


  loadct,0,/silent
  if ct6 gt 0 then begin
;     multiplot,/doyaxis
     !P.position = [0.7,0.3,0.95,0.5]
     ydat = spectot
     ymod = modtot
     yran = [min([ydat[i6],ymod[i6]]),max([ydat[i6],ymod[i6]])]
     plot,wave,ydat,xran=xran6,yran=yran,/xsty,/ysty,pos=!P.position,$
          xtickn=replicate(' ',5),/noerase
     loadct,13,/silent
     oplot,wave,ymod,color=255
     if keyword_set(plotilines) then begin
        for i=1,ncomp do begin
           flux = componeline(instr,'[SII]6716',i,velsig=velsig)
           oplot,wave,yran[0]+flux/norm,color=75,linesty=2
           flux = componeline(instr,'[SII]6731',i,velsig=velsig)
           oplot,wave,yran[0]+flux/norm,color=75,linesty=2
        endfor
     endif
     loadct,0,/silent
     xyouts,xran6[0]+(xran6[1]-xran6[0])*0.05d,$
            yran[0]+(yran[1]-yran[0])*0.85d,$
            lab[5],charsize=1.5,charthick=2
;     multiplot,/doyaxis,/doxaxis
     !P.position = [0.7,0.1,0.95,0.3]
     ydat = specstars
     ymod = modstars
     yran = [min([ydat[i6],ymod[i6]]),max([ydat[i6],ymod[i6]])]
     plot,wave,ydat,xran=xran6,yran=yran,/xsty,/ysty,$
          pos=!P.position,/noerase
     loadct,13,/silent
     oplot,wave,ymod,color=255
  endif
;;   endif else begin
;;      multiplot
;;      multiplot
;;   endelse

;;   loadct,0,/silent
;;   multiplot,/reset

  loadct,0,/silent
  tit = 'EMISSION LINE FITS'
  xtit = 'Observed Wavelength (!3' + STRING(197B) + '!X)'
  ytit = textoidl('Normalized F_\lambda')
  xyouts,0.5,0.97,tit,/norm,align=0.5,charsize=2,charthick=2
  xyouts,0.5,0.02,xtit,/norm,align=0.5,charsize=2,charthick=2
  xyouts,0.03,0.5,ytit,/norm,align=0.5,charsize=2,orient=90,charthick=2

  tmpfile = outfile
  if (dops) then device,/close_file $
  else img = tvread(filename=tmpfile,/jpeg,/nodialog,quality=100)

  !X.tickinterval=defaultXtickint
  !X.minor=defaultXminor
  
end
