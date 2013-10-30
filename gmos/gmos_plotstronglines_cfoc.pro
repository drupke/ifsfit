;
; History
;  10jul17  DSNR  re-written for new GMOS data, 
;                 based on new PLOTSTRONGLINES routine
;
; Plot strong emission line fits.
;

function componeline,instr,line,comp,center=center
  
  iline = where(instr.linelabel eq line,ct)
  ppoff = instr.param[0]
  ncomp = instr.param[1]
  ppoff0 = ppoff - (ncomp-1)

  nline = n_elements(instr.linelabel)
  indices = ppoff+(comp-1)*nline*3+iline*3
  indices = indices[0] + indgen(3)
;  print,comp,indices[0],instr.param[indices]
  flux = gaussian(instr.wave,instr.param[indices],/double)

; central wavelength
  tmp = instr.param[indices]
  center = tmp[1]

  return,flux

end

pro gmos_plotstronglines_cfoc,instr,outfile,ps=ps,zbuf=zbuf,comp=comp,$
                         plotilines=plotilines,oplothan2=oplothan2,$
                         xran=xran,yran=yran,cfoc=cfoc,woff=woff

  if ~ keyword_set(oplothan2) then begin
     dops=0
     dozbuf=0
     if keyword_set(ps) then dops=1
     if keyword_set(zbuf) then dozbuf=1
     if ~ keyword_set(comp) then comp=1

     if (dops) then begin
        set_plot,'ps',/copy,/interpolate
        device,filename=outfile+'.eps',/encapsulated,/inches,$
               xsize=10,ysize=7.5,bits_per_pixel=8,/color
        !P.charsize=1
        !P.charthick=1
     endif else if (dozbuf) then begin
        set_plot,'Z'
        device,decomposed=0,set_resolution=[1920,960],set_pixel_depth=24
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
     !X.minor=5

  endif

  if keyword_set(plotilines) then begin
     ncomp = instr.param[1]
     colors = [255,75,125]
  endif

  wave = instr.wave

; polynomial near Ha/[NII]
  ypoly = dblarr(n_elements(wave))
  npoly = instr.param[4]
  wavelo = instr.param[5]
  wavehi = instr.param[6]
  n2hawave = where(wave ge wavelo AND wave le wavehi)
  ypoly[n2hawave] = poly(wave[n2hawave]-mean(wave[n2hawave]),instr.param[7:7+npoly])

  spectot = instr.spec
  specstars = instr.spec - instr.specfit + ypoly
  speclines = instr.spec_nocnt - ypoly
  specerr = instr.spec_err
  modtot = instr.specfit + (instr.spec - instr.spec_nocnt)
  modstars = instr.spec - instr.spec_nocnt + ypoly
  modlines = instr.specfit - ypoly

  norm = max(modstars)

  if keyword_set(cfoc) then begin
     speclines_cfoc = speclines
     modlines_cfoc = modlines
     for i=1,ncomp do begin
        if i ne cfoc then begin
           flux = componeline(instr,'[NII]6548',i) + $
                  componeline(instr,'Halpha',i) + $
                  componeline(instr,'[NII]6583',i) + $
                  componeline(instr,'[OI]6300',i) + $
                  componeline(instr,'[OI]6364',i)
           speclines_cfoc -= flux
           modlines_cfoc -= flux
        endif
     endfor
     speclines_cfoc /= norm
     modlines_cfoc /= norm
  endif

  spectot /= norm
  specstars /= norm
  speclines /= norm
  specerr /= norm
  modtot /= norm
  modstars /= norm
  modlines /= norm
     
  if keyword_set(comp) then zbase = instr.z.gas[comp-1] $
  else zbase = instr.z.star

  lab = textoidl(['H\alpha/[NII]','[OI]'])
  off = [-15d,15d]
  xran1 = ([6548.0d,6583.4d]+6d*off) * (1d + zbase)
  xran2 = ([6300.3d,6363.7d]+6d*off) * (1d + zbase)
  if ~ keyword_set(woff) then woff = 3d*off
  xran1ss = ([6548.0d,6583.4d]+woff) * (1d + zbase)
  xran2ss = ([6300.3d,6363.7d]+woff) * (1d + zbase)
  telran1 = [6865d,6925d]
  telran2 = [6275d,6285d]
  i1 = where(wave gt xran1[0] AND wave lt xran1[1],ct1)
  i2 = where(wave gt xran2[0] AND wave lt xran2[1],ct2)
  i1ss = where(wave gt xran1ss[0] AND wave lt xran1ss[1],ct1ss)
  i2ss = where(wave gt xran2ss[0] AND wave lt xran2ss[1],ct2ss)

; Halpha / [NII]

  if ~ keyword_set(oplothan2) then begin

     pxl = [0,0.06,0.37]
     pxu = [0,0.33,0.95]
     pyl = [0,0.05,0.15,0.55,0.65]
     pyu = [0,0.15,0.45,0.65,0.95]

     loadct,0,/silent
     ydat = spectot
     ymod = modtot
     yran = [min([ydat[i1],ymod[i1]]),max([ydat[i1],ymod[i1]])]
     ip = [1,4]
     plot,wave,ydat,xran=xran1,yran=yran,/xsty,/ysty,$
          position=[pxl[ip[0]],pyl[ip[1]],pxu[ip[0]],pyu[ip[1]]],$
          xtickn=replicate(' ',60)
     loadct,13,/silent
     polyfill,[telran1,reverse(telran1)],[yran[0],yran[0],yran[1],yran[1]],$
              /line_fill,spacing=0.5d,color=125,noclip=0,$
              clip=[xran1[0],yran[0],xran1[1],yran[1]]
     polyfill,[telran2,reverse(telran2)],[yran[0],yran[0],yran[1],yran[1]],$
              /line_fill,spacing=0.5d,color=125,noclip=0,$
              clip=[xran1[0],yran[0],xran1[1],yran[1]]
     oplot,wave,ymod,color=75,thick=4
     oplot,wave,modstars,color=255,thick=5
     loadct,0,/silent
     xyouts,xran1[0]+(xran1[1]-xran1[0])*0.05d,$
            yran[0]+(yran[1]-yran[0])*0.85d,$
            lab[0],charsize=1.5,charthick=2
     ydat = specerr
     yran = [min(ydat[i1]),max(ydat[i1])]
     ip = [1,3]
     plot,wave,ydat,xran=xran1,yran=yran,/xsty,/ysty,/noerase,$
          position=[pxl[ip[0]],pyl[ip[1]],pxu[ip[0]],pyu[ip[1]]]

     if keyword_set(cfoc) then begin
        ydat = speclines_cfoc
        ymod = modlines_cfoc
     endif else begin
        ydat = speclines
        ymod = modlines
     endelse
     yran = [min([ydat[i1ss],ymod[i1ss]]),max([ydat[i1ss],ymod[i1ss]])]
     ip = [2,4]
     plot,wave,ydat,xran=xran1ss,yran=yran,/xsty,/ysty,/noerase,$
          position=[pxl[ip[0]],pyl[ip[1]],pxu[ip[0]],pyu[ip[1]]],$
          xtickn=replicate(' ',60)
     loadct,13,/silent
     ; polyfill,[telran1,reverse(telran1)],[yran[0],yran[0],yran[1],yran[1]],$
     ;          /line_fill,spacing=0.5d,color=125,noclip=0,$
     ;          clip=[xran1ss[0],yran[0],xran1ss[1],yran[1]]
     ; polyfill,[telran2,reverse(telran2)],[yran[0],yran[0],yran[1],yran[1]],$
     ;          /line_fill,spacing=0.5d,color=125,noclip=0,$
     ;          clip=[xran1ss[0],yran[0],xran1ss[1],yran[1]]
     if keyword_set(plotilines) then begin
        for i=1,ncomp do begin
           flux = componeline(instr,'[NII]6548',i,center=clam)
           oplot,wave,flux/norm,color=colors[i-1],thick=4
           oplot,[clam,clam],yran,linesty=2,color=colors[i-1],thick=2
           flux = componeline(instr,'Halpha',i,center=clam)
           oplot,wave,flux/norm,color=colors[i-1],thick=4
           oplot,[clam,clam],yran,linesty=2,color=colors[i-1],thick=2
           flux = componeline(instr,'[NII]6583',i,center=clam)
           oplot,wave,flux/norm,color=colors[i-1],thick=4
           oplot,[clam,clam],yran,linesty=2,color=colors[i-1],thick=2
        endfor
     endif
     loadct,0,/silent
     ydat = specerr
     yran = [min(ydat[i1ss]),max(ydat[i1ss])]
     ip = [2,3]
     plot,wave,ydat,xran=xran1ss,yran=yran,/xsty,/ysty,/noerase,$
          position=[pxl[ip[0]],pyl[ip[1]],pxu[ip[0]],pyu[ip[1]]]

; [OI]

     ydat = spectot
     ymod = modtot
     yran = [min([ydat[i2],ymod[i2]]),max([ydat[i2],ymod[i2]])]
     ip = [1,2]
     plot,wave,ydat,xran=xran2,yran=yran,/xsty,/ysty,/noerase,$
          position=[pxl[ip[0]],pyl[ip[1]],pxu[ip[0]],pyu[ip[1]]],$
          xtickn=replicate(' ',60)
     loadct,13,/silent
     polyfill,[telran1,reverse(telran1)],[yran[0],yran[0],yran[1],yran[1]],$
              /line_fill,spacing=0.5d,color=125,noclip=0,$
              clip=[xran2[0],yran[0],xran2[1],yran[1]]
     polyfill,[telran2,reverse(telran2)],[yran[0],yran[0],yran[1],yran[1]],$
              /line_fill,spacing=0.5d,color=125,noclip=0,$
              clip=[xran2[0],yran[0],xran2[1],yran[1]]
     oplot,wave,ymod,color=75,thick=4
     oplot,wave,modstars,color=255,thick=4
     loadct,0,/silent
     xyouts,xran2[0]+(xran2[1]-xran2[0])*0.05d,$
            yran[0]+(yran[1]-yran[0])*0.85d,$
            lab[1],charsize=1.5,charthick=2
     ydat = specerr
     yran = [min(ydat[i2]),max(ydat[i2])]
     ip = [1,1]
     plot,wave,ydat,xran=xran2,yran=yran,/xsty,/ysty,/noerase,$
          position=[pxl[ip[0]],pyl[ip[1]],pxu[ip[0]],pyu[ip[1]]]

     if keyword_set(cfoc) then begin
        ydat = speclines_cfoc
        ymod = modlines_cfoc
     endif else begin
        ydat = speclines
        ymod = modlines
     endelse
     yran = [min([ydat[i2],ymod[i2]]),max([ydat[i2],ymod[i2]])] 
     ip = [2,2]
     plot,wave,ydat,xran=xran2ss,yran=yran,/xsty,/ysty,/noerase,$
          position=[pxl[ip[0]],pyl[ip[1]],pxu[ip[0]],pyu[ip[1]]],$
          xtickn=replicate(' ',60)
     loadct,13,/silent
     polyfill,[telran1,reverse(telran1)],[yran[0],yran[0],yran[1],yran[1]],$
              /line_fill,spacing=0.5d,color=125,noclip=0,$
              clip=[xran2ss[0],yran[0],xran2ss[1],yran[1]]
     polyfill,[telran2,reverse(telran2)],[yran[0],yran[0],yran[1],yran[1]],$
              /line_fill,spacing=0.5d,color=125,noclip=0,$
              clip=[xran2ss[0],yran[0],xran2ss[1],yran[1]]
     if keyword_set(plotilines) then begin
        for i=1,ncomp do begin
           flux = componeline(instr,'[OI]6300',i,center=clam)
           oplot,wave,flux/norm,color=colors[i-1],thick=4
           oplot,[clam,clam],yran,linesty=2,color=colors[i-1],thick=2
           flux = componeline(instr,'[OI]6364',i,center=clam)
           oplot,wave,flux/norm,color=colors[i-1],thick=4
           oplot,[clam,clam],yran,linesty=2,color=colors[i-1],thick=2
        endfor
     endif
     loadct,0,/silent
     ydat = specerr
     yran = [min(ydat[i2ss]),max(ydat[i2ss])]
     ip = [2,1]
     plot,wave,ydat,xran=xran2ss,yran=yran,/xsty,/ysty,/noerase,$
          position=[pxl[ip[0]],pyl[ip[1]],pxu[ip[0]],pyu[ip[1]]]

; Finish

     tmpfile = outfile
     if (dops) then device,/close_file $
     else img = cgsnapshot(filename=tmpfile,/jpeg,/nodialog,quality=100)

     !X.tickinterval=defaultXtickint
     !X.minor=defaultXminor

  endif else begin

     if ~keyword_set(xran) then xran = xran1s

     ydat = spectot
     ymod = modtot
     yran = [min([ydat[i1],ymod[i1]]),max([ydat[i1],ymod[i1]])]
     ; ydat = speclines
     ; ymod = modlines
     ; yran = [min([ydat[i1s],ymod[i1s]]),max([ydat[i1s],ymod[i1s]])]

     plot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty
     loadct,13,/silent
     oplot,wave,ymod,color=255,thick=4
     oplot,wave,modstars,color=75,thick=4

  endelse

  
end
