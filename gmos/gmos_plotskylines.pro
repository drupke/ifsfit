;
; History
;  10jul17  DSNR  re-written for new GMOS data, 
;                 based on new PLOTSTRONGLINES routine
;
; Plot strong emission line fits.
;

pro gmos_plotskylines,instr,outfile,ps=ps,zbuf=zbuf,comp=comp,$
                      plotilines=plotilines,oplothan2=oplothan2,$
                      xran=xran,yran=yran

  if ~ keyword_set(oplothan2) then begin
     dops=0
     dozbuf=0
     if keyword_set(ps) then dops=1
     if keyword_set(zbuf) then dozbuf=1
     if ~ keyword_set(comp) then comp=1

     if (dops) then begin
        set_plot,'ps',/copy,/interpolate
        device,filename=outfile+'.eps',/encapsulated,/inches,$
               xsize=7.5,ysize=7.5,bits_per_pixel=8,/color
        !P.charsize=1
        !P.charthick=1
     endif else if (dozbuf) then begin
        set_plot,'Z'
        device,decomposed=0,set_resolution=[640,640],set_pixel_depth=24
        !P.charsize=1
        !P.charthick=1
        erase
     endif else begin
        set_plot,'X'
        device,decomposed=0
        window,xsize=640,ysize=640,xpos=0,ypos=0,retain=2
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

  spectot = instr.spec
  specstars = instr.spec - instr.specfit
  speclines = instr.spec_nocnt
  specerr = instr.spec_err
  modtot = instr.specfit + (instr.spec - instr.spec_nocnt)
  modstars = instr.spec - instr.spec_nocnt
  modlines = instr.specfit

  norm = max(modstars)
  spectot /= norm
  specstars /= norm
  speclines /= norm
  specerr /= norm
  modtot /= norm
  modstars /= norm
  modlines /= norm
     
  if keyword_set(comp) then zbase = instr.z.gas[comp-1] $
  else zbase = instr.z.star

  lab = textoidl(['Na D','[OI]'])
  off = [-25d,25d]
  xran1 = ([5889.95d,5895.92d]+off) * (1d + zbase)
  xran2 = ([6300.3d,6363.7d]+off) * (1d + zbase)
  i1 = where(wave gt xran1[0] AND wave lt xran1[1],ct1)
  i2 = where(wave gt xran2[0] AND wave lt xran2[1],ct2)

; Na D

  loadct,0,/silent
  ydat = spectot
  ymod = modtot
  yran = [min([ydat[i1],ymod[i1]]),max([ydat[i1],ymod[i1]])]
  plot,wave,ydat,xran=xran1,yran=yran,/xsty,/ysty,$
       position=[0.1,0.55,0.9,0.95]
  loadct,13,/silent
  oplot,wave,specerr
  oplot,wave,ymod,color=75,thick=4
  oplot,wave,modstars,color=255,thick=5
  loadct,0,/silent
  xyouts,xran1[0]+(xran1[1]-xran1[0])*0.05d,$
         yran[0]+(yran[1]-yran[0])*0.85d,$
         lab[0],charsize=1.5,charthick=2

; [OI]
  
  ydat = spectot
  ymod = modtot
  yran = [min([ydat[i2],ymod[i2]]),max([ydat[i2],ymod[i2]])]
  plot,wave,ydat,xran=xran2,yran=yran,/xsty,/ysty,/noerase,$
       position=[0.1,0.05,0.9,0.45]
  loadct,13,/silent
  oplot,wave,specerr
  oplot,wave,ymod,color=75,thick=4
  oplot,wave,modstars,color=255,thick=4
  loadct,0,/silent
  xyouts,xran2[0]+(xran2[1]-xran2[0])*0.05d,$
         yran[0]+(yran[1]-yran[0])*0.85d,$
         lab[1],charsize=1.5,charthick=2

; Finish

  tmpfile = outfile
  if (dops) then device,/close_file $
  else img = tvread(filename=tmpfile,/jpeg,/nodialog,quality=100)
  
  !X.tickinterval=defaultXtickint
  !X.minor=defaultXminor

end
