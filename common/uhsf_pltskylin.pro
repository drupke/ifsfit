;
; History
;  10jul17  DSNR  re-written for new GMOS data, 
;                 based on new PLOTSTRONGLINES routine
;  13sep10  DSNR  added ability to plot [NI] line
;
; Plot strong emission line fits.
;

pro uhsf_pltskylin,instr,outfile,ps=ps,zbuf=zbuf,comp=comp,$
                   plotilines=plotilines,xran=xran,yran=yran,$
                   outlines=outlines,argslinelist=argslinelist

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

; Get linelist
  if keyword_set(argslinelist) then $
     linelist = call_function('gmos_initlinelist',_extra=argslinelist) $
  else linelist = gmos_initlinelist()


  off = [-100d,100d]
  for i=0,1 do begin

     lab = outlines[i]
     iwave = where(linelist.label eq outlines[i])
     tmpwave = linelist.wave[iwave]
     xran = (tmpwave[0]+off) * (1d + zbase)
     ind = where(wave gt xran AND wave lt xran,ct)

     ydat = spectot
     ymod = modtot
     yran = [min([ydat[ind],ymod[ind]]),max([ydat[ind],ymod[ind]])]
     cgplot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty,layout=[1,2,i+1],$
            color='White',axiscolor='White'
     cgoplot,wave,specerr,color='White'
     cgoplot,wave,ymod,color='Blue',thick=4
     cgoplot,wave,modstars,color='Red',thick=5
     cgtext,xran[0]+(xran[1]-xran[0])*0.05d,$
            yran[0]+(yran[1]-yran[0])*0.8d,$
            lab,charsize=1.5,charthick=2

     if i eq 0 AND n_elements(outlines) eq 1 then goto,finish

  endfor

; Finish

finish:

  tmpfile = outfile
  if (dops) then device,/close_file $
  else img = cgsnapshot(filename=tmpfile,/jpeg,/nodialog,quality=100)

end
