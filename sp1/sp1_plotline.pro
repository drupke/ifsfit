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

pro sp1_plotline,initfile,linelab,ps=ps,comp=comp,woff=woff,cfoc=cfoc

  specresthresh = 10d

  specfile=''
  errfile=''
  outfile=''
  fitranstr=''
  mcompstr=''
  vacair=''
  sigstr=''
  dispstr=''
  strongstr=''

  openu,ilun,initfile,/get_lun
  readf,ilun,specfile,format='(A0)'
  readf,ilun,errfile,format='(A0)'
  readf,ilun,outfile,format='(A0)'
  readf,ilun,fitranstr,format='(A0)'
  readf,ilun,zguess,format='(D)'
  readf,ilun,mcompstr,format='(A0)'
  readf,ilun,dzstel,format='(D)'
  readf,ilun,vdisp,format='(D)'
  readf,ilun,specres,format='(D)'
  readf,ilun,vacair,format='(A0)'
  readf,ilun,sigstr,format='(A0)'
  readf,ilun,dispstr,format='(A0)'
  readf,ilun,strongstr,format='(A0)'
  free_lun,ilun

  infile = outfile+'.xdr'
  outfile += '_'+linelab
  if keyword_set(cfoc) then outfile += '_cfoc'+string(cfoc,format='(I0)')
  if specres gt specresthresh then velsig=1 else velsig=0

  restore,file=infile
  ncomp = struct.param[1]
  linepars = sepfitpars(struct.param,struct.perror)
  nlines = n_elements(linepars.flux)

  dops=0
  if keyword_set(ps) then dops=1
  if ~ keyword_set(comp) then comp=1

  if (dops) then begin
     set_plot,'ps',/copy,/interpolate
     device,filename=outfile+'.eps',/encapsulated,/inches,$
            xsize=10,ysize=7.5,bits_per_pixel=8,/color
     !P.charsize=1
     !P.charthick=2
     !P.thick=2
  endif else begin
     set_plot,'Z'
     device,decomposed=0,set_resolution=[1280,960],set_pixel_depth=24
     !P.charsize=1
     !P.charthick=1
     erase
  endelse

  defaultXtickint=!X.tickinterval
  defaultXminor=!X.minor

  wave = struct.wave
  spectot = struct.spec
  specstars = struct.spec - struct.specfit
  speclines = struct.spec_nocnt
  modtot = struct.specfit + (struct.spec - struct.spec_nocnt)
  modstars = struct.spec - struct.spec_nocnt
  modlines = struct.specfit

  norm = max(modstars)

  if keyword_set(cfoc) then begin
     spectot_cfoc = spectot
     modtot_cfoc = modtot
     for i=1,ncomp do begin
        if i ne cfoc then begin
           if linelab eq 'han2' then begin
              flux = componeline(struct,'[NII]6548',i,velsig=velsig) + $
                     componeline(struct,'Halpha',i,velsig=velsig) + $
                     componeline(struct,'[NII]6583',i,velsig=velsig)
           endif else if linelab eq 's2' then begin
              flux = componeline(struct,'[SII]6716',i,velsig=velsig) + $
                     componeline(struct,'[SII]6731',i,velsig=velsig)
           endif
           spectot_cfoc -= flux
           modtot_cfoc -= flux
        endif
     endfor
     spectot_cfoc /= norm
     modtot_cfoc /= norm
  endif

  spectot /= norm
  specstars /= norm
  speclines /= norm
  modtot /= norm
  modstars /= norm
  modlines /= norm

  if keyword_set(comp) then zbase = struct.z.gas[comp-1] $
  else zbase = struct.z.star

  if ~keyword_set(woff) then woff = [-50d,50d]
  if linelab eq 'han2' then begin
     xran = ([6548.0d,6583.4d]+woff) * (1d + zbase)
  endif else if linelab eq 's2' then begin
     xran = ([6716.3d,6730.7d]+woff) * (1d + zbase)
  endif
  ind = where(wave gt xran[0] AND wave lt xran[1],ct)

  xtit = 'Observed Wavelength (!3' + STRING(197B) + '!X)'
  ytit = textoidl('Normalized F_\lambda')

  loadct,0,/silent
  !X.tickinterval=10
  !X.minor=5

  if keyword_set(cfoc) then begin
     ydat = spectot_cfoc
     ymod = modtot_cfoc
  endif else begin
     ydat = spectot
     ymod = modtot
  endelse
  yran = [min([ydat[ind],ymod[ind]]),max([ydat[ind],ymod[ind]])]
  plot,wave,ydat,xran=xran,yran=yran,/xsty,/ysty,xtit=xtit,ytit=ytit
  loadct,13,/silent
  oplot,wave,ymod,color=255
  for i=1,ncomp do begin
     if linelab eq 'han2' then begin
        flux = componeline(struct,'Halpha',i,velsig=velsig)
        oplot,wave,yran[0]+flux/norm,color=75,linesty=2
        flux = componeline(struct,'[NII]6548',i,velsig=velsig)
        oplot,wave,yran[0]+flux/norm,color=75,linesty=2
        flux = componeline(struct,'[NII]6583',i,velsig=velsig)
        oplot,wave,yran[0]+flux/norm,color=75,linesty=2
     endif else if linelab eq 's2' then begin
        flux = componeline(struct,'[SII]6716',i,velsig=velsig)
        oplot,wave,yran[0]+flux/norm,color=75,linesty=2
        flux = componeline(struct,'[SII]6731',i,velsig=velsig)
        oplot,wave,yran[0]+flux/norm,color=75,linesty=2 
    endif
  endfor

  tmpfile = outfile
  if (dops) then device,/close_file $
  else img = tvread(filename=tmpfile,/jpeg,/nodialog,quality=100)

  !X.tickinterval=defaultXtickint
  !X.minor=defaultXminor
  
end
