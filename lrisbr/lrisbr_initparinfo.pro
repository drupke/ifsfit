function lrisbr_initparinfo,linename,linelambda,linelambdaz,initflux,$
                            initsig,z,tielambda=tielambda,zfix=zfix,$
                            bluered=bluered,floatcomp=floatcomp
;
; 10mar18  DSNR  created
;

; Number of velocity components
  ncomp = n_elements(z.gas)
  zbase = z.gas
  ztilt = z.redtilt
  zint = z.redint
  if keyword_set(bluered) then begin
     if bluered eq 'blue' then begin
        ncomp = n_elements(z.blue)
        zbase = z.blue
        ztilt = z.bluetilt
        zint = z.blueint
     endif
  endif

; Number of initial parameters before Gaussian parameters begin
  ppoff0 = 9
  ppoff = ppoff0 + ncomp-1

; Number of emission lines to fit
  nline = N_ELEMENTS(linename)

  parinfo = REPLICATE({value:0d, fixed:0, limited:[0B,0B], tied:'', $
                       limits:[0D,0D], step:0d, mpprint:0b}, $
                      ppoff+ncomp*(nline*3))

  itie = WHERE(linename eq tielambda)

; Number of initial parameters before Gaussian parameters begin
  parinfo[0].value = ppoff
  parinfo[0].fixed = 1B

; Number of velocity components
  parinfo[1].value = ncomp
  parinfo[1].fixed = 1B  

; redshift  
  parinfo[2].value   = zbase[0]
  parinfo[2].limited = [1B,1B]
  parinfo[2].limits  = [zbase[0]-0.003,zbase[0]+0.003]
  parinfo[2].step = 0.0001d
  parinfo[2].mpprint = 1b

; redshift tilt  
  parinfo[3].value   = ztilt
  parinfo[3].limited = [1B,1B]
  parinfo[3].limits  = [-5d-7,5d-7]
  parinfo[3].mpprint = 1b
  if keyword_set(zfix) then parinfo[3].fixed = 1B

; [OII] ratio
  parinfo[4].value   = 1.3d
  parinfo[4].limited = [1B,1B]
  parinfo[4].limits  = [0.27d,1.48d]
  
; [SII] ratio
  parinfo[5].value   = 1.3d
  parinfo[5].limited = [1B,1B]
  parinfo[5].limits  = [0.44d,1.43d]

; Floating broad component
  if keyword_set(floatcomp) then begin
;    Flux estimate
     parinfo[6].value = floatcomp[0]
     parinfo[6].limited = [1B,1B]
     parinfo[6].limits = floatcomp[0]*[0.3d,3d]
     parinfo[6].mpprint = 1b
;     parinfo[14].relstep = 0.01d
;    Location estimate, in A
     parinfo[7].value = floatcomp[1]
     parinfo[7].limited = [1B,1B]
     parinfo[7].limits = floatcomp[1]-[10d,-10d]
     parinfo[7].mpprint = 1b
;    Width estimate
     parinfo[8].value = floatcomp[2]
     parinfo[8].limited = [1B,1B]
     parinfo[8].limits = floatcomp[2]*[0.6d,1.4d]
     parinfo[8].mpprint = 1b
;     parinfo[16].relstep = 0.01d
  endif else begin
     parinfo[6].value = 0d
     parinfo[6].fixed = 1B
     parinfo[7].value = 6563d
     parinfo[7].fixed = 1B
     parinfo[8].value = 1d
     parinfo[8].fixed = 1B
  endelse
  
; main velocity component
  for iline=0, nline-1 do begin

     ifoff = ppoff+iline*3
     iwoff = ifoff+1
     isoff = ifoff+2
      
     parinfo[ifoff].value = initflux[iline]
     parinfo[ifoff].limited[0] = 1B
     parinfo[ifoff].limits[0]  = 0d
    
     parinfo[iwoff].value = linelambdaz[iline,0]
     parinfo[iwoff].tied = string(linelambda[iline],format='(D0.2)') + $
                           '*(1d + P[2] + P[3]*(' + $
                           string(linelambda[iline],'-',zint,$
                                  format='(D0.2,A0,D0.2)')+'))'
;     parinfo[iwoff].tied = string(linelambda[iline],$
;                                  format='(D0.2)')+'*(1d + P[2])'
    
     parinfo[isoff].value = initsig[iline]
     parinfo[isoff].limited = [1B,1B]
     parinfo[isoff].limits = [1d, 15d]
     parinfo[isoff].tied = STRING("P[",ppoff+2+itie*3,"]",format='(A,I0,A)')
     if (linename[iline] eq tielambda) then parinfo[isoff].tied = ''
     
  endfor

  linea = where(linename eq '[OII]3729')
  lineb = where(linename eq '[OII]3726')
  parinfo[ppoff+lineb*3].tied = 'P['+string(ppoff+linea*3,format='(I0)')+$
                                ']/P[4]'
  linea = where(linename eq '[NeIII]3967',cta)
  lineb = where(linename eq '[NeIII]3869',ctb)
  if (cta gt 0 AND ctb gt 0) then $
     parinfo[ppoff+linea*3].tied = 'P['+string(ppoff+lineb*3,format='(I0)')+$
                                   ']/3.1d'
  linea = where(linename eq '[OIII]4959')
  lineb = where(linename eq '[OIII]5007')
  parinfo[ppoff+linea*3].tied = 'P['+string(ppoff+lineb*3,format='(I0)')+']/3.0d'
  linea = where(linename eq '[OI]6364')
  lineb = where(linename eq '[OI]6300')
  parinfo[ppoff+linea*3].tied = 'P['+string(ppoff+lineb*3,format='(I0)')+']/3.0d'
  linea = where(linename eq '[NII]6548')
  lineb = where(linename eq '[NII]6583')
  parinfo[ppoff+linea*3].tied = 'P['+string(ppoff+lineb*3,format='(I0)')+']/3.0d'
  linea = where(linename eq '[SII]6716')
  lineb = where(linename eq '[SII]6731')
  parinfo[ppoff+lineb*3].tied = 'P['+string(ppoff+linea*3,format='(I0)')+$
                                ']/P[5]'

; other velocity components
  for i=1,ncomp-1 do begin

;    redshift
     parinfo[ppoff0-1+i].value = zbase[i]
     parinfo[ppoff0-1+i].limited = [1B,1B]
     parinfo[ppoff0-1+i].limits = [zbase[i]-0.003,zbase[i]+0.003]
     parinfo[ppoff0-1+i].step = 0.0001d
     parinfo[ppoff0-1+i].mpprint = 1b

     foff = ppoff+i*nline*3
     woff = foff+1
     soff = foff+2

;    line parameters
     for iline=0, nline-1 do begin

        ifoff = foff + iline*3
        iwoff = woff + iline*3
        isoff = soff + iline*3

        parinfo[ifoff].value = initflux[iline,i]
        parinfo[ifoff].limited[0] = 1B
        parinfo[ifoff].limits[0]  = 0d
    
        parinfo[iwoff].value = linelambdaz[iline,i]
        parinfo[iwoff].tied = string(linelambda[iline],$
                                     '*(1d + P[',$
                                     ppoff0-1+i,$
                                     ']+P[3]*(',$
                                     linelambda[iline],'-',zint,$
                                     '))',$
                                     format='(D0.2,A0,I0,A0,D0.2,A0,D0.2,A0)')
;        parinfo[iwoff].tied = $
;           string(linelambda[iline],'*(1d + P[',ppoff0-1+i,'])',$
;                  format='(D0.2,A0,I0,A0)')
    
        parinfo[isoff].value = initsig[iline,i]
        parinfo[isoff].limited = [1B,1B]
        parinfo[isoff].limits = [1d,10d]
        parinfo[isoff].tied = $
           STRING("P[",soff+itie*3,"]",format='(A0,I0,A0)')
        if (linename[iline] eq tielambda) then parinfo[isoff].tied = ''

     endfor

     linea = where(linename eq '[OII]3729')
     lineb = where(linename eq '[OII]3726')
     parinfo[foff+lineb*3].tied = 'P['+string(foff+linea*3,$
                                              format='(I0)')+']/P[4]'
     linea = where(linename eq '[NeIII]3967',cta)
     lineb = where(linename eq '[NeIII]3869',ctb)
     if (cta gt 0 AND ctb gt 0) then $
        parinfo[foff+linea*3].tied = 'P['+string(foff+lineb*3,format='(I0)')+$
                                     ']/3.1d'
     linea = where(linename eq '[OIII]4959')
     lineb = where(linename eq '[OIII]5007')
     parinfo[foff+linea*3].tied = 'P['+string(foff+lineb*3,format='(I0)')+']/3.0d'
     linea = where(linename eq '[OI]6364')
     lineb = where(linename eq '[OI]6300')
     parinfo[foff+linea*3].tied = 'P['+string(foff+lineb*3,format='(I0)')+']/3.0d'
     linea = where(linename eq '[NII]6548')
     lineb = where(linename eq '[NII]6583')
     parinfo[foff+linea*3].tied = 'P['+string(foff+lineb*3,format='(I0)')+']/3.0d'
     linea = where(linename eq '[SII]6716')
     lineb = where(linename eq '[SII]6731')
     parinfo[foff+lineb*3].tied = 'P['+string(foff+linea*3,$
                                              format='(I0)')+']/P[5]'
     
  endfor

  return, parinfo

end
