function sp1_initparinfo,linename,linelambda,linelambdaz,$
                         initflux,initsig,z,zerolines=zerolines,$
                         sigma=sigma,specres=specres,$
                         specrthresh=specrthresh
;
; 09nov24  DSNR  created
; 11mar07  DSNR  changed siglim from 2000d to 600d
;

; Sigma limit.  10 A in sigma is ~1000 km/s FWHM at 6563 A.
  if specres le specrthresh then siglim = [specres,20d] $
;  else siglim = [specres,2000d]
  else siglim = [specres,600d]
  if ~ keyword_set(sigma) then sigma = ''

; Number of velocity components
  ncomp = n_elements(z.gas)
  zbase = z.gas

; Number of initial parameters before Gaussian parameters begin
  ppoff0 = 5
  ppoff = ppoff0 + ncomp-1

; Number of emission lines to fit
  nline = N_ELEMENTS(linename)

  parinfo = REPLICATE({value:0d, fixed:0, limited:[0,0], tied:'', $
                       limits:[0.D,0], step:0d, mpprint:0b}, $
                      ppoff+ncomp*(nline*3))

  halpha = WHERE(linename eq 'Halpha')

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

; [OII] ratio
  parinfo[3].value   = 1.3d
  parinfo[3].limited = [1B,1B]
  parinfo[3].limits  = [0.27d,1.48d]

; [SII] ratio
  parinfo[4].value   = 1.3d
  parinfo[4].limited = [1B,1B]
  parinfo[4].limits  = [0.44d,1.43d]
  
; main velocity component
  for iline=0, nline-1 do begin

     ifoff = ppoff+iline*3
     iwoff = ifoff+1
     isoff = ifoff+2
      
     parinfo[ifoff].value = initflux[iline]
     parinfo[ifoff].limited[0] = 1B
     parinfo[ifoff].limits[0]  = [0d] ;;, 2d*MAX(flux)]
     if keyword_set(zerolines) then begin
        if zerolines[iline,0] eq 1 then begin
           parinfo[ifoff].value = 0d
           parinfo[ifoff].fixed = 1B
        endif
     endif
    
     parinfo[iwoff].value = linelambdaz[iline,0]
     parinfo[iwoff].tied = string(linelambda[iline],$
                                  format='(D0.2)')+$
                           '*(1d +P[2])'

     parinfo[isoff].value = initsig[iline,0]
     parinfo[isoff].limited = [1B,1B]
     parinfo[isoff].limits = siglim
     if sigma eq 'fix' then begin
        parinfo[isoff].tied = STRING('P[',ppoff+2+halpha*3,']',$
                                     format='(A,I0,A)')
        if (linename[iline] eq 'Halpha') then begin
           parinfo[isoff].tied = ''
           parinfo[isoff].mpprint = 1b
        endif
     endif
    
  endfor

  linea = where(linename eq '[OII]3729')
  lineb = where(linename eq '[OII]3726')
  parinfo[ppoff+lineb*3].tied = 'P['+string(ppoff+linea*3,$
                                            format='(I0)')+']/P[3]'
  if sigma ne 'fix' then $
     parinfo[ppoff+lineb*3+2].tied = 'P['+string(ppoff+linea*3+2,$
                                                 format='(I0)')+']'
     
  linea = where(linename eq '[NeIII]3967',cta)
  lineb = where(linename eq '[NeIII]3869',ctb)
  if (cta gt 0 AND ctb gt 0) then begin
     parinfo[ppoff+linea*3].tied = 'P['+string(ppoff+lineb*3,$
                                               format='(I0)')+']/3.1d'
     if sigma ne 'fix' then $
        parinfo[ppoff+linea*3+2].tied = 'P['+string(ppoff+lineb*3+2,$
                                                    format='(I0)')+']'
  endif

  linea = where(linename eq '[OIII]4959')
  lineb = where(linename eq '[OIII]5007')
  parinfo[ppoff+linea*3].tied = 'P['+string(ppoff+lineb*3,$
                                            format='(I0)')+']/3.0d'
  if sigma ne 'fix' then $
     parinfo[ppoff+linea*3+2].tied = 'P['+string(ppoff+lineb*3+2,$
                                                 format='(I0)')+']'
  linea = where(linename eq '[OI]6364')
  lineb = where(linename eq '[OI]6300')
  parinfo[ppoff+linea*3].tied = 'P['+string(ppoff+lineb*3,$
                                            format='(I0)')+']/3.0d'
  if sigma ne 'fix' then $
     parinfo[ppoff+linea*3+2].tied = 'P['+string(ppoff+lineb*3+2,$
                                                 format='(I0)')+']'
  linea = where(linename eq '[NII]6548')
  lineb = where(linename eq '[NII]6583')
  parinfo[ppoff+linea*3].tied = 'P['+string(ppoff+lineb*3,$
                                            format='(I0)')+']/3.0d'
  if sigma ne 'fix' then $
     parinfo[ppoff+linea*3+2].tied = 'P['+string(ppoff+lineb*3+2,$
                                                 format='(I0)')+']'
  linea = where(linename eq '[SII]6716')
  lineb = where(linename eq '[SII]6731')
  parinfo[ppoff+lineb*3].tied = 'P['+string(ppoff+linea*3,$
                                            format='(I0)')+']/P[4]'
  if sigma ne 'fix' then $
     parinfo[ppoff+lineb*3+2].tied = 'P['+string(ppoff+linea*3+2,$
                                                 format='(I0)')+']'

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
        parinfo[ifoff].limits[0]  = [0d]
        if keyword_set(zerolines) then begin
           if zerolines[iline,i] eq 1 then begin
              parinfo[ifoff].value = 0d
              parinfo[ifoff].fixed = 1B
           endif
        endif
    
        parinfo[iwoff].value = linelambdaz[iline,i]
        parinfo[iwoff].tied = string(linelambda[iline],$
                                     '*(1d + P[',ppoff0-1+i,'])',$
                                     format='(D0.2,A0,I0,A0)')

        parinfo[isoff].value = initsig[iline,i]
        parinfo[isoff].limited = [1B,1B]
        parinfo[isoff].limits = siglim
        if sigma eq 'fix' then begin
           parinfo[isoff].tied = STRING('P[',soff+halpha*3,']',$
                                        format='(A,I0,A)')
           if (linename[iline] eq 'Halpha') then begin
              parinfo[isoff].tied = ''
              parinfo[isoff].mpprint = 1b
           endif
        endif
        
     endfor

     linea = where(linename eq '[OII]3729')
     lineb = where(linename eq '[OII]3726')
     parinfo[foff+lineb*3].tied = 'P['+string(foff+linea*3,$
                                               format='(I0)')+']/P[3]'
     if sigma ne 'fix' then $
        parinfo[foff+lineb*3+2].tied = 'P['+string(foff+linea*3+2,$
                                                    format='(I0)')+']'
     
     linea = where(linename eq '[NeIII]3967',cta)
     lineb = where(linename eq '[NeIII]3869',ctb)
     if (cta gt 0 AND ctb gt 0) then begin
        parinfo[foff+linea*3].tied = 'P['+string(foff+lineb*3,$
                                                  format='(I0)')+']/3.1d'
        if sigma ne 'fix' then $
           parinfo[foff+linea*3+2].tied = 'P['+string(foff+lineb*3+2,$
                                                       format='(I0)')+']'
     endif

     linea = where(linename eq '[OIII]4959')
     lineb = where(linename eq '[OIII]5007')
     parinfo[foff+linea*3].tied = 'P['+string(foff+lineb*3,$
                                               format='(I0)')+']/3.0d'
     if sigma ne 'fix' then $
        parinfo[foff+linea*3+2].tied = 'P['+string(foff+lineb*3+2,$
                                                    format='(I0)')+']'
     linea = where(linename eq '[OI]6364')
     lineb = where(linename eq '[OI]6300')
     parinfo[foff+linea*3].tied = 'P['+string(foff+lineb*3,$
                                               format='(I0)')+']/3.0d'
     if sigma ne 'fix' then $
        parinfo[foff+linea*3+2].tied = 'P['+string(foff+lineb*3+2,$
                                                    format='(I0)')+']'
     linea = where(linename eq '[NII]6548')
     lineb = where(linename eq '[NII]6583')
     parinfo[foff+linea*3].tied = 'P['+string(foff+lineb*3,$
                                               format='(I0)')+']/3.0d'
     if sigma ne 'fix' then $
        parinfo[foff+linea*3+2].tied = 'P['+string(foff+lineb*3+2,$
                                                    format='(I0)')+']'
     linea = where(linename eq '[SII]6716')
     lineb = where(linename eq '[SII]6731')
     parinfo[foff+lineb*3].tied = 'P['+string(foff+linea*3,$
                                               format='(I0)')+']/P[4]'
     if sigma ne 'fix' then $
        parinfo[foff+lineb*3+2].tied = 'P['+string(foff+linea*3+2,$
                                                    format='(I0)')+']'


  endfor

  return, parinfo

end
