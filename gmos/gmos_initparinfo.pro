function gmos_initparinfo,linename,linelambda,linelambdaz,initflux,$
                          initsig,z,cont_sig=cont_sig,siglim=siglim,$
                          zfix=zfix,n2hafix=n2hafix,sigfix=sigfix,$
                          polypars=polypars
;
; 09jun01  DSNR  created
; 09jun08  DSNR  added multiple components
; 10may27  DSNR  re-written to fit in observed frame
;

; Sigma limits. 0.7-0.8 A is approximately resolution for GMOS B600
; grating.
  if ~ keyword_set(siglim) then siglim=[0.7d,10d]

; Number of velocity components
  ncomp = n_elements(z.gas)

; Polynomial parameters for continuum fit tweaking
  if keyword_set(polypars) then begin
     npoly=polypars[0]
     wavelo=polypars[1]
     wavehi=polypars[2]
  endif else begin
     npoly = 0
     wavelo=1d
     wavehi=1d
  endelse

; Number of initial parameters before Gaussian parameters begin
  ppoff0 = 4+(npoly+4)
  ppoff = ppoff0 + ncomp-1

; Number of emission lines to fit
  nline = N_ELEMENTS(linename)

  parinfo = REPLICATE({value:0d, fixed:0b, limited:[0B,0B], tied:'', $
                       limits:[0d,0d], step:0d, mpprint:0b, mpside:2}, $
                      ppoff+ncomp*(nline*3))

  halpha = WHERE(linename eq 'Halpha')

; Number of initial parameters before Gaussian parameters begin
  parinfo[0].value = ppoff
  parinfo[0].fixed = 1B

; Number of velocity components
  parinfo[1].value = ncomp
  parinfo[1].fixed = 1B

; redshift
  parinfo[2].value   = z.gas[0]
  parinfo[2].limited = [1B,1B]
  parinfo[2].limits  = [z.gas[0]-0.003,z.gas[0]+0.003]
;  parinfo[2].step = 0.0001d
  parinfo[2].mpprint = 1b
  if keyword_set(zfix) then if zfix[0] then parinfo[2].fixed=1B

; [SII] ratio
  parinfo[3].value   = 1.3d
  parinfo[3].limited = [1B,1B]
  parinfo[3].limits  = [0.44d,1.43d]

; Underlying polynomial
  parinfo[4].value = npoly
  parinfo[4].fixed = 1b
  parinfo[5].value = wavelo
  parinfo[5].fixed = 1b
  parinfo[6].value = wavehi
  parinfo[6].fixed = 1b
  for i=0,npoly do parinfo[7+i].value = 0d

; main velocity component
  for iline=0, nline-1 do begin

     ifoff = ppoff+iline*3
     iwoff = ifoff+1
     isoff = ifoff+2
      
     parinfo[ifoff].value = initflux[iline]
     parinfo[ifoff].limited[0] = 1B
     parinfo[ifoff].limits[0]  = [0d] ;;, 2d*MAX(flux)]
    
     parinfo[iwoff].value = linelambdaz[iline,0]
     parinfo[iwoff].tied = string(linelambda[iline],format='(D0.2)')+$
                           '*(1d +P[2])'
    
     parinfo[isoff].value = initsig[iline]
     parinfo[isoff].limited = [1B,1B]
     parinfo[isoff].limits = siglim
     parinfo[isoff].tied = STRING('P[',ppoff+2+halpha*3,']',format='(A,I0,A)')
     if (linename[iline] eq 'Halpha') then parinfo[isoff].tied = ''
     if keyword_set(sigfix) then if sigfix[0] then parinfo[isoff].fixed=1B
    
  endfor

  linea = where(linename eq '[SII]6716')
  lineb = where(linename eq '[SII]6731')
  parinfo[ppoff+lineb*3].tied = 'P['+string(ppoff+linea*3,format='(I0)')+$
                                ']/P[3]'
  linea = where(linename eq '[OI]6364')
  lineb = where(linename eq '[OI]6300')
  parinfo[ppoff+linea*3].tied = 'P['+string(ppoff+lineb*3,format='(I0)')+']/3.0d'
  linea = where(linename eq '[NII]6548')
  lineb = where(linename eq '[NII]6583')
  parinfo[ppoff+linea*3].tied = 'P['+string(ppoff+lineb*3,format='(I0)')+']/3.0d'
  if keyword_set(n2hafix) then begin
     if n2hafix[0] gt 0 then $
        parinfo[ppoff+lineb*3].tied = $
     'P['+string(ppoff+halpha*3,format='(I0)')+']*'+$
     string(n2hafix[0],format='(D0.4)')
  endif

; other velocity components
  for i=1,ncomp-1 do begin

;    redshift
     parinfo[ppoff0-1+i].value = z.gas[i]
     parinfo[ppoff0-1+i].limited = [1B,1B]
     parinfo[ppoff0-1+i].limits = [z.gas[i]-0.003,z.gas[i]+0.003]
;     parinfo[ppoff0-1+i].step = 0.0001d
     parinfo[ppoff0-1+i].mpprint = 1b
     if keyword_set(zfix) then if zfix[i] then parinfo[ppoff0-1+i].fixed=1B

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
        parinfo[iwoff].tied = $
           string(linelambda[iline],'*(1d + P[',ppoff0-1+i,'])',$
                  format='(D0.2,A0,I0,A0)')
    
        parinfo[isoff].value = initsig[iline,i]
        parinfo[isoff].limited = [1B,1B]
        parinfo[isoff].limits = siglim
        parinfo[isoff].tied = $
           STRING("P[",soff+halpha*3,"]",format='(A0,I0,A0)')
        if (linename[iline] eq 'Halpha') then parinfo[isoff].tied = ''
        if keyword_set(sigfix) then if sigfix[i] then $
           parinfo[isoff].fixed=1B

     endfor

     linea = where(linename eq '[SII]6716')
     lineb = where(linename eq '[SII]6731')
     parinfo[foff+lineb*3].tied = 'P['+string(foff+linea*3,$
                                               format='(I0)')+']/P[3]'
     linea = where(linename eq '[OI]6364')
     lineb = where(linename eq '[OI]6300')
     parinfo[foff+linea*3].tied = 'P['+string(foff+lineb*3,format='(I0)')+']/3.0d'
     linea = where(linename eq '[NII]6548')
     lineb = where(linename eq '[NII]6583')
     parinfo[foff+linea*3].tied = 'P['+string(foff+lineb*3,format='(I0)')+']/3.0d'
     if keyword_set(n2hafix) then $
        if n2hafix[i] gt 0 then $
           parinfo[foff+lineb*3].tied = $
        'P['+string(foff+halpha*3,format='(I0)')+']*'+$
        string(n2hafix[i],format='(D0.4)')

  endfor

  return, parinfo

end
