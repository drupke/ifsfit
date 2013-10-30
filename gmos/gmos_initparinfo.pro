function gmos_initparinfo,linename,linelambda,initflux,initsig,$
                          mcomp
;
; 09jun01  DSNR  created
; 09jun08  DSNR  added multiple components
;

; Sigma limit.  0.8d is approximately resolution for GMOS R600
; grating.  10 A in sigma is ~1000 km/s FWHM.
  siglim = [0d,20d]

; Number of velocity components
  ncomp = n_elements(mcomp)

; Number of initial parameters before Gaussian parameters begin
  ppoff0 = 3
  ppoff = ppoff0 + 2*(ncomp-1)

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
  parinfo[2].value   = 0d
  parinfo[2].limited = [1B,1B]
  parinfo[2].limits  = [-0.003,0.003]
  parinfo[2].step = 0.0001d
  parinfo[2].mpprint = 1b

;; ; [SII] ratio
;;   parinfo[3].value   = 1.3d
;;   parinfo[3].limited = [1B,1B]
;;   parinfo[3].limits  = [0.44d,1.43d]

;; ; [OI] z offset
;;   parinfo[3].value = [0d]
;;   parinfo[3].limited = [1B,1B]
;;   parinfo[3].limits  = [-0.003,0.003]
;;   parinfo[3].step = 0.0001d
;;   parinfo[3].mpprint = 1b
  
; main velocity component
  for iline=0, nline-1 do begin

     ifoff = ppoff+iline*3
     iwoff = ifoff+1
     isoff = ifoff+2
      
     parinfo[ifoff].value = initflux[iline]
     parinfo[ifoff].limited[0] = 1B
     parinfo[ifoff].limits[0]  = [0d] ;;, 2d*MAX(flux)]
    
     parinfo[iwoff].value = linelambda[iline]
     parinfo[iwoff].tied = string(linelambda[iline],format='(D0.2)')+$
                           '*(1d +P[2])'
;;      if (linename[iline] eq '[OI]6300' OR linename[iline] eq '[OI]6364') then $
;;         parinfo[iwoff].tied = string(linelambda[iline],format='(D0.2)')+$
;;                               '*(1d +P[2]+P[3])'
    
     parinfo[isoff].value = initsig[iline]
     parinfo[isoff].limited = [1B,1B]
     parinfo[isoff].limits = siglim
     parinfo[isoff].tied = STRING('P[',ppoff+2+halpha*3,']',format='(A,I0,A)')
     if (linename[iline] eq 'Halpha') then parinfo[isoff].tied = ''
    
  endfor

;;   linea = where(linename eq '[SII]6716')
;;   lineb = where(linename eq '[SII]6731')
;;   parinfo[ppoff+lineb*3].tied = 'P['+string(ppoff+linea*3,format='(I0)')+$
;;                                 ']/P[3]'
  linea = where(linename eq '[OI]6364')
  lineb = where(linename eq '[OI]6300')
  parinfo[ppoff+linea*3].tied = 'P['+string(ppoff+lineb*3,format='(I0)')+']/3.0d'
  linea = where(linename eq '[NII]6548')
  lineb = where(linename eq '[NII]6583')
  parinfo[ppoff+linea*3].tied = 'P['+string(ppoff+lineb*3,format='(I0)')+']/3.0d'

; other velocity components
  for i=1,ncomp-1 do begin

     izoff1 = ppoff0+2*i-2
     izoff2 = ppoff0+2*i-1

;    redshift diff
     parinfo[izoff1].value = mcomp[i]-mcomp[0]
     parinfo[izoff1].limited = [1B,1B]
     parinfo[izoff1].limits = [-0.002,0.002]
     parinfo[izoff1].step = 0.00001d

;    redshift
     parinfo[izoff2].value = mcomp[i]
     parinfo[izoff2].tied = 'P['+string(izoff1,format='(I0)')+']+P[2]'
     parinfo[izoff2].mpprint = 1b

     foff = ppoff+i*nline*3
     woff = foff+1
     soff = foff+2

;    line parameters
     for iline=0, nline-1 do begin

        ifoff = foff + iline*3
        iwoff = woff + iline*3
        isoff = soff + iline*3

        parinfo[ifoff].value = initflux[iline]*0.2d
        parinfo[ifoff].limited[0] = 1B
        parinfo[ifoff].limits[0]  = [0d]
    
        parinfo[iwoff].value = linelambda[iline]*(1d + mcomp[i])
        parinfo[iwoff].tied = string(linelambda[iline],'*(1d + P[',ppoff0-1+i,'])',$
                                     format='(D0.2,A0,I0,A0)')
;;         if (linename[iline] eq '[OI]6300' OR linename[iline] eq '[OI]6364') then $
;;            parinfo[iwoff].tied = string(linelambda[iline],'*(1d + P[',ppoff0-1+i,'] + P[3])',$
;;                                         format='(D0.2,A0,I0,A0)')
    
        parinfo[isoff].value = initsig[iline]
        parinfo[isoff].limited = [1B,1B]
        parinfo[isoff].limits = siglim
        parinfo[isoff].tied = STRING("P[",soff+halpha*3,"]",format='(A0,I0,A0)')
        if (linename[iline] eq 'Halpha') then parinfo[isoff].tied = ''

     endfor

;;      linea = where(linename eq '[SII]6716')
;;      lineb = where(linename eq '[SII]6731')
;;      parinfo[foff+lineb*3].tied = 'P['+string(foff+linea*3,$
;;                                                format='(I0)')+']/P[3]'
     linea = where(linename eq '[OI]6364')
     lineb = where(linename eq '[OI]6300')
     parinfo[foff+linea*3].tied = 'P['+string(foff+lineb*3,format='(I0)')+']/3.0d'
     linea = where(linename eq '[NII]6548')
     lineb = where(linename eq '[NII]6583')
     parinfo[foff+linea*3].tied = 'P['+string(foff+lineb*3,format='(I0)')+']/3.0d'

  endfor

  return, parinfo

end
