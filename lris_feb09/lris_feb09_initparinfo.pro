function lris_feb09_initparinfo,linename,linelambda,linelambdaz,$
                                initflux,initsig,z,$
                                stitchwave=stitchwave,zfix=zfix,$
                                floatcomp=floatcomp,initbcomp=initbcomp,$
                                initrcomp=initrcomp
;
; 09jun01  DSNR  created
; 09jun08  DSNR  added multiple components
; 09sep15  DSNR  modified for 09 Feb. data
; 10apr21  DSNR  re-written to fit in observed frame
;

  if ~ keyword_set(stitchwave) then stitchwave = 5600

; Number of velocity components
  ncomp = n_elements(z.gas)

; Number of initial parameters before Gaussian parameters begin
  ppoff0 = 16
  ppoff = ppoff0 + ncomp-1

; Number of emission lines to fit
  nline = N_ELEMENTS(linename)

  parinfo = REPLICATE({value:0d, fixed:0, limited:[0B,0B], tied:'', $
                       limits:[0D,0D], step:0d, relstep:0d, mpprint:0b}, $
                      ppoff+ncomp*(nline*3))

  halpha = WHERE(linename eq 'Halpha')

; Number of initial parameters before Gaussian parameters begin
  parinfo[0].value = ppoff
  parinfo[0].fixed = 1B

; Number of velocity components
  parinfo[1].value = ncomp
  parinfo[1].fixed = 1B  

; blue side redshift
  parinfo[2].value   = z.blue[0]
  parinfo[2].limited = [1B,1B]
  parinfo[2].limits  = [z.blue[0]-0.003,z.blue[0]+0.003]
  parinfo[2].step = 0.0001d
  parinfo[2].mpprint = 1b
  if keyword_set(zfix) then parinfo[2].tied = 'P[3]'

; red side redshift  
  parinfo[3].value   = z.gas[0]
  parinfo[3].limited = [1B,1B]
  parinfo[3].limits  = [z.gas[0]-0.003,z.gas[0]+0.003]
  parinfo[3].step = 0.0001d
  parinfo[3].mpprint = 1b

; blue side tilt  
  intercept = z.bluetilt
  parinfo[4].value   = 0d
  parinfo[4].limited = [1B,1B]
  parinfo[4].limits  = [-5d-7,5d-7]
  parinfo[4].mpprint = 1b
  if keyword_set(zfix) then parinfo[4].fixed = 1B

;
; References for relative line intensities of p1/p3 ions:
;  Osterbrock 1989 (AGN^2)
;  Osterbrock & Ferland 2005 (AGN^2), Figure 5.8
;
; [OII] ratio
  parinfo[5].value   = 1.3d
  parinfo[5].limited = [1B,1B]
  parinfo[5].limits  = [0.27d,1.48d]
; [SII] ratio
  parinfo[6].value   = 1.3d
  parinfo[6].limited = [1B,1B]
  parinfo[6].limits  = [0.44d,1.43d]

; Ratio between sigmas (blue over red)
  parinfo[7].value = 1d
  parinfo[7].limited = [1B,1B]
  parinfo[7].limits = [1d,1.4d]

; Broad component parameters:
; Ratio between broad component peak flux and narrow component peak flux
  parinfo[8].value = initbcomp[0] ;0.03d
  parinfo[8].limited = [1B,1B]
  parinfo[8].limits = [0d,0.2d]
; Broad component sigma, in A
  parinfo[9].value = initbcomp[1] ;5d
  parinfo[9].limited = [1B,1B]
  parinfo[9].limits = [8.5d,50d]

; Red component parameters:
; Difference between wavelength of red component and narrow component
  parinfo[10].value = initrcomp[0] ;3d
  parinfo[10].limited = [1B,1B]
  parinfo[10].limits = [1d,10d]
; Ratio between red component peak flux and narrow component peak flux
  parinfo[11].value = initrcomp[1] ;0.1d
  parinfo[11].limited = [1B,1B]
  parinfo[11].limits = [0d,0.5d]
; Red component sigma, in A
  parinfo[12].value = initrcomp[2] ;1d
  parinfo[12].limited = [1B,1B]
  parinfo[12].limits = [0.5d,5d]

; Floating broad component
  if keyword_set(floatcomp) then begin
;    Flux estimate
     parinfo[13].value = floatcomp[0]
     parinfo[13].limited = [1B,1B]
     parinfo[13].limits = floatcomp[0]*[0.1d,10d]
     parinfo[13].mpprint = 1b
;    Location estimate, in A
     parinfo[14].value = floatcomp[1]
     parinfo[14].limited = [1B,1B]
     parinfo[14].limits = floatcomp[1]-[15d,-15d]
     parinfo[14].mpprint = 1b
;    Width estimate
     parinfo[15].value = floatcomp[2]
     parinfo[15].limited = [1B,1B]
     parinfo[15].limits = floatcomp[2]*[0.5d,2d]
     parinfo[15].mpprint = 1b
  endif else begin
     parinfo[13].value = 0d
     parinfo[13].fixed = 1B
     parinfo[14].value = 6563d
     parinfo[14].fixed = 1B
     parinfo[15].value = 1d
     parinfo[15].fixed = 1B
  endelse


;
; main velocity component
;
  for iline=0, nline-1 do begin

    ifoff = ppoff+iline*3
    iwoff = ifoff+1
    isoff = ifoff+2
      
    parinfo[ifoff].value = initflux[iline]
    parinfo[ifoff].limited[0] = 1B
    parinfo[ifoff].limits[0]  = 0d
    
    parinfo[iwoff].value = linelambdaz[iline,0]
    if (linelambda[iline] lt stitchwave) then parinfo[iwoff].tied = $
       string(linelambda[iline],format='(D0.2)')+'*(1d + P[2] + P[4]*('+$
       string(linelambda[iline],'-',z.blueint,$
              format='(D0.2,A0,D0.2)')+'))' $
    else $
       parinfo[iwoff].tied = string(linelambda[iline],$
                                    format='(D0.2)')+'*(1d + P[3])'
    
    parinfo[isoff].value = initsig[iline]
    parinfo[isoff].limited = [1B,1B]
    parinfo[isoff].limits = [1d, 10d]
    if (linelambda[iline] lt stitchwave) then $
       parinfo[isoff].tied = STRING("P[7]*P[",ppoff+2+halpha*3,"]",$
       format='(A,I0,A)') $
    else if (linelambda[iline] gt stitchwave) then $
       parinfo[isoff].tied = STRING("P[",ppoff+2+halpha*3,"]",$
                                    format='(A,I0,A)')
    if (linename[iline] eq 'Halpha') then parinfo[isoff].tied = ''
    
  endfor

  linea = where(linename eq '[OII]3729')
  lineb = where(linename eq '[OII]3726')
  parinfo[ppoff+lineb*3].tied = 'P['+string(ppoff+linea*3,$
                                            format='(I0)')+$
                                ']/P[5]'
  linea = where(linename eq '[SII]6716')
  lineb = where(linename eq '[SII]6731')
  parinfo[ppoff+lineb*3].tied = 'P['+string(ppoff+linea*3,$
                                            format='(I0)')+$
                                ']/P[6]'

;
; References for relative line intensities of p2/p4 ions:
;  Storrey & Zeippen 2000
;  Leisy & Dennefeld 1996, as quoted by SZ00
;  Osterbrock 1989 (AGN^2), Tables 3.8, 3.10
;  Osterbrock & Ferland 2005 (AGN^2)
; #s are approximations; could refine later.  Exact obs/th values hard
; to track down.
;
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

; other velocity components
  for i=1,ncomp-1 do begin

;    redshift
     parinfo[ppoff0-1+i].value = z.gas[i]
     parinfo[ppoff0-1+i].limited = [1B,1B]
     parinfo[ppoff0-1+i].limits = [z.gas[i]-0.003,z.gas[i]+0.003]
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
        if (linelambda[iline] lt stitchwave) then $
           parinfo[iwoff].tied = $
           string(linelambda[iline],'*(1d +P[',ppoff0-1+i,$
                  ']+P[2]-P[3]+P[4]*(',linelambda[iline],'-',z.blueint,$
                  '))',format='(D0.2,A0,I0,A0,D0.2,A0,D0.2,A0)') $
        else $
          parinfo[iwoff].tied = $
            string(linelambda[iline],'*(1d + P[',ppoff0-1+i,'])',$
            format='(D0.2,A0,I0,A0)')
    
        parinfo[isoff].value = initsig[iline,i]
        parinfo[isoff].limited = [1B,1B]
        parinfo[isoff].limits = [1d,10d]
        if (linelambda[iline] lt stitchwave) then parinfo[isoff].tied = $
          STRING("P[7]*P[",soff+halpha*3,"]",format='(A0,I0,A0)') $
        else if (linelambda[iline] gt stitchwave) then parinfo[isoff].tied = $
          STRING("P[",soff+halpha*3,"]",format='(A0,I0,A0)')
        if (linename[iline] eq 'Halpha') then parinfo[isoff].tied = ''

     endfor

     linea = where(linename eq '[OII]3729')
     lineb = where(linename eq '[OII]3726')
     parinfo[foff+lineb*3].tied = 'P['+string(foff+linea*3,$
                                              format='(I0)')+']/P[5]'
     linea = where(linename eq '[SII]6716')
     lineb = where(linename eq '[SII]6731')
     parinfo[foff+lineb*3].tied = 'P['+string(foff+linea*3,$
                                              format='(I0)')+']/P[6]'

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

  endfor

; ; Check for bad parinfo:
;  for i=0,n_elements(parinfo)-1 do begin
;     if ((parinfo[i].limited[0] AND parinfo[i].value lt parinfo[i].limits[0]) OR $
;     (parinfo[i].limited[1] AND parinfo[i].value gt parinfo[i].limits[1])) then $
;     stop,i,parinfo[i].value,parinfo[i].limits
;  endfor

  return, parinfo

end
