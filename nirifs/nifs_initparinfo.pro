function nifs_initparinfo,linename,linelambda,linelambdaz,initflux,$
                          initsig,z,siglim=siglim,$
                          zfix=zfix,sigfix=sigfix,doubleline=doubleline
;
; 13mar05  DSNR  created
;

; Sigma = 1.6 A corresponds to R = 5290 (from NIFS website) at 2um
; Sigma = 85 A corresponds to 3000 km/s FWHM at 2um
  if ~ keyword_set(siglim) then siglim=[1.6d,85d]

; Number of velocity components
  ncomp = n_elements(z.gas)

; Number of initial parameters before Gaussian parameters begin
  ppoff0 = 3
  ppoff = ppoff0 + ncomp-1

; Number of emission lines to fit
  nline = N_ELEMENTS(linename)

  parinfo = REPLICATE({value:0d, fixed:0b, limited:[0B,0B], tied:'', $
                       limits:[0d,0d], step:0d, mpprint:0b, mpside:2}, $
                      ppoff+ncomp*(nline*3))

  paa = WHERE(linename eq 'Paa')
  h2s1 = WHERE(linename eq 'H2_10_S1')

; Number of initial parameters before Gaussian parameters begin
  parinfo[0].value = ppoff
  parinfo[0].fixed = 1B

; Number of velocity components
  parinfo[1].value = ncomp
  parinfo[1].fixed = 1B

; redshift of narrow component of recombination lines
  parinfo[2].value   = z.gas[0]
  parinfo[2].limited = [1B,1B]
  parinfo[2].limits  = [z.gas[0]-0.003,z.gas[0]+0.003]
  parinfo[2].mpprint = 1b
  if keyword_set(zfix) then if zfix[0] then parinfo[2].fixed=1B

;; ; H2 1-0 S(3)/S(1)
;;   parinfo[3].value   = 1d
;;   parinfo[3].limited = [1B,1B]
;;   parinfo[3].limits  = [0,1.906d]
;;   parinfo[3].mpprint = 1b

;; ; H2 1-0 S(3)/S()2
;;   parinfo[4].value   = 2d
;;   parinfo[4].limited = [1B,1B]
;;   parinfo[4].limits  = [0,3.878d]
;;   parinfo[4].mpprint = 1b

;; ; H2 1-0 S(3)/S(1)
;;   parinfo[5].value   = 1d
;;   parinfo[5].limited = [1B,1B]
;;   parinfo[5].limits  = [0,1.906d]
;;   parinfo[5].mpprint = 1b

;; ; H2 1-0 S(3)/S(2)
;;   parinfo[6].value   = 2d
;;   parinfo[6].limited = [1B,1B]
;;   parinfo[6].limits  = [0,3.878d]
;;   parinfo[6].mpprint = 1b

; main velocity component
  for iline=0, nline-1 do begin

     ifoff = ppoff+iline*3
     iwoff = ifoff+1
     isoff = ifoff+2
      
     parinfo[ifoff].value = initflux[iline]
     parinfo[ifoff].limited[0] = 1B
     parinfo[ifoff].limits[0]  = [0d]
     if (linename[iline] eq 'H2_10_S1' OR $
         linename[iline] eq 'H2_10_S2' OR $
         linename[iline] eq 'H2_10_S3' OR $
         linename[iline] eq 'H2_10_S4' OR $
         linename[iline] eq 'H2_10_S5') then begin
        parinfo[ifoff].value = 0d
        parinfo[ifoff].fixed = 1b
     endif
    
     parinfo[iwoff].value = linelambdaz[iline,0]
     parinfo[iwoff].tied = string(linelambda[iline],format='(D0.2)')+$
                           '*(1d +P[2])'
    
     parinfo[isoff].value = initsig[iline]
     parinfo[isoff].limited = [1B,1B]
     parinfo[isoff].limits = siglim
     parinfo[isoff].tied = STRING('P[',ppoff+2+paa*3,']',format='(A,I0,A)')
     if (linename[iline] eq 'Paa') then parinfo[isoff].tied = ''
     if keyword_set(sigfix) then if sigfix[0] then parinfo[isoff].fixed=1B
    
  endfor

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
        if doubleline eq 'paa' AND i eq 2 then begin
           if (linename[iline] eq 'H2_10_S1' OR $
               linename[iline] eq 'H2_10_S2' OR $
               linename[iline] eq 'H2_10_S3' OR $
               linename[iline] eq 'H2_10_S4' OR $
               linename[iline] eq 'H2_10_S5') then begin
              parinfo[ifoff].value = 0d
              parinfo[ifoff].fixed = 1b
           endif
        endif else begin
           if (linename[iline] ne 'H2_10_S1' AND $
               linename[iline] ne 'H2_10_S2' AND $
               linename[iline] ne 'H2_10_S3' AND $
               linename[iline] ne 'H2_10_S4' AND $
               linename[iline] ne 'H2_10_S5') then begin
              parinfo[ifoff].value = 0d
              parinfo[ifoff].fixed = 1b
           endif

; To remove consideration of H_2 S(2) in finding second component ...
           if i eq 2 AND linename[iline] eq 'H2_10_S2' then begin
              parinfo[ifoff].value = 0d
              parinfo[ifoff].fixed = 1b
           endif

        endelse
    
        parinfo[iwoff].value = linelambdaz[iline,i]
        parinfo[iwoff].tied = $
           string(linelambda[iline],'*(1d + P[',ppoff0-1+i,'])',$
                  format='(D0.2,A0,I0,A0)')
    
        parinfo[isoff].value = initsig[iline,i]
        parinfo[isoff].limited = [1B,1B]
        parinfo[isoff].limits = siglim
        if doubleline eq 'paa' AND i eq 2 then begin
           parinfo[isoff].tied = STRING('P[',soff+paa*3,']',format='(A,I0,A)')
           if (linename[iline] eq 'Paa') then parinfo[isoff].tied = ''
        endif else begin
           parinfo[isoff].tied = $
              STRING("P[",soff+h2s1*3,"]",format='(A0,I0,A0)')
           if (linename[iline] eq 'H2_10_S1') then parinfo[isoff].tied = ''
        endelse
        if keyword_set(sigfix) then if sigfix[i] then $
           parinfo[isoff].fixed=1B

     endfor

     ;; if i eq 1 then begin
     ;;    p_h2s3s1 = '3'
     ;;    p_h2s3s2 = '4'
     ;; endif else begin 
     ;;    p_h2s3s1 = '5'
     ;;    p_h2s3s2 = '6'
     ;; endelse
     ;; linea = where(linename eq 'H2_10_S3')
     ;; lineb = where(linename eq 'H2_10_S1')
     ;; parinfo[foff+lineb*3].tied = 'P['+string(foff+linea*3,$
     ;;                                          format='(I0)')+$
     ;;                              ']/P['+p_h2s3s1+']'
     ;; linea = where(linename eq 'H2_10_S3')
     ;; lineb = where(linename eq 'H2_10_S2')
     ;; parinfo[foff+lineb*3].tied = 'P['+string(foff+linea*3,$
     ;;                                          format='(I0)')+$
     ;;                              ']/P['+p_h2s3s2+']'
     
  endfor

  return, parinfo

end
