function uhsf_gm_initpar,linename,linelambda,linelambdaz,linetie,$
                         initflux,initsig,ncomp,$
                         siglim=siglim,sigfix=sigfix,zfix=zfix

;
; 09jun01  DSNR  created
; 09jun08  DSNR  added multiple components
; 10may27  DSNR  re-written to fit in observed frame
; 13sep13  DSNR  re-written to allow more than one common redshift
;

  c = 299792.458d

; Sigma limits
  gmosres = 3000d
  if ~ keyword_set(siglim) then siglim=[299792d/gmosres/2.35d,2000d]

; Number of emission lines to fit
  nline = N_ELEMENTS(linename)

  maxncomp = max(ncomp)
  
; Number of initial parameters before Gaussian parameters begin
  lratlim = 4
  ppoff0 = 2
  ppoff = ppoff0 + maxncomp*lratlim
  
  parinfo = REPLICATE({value:0d, fixed:0b, limited:[0B,0B], tied:'', $
                       limits:[0d,0d], step:0d, mpprint:0b, mpside:2}, $
                      ppoff+maxncomp*(nline*3))

; Number of initial parameters before Gaussian parameters begin
  parinfo[0].value = ppoff
  parinfo[0].fixed = 1B

; Number of velocity components
  parinfo[1].value = maxncomp
  parinfo[1].fixed = 1B

; [SII] ratio
  ilratlim = 0
  is2a = where(linename eq '[SII]6716')
  is2b = where(linename eq '[SII]6731')
  tmp_ncomp = ncomp[is2a]
  tmp_ncomp = tmp_ncomp[0]
  if tmp_ncomp gt 0 then begin
     ip1 = ppoff0 + ilratlim*maxncomp
     ip2 = ip1+tmp_ncomp-1
     parinfo[ip1:ip2].value = transpose(initflux[is2a,0:tmp_ncomp-1]/$
                                        initflux[is2b,0:tmp_ncomp-1])
     parinfo[ip1:ip2].limited = rebin([1b,1b],2,tmp_ncomp)
     parinfo[ip1:ip2].limits  = rebin([0.44d,1.43d],2,tmp_ncomp)
     for i=0,tmp_ncomp-1 do begin
;       case of both lines getting zero-ed        
        if finite(parinfo[ip1+i].value,/nan) then parinfo[ip1+i].value = $
           (parinfo[ip1+i].limits[0] + parinfo[ip1+i].limits[0])/2d
;       case of pegging at or exceeding upper limit
        if parinfo[ip1+i].value ge parinfo[ip1+i].limits[1] then $
           parinfo[ip1+i].value = parinfo[ip1+i].limits[1] - $
                                  (parinfo[ip1+i].limits[1] - $
                                   parinfo[ip1+i].limits[0])*0.1
;       case of pegging at or dipping below lower limit
        if parinfo[ip1+i].value le parinfo[ip1+i].limits[0] then $
           parinfo[ip1+i].value = parinfo[ip1+i].limits[0] + $
                                  (parinfo[ip1+i].limits[1] - $
                                   parinfo[ip1+i].limits[0])*0.1
     endfor
  endif
  
; [NI] ratio
; See Ferland+12 for collisional case, Bautista99 for other cases 
  ilratlim = 1
  in1a = where(linename eq '[NI]5198')
  in1b = where(linename eq '[NI]5200')
  tmp_ncomp = ncomp[in1a]
  tmp_ncomp = tmp_ncomp[0]
  if ncomp[in1a] gt 0 then begin
     ip1 = ppoff0 + ilratlim*maxncomp
     ip2 = ip1+tmp_ncomp-1
     parinfo[ip1:ip2].value = transpose(initflux[in1b,0:tmp_ncomp-1]/$
                                        initflux[in1a,0:tmp_ncomp-1])
     parinfo[ip1:ip2].limited = rebin([1b,1b],2,tmp_ncomp)
     parinfo[ip1:ip2].limits  = rebin([0.6d,3d],2,tmp_ncomp)
     for i=0,tmp_ncomp-1 do begin
        if finite(parinfo[ip1+i].value,/nan) then parinfo[ip1+i].value = $
           (parinfo[ip1+i].limits[0] + parinfo[ip1+i].limits[0])/2d
        if parinfo[ip1+i].value ge parinfo[ip1+i].limits[1] then $
           parinfo[ip1+i].value = parinfo[ip1+i].limits[1] - $
                                  (parinfo[ip1+i].limits[1] - $
                                   parinfo[ip1+i].limits[0])*0.1
        if parinfo[ip1+i].value le parinfo[ip1+i].limits[0] then $
           parinfo[ip1+i].value = parinfo[ip1+i].limits[0] + $
                                  (parinfo[ip1+i].limits[1] - $
                                   parinfo[ip1+i].limits[0])*0.1
     endfor
  endif

; [NII]/Ha ratio
  ilratlim = 2
  iha = where(linename eq 'Halpha')
  in2 = where(linename eq '[NII]6583')
  tmp_ncomp = ncomp[iha]
  tmp_ncomp = tmp_ncomp[0]
  if tmp_ncomp gt 0 then begin
     ip1 = ppoff0 + ilratlim*maxncomp
     ip2 = ip1 + tmp_ncomp - 1
     parinfo[ip1:ip2].value = transpose(initflux[in2,0:tmp_ncomp-1]/$
                                        initflux[iha,0:tmp_ncomp-1])
     parinfo[ip1:ip2].limited = rebin([0B,1B],2,tmp_ncomp)
     parinfo[ip1:ip2].limits  = rebin([0d,4d],2,tmp_ncomp)
     for i=0,tmp_ncomp-1 do begin
        if finite(parinfo[ip1+i].value,/nan) then parinfo[ip1+i].value = $
           (parinfo[ip1+i].limits[0] + parinfo[ip1+i].limits[0])/2d
        if parinfo[ip1+i].value ge parinfo[ip1+i].limits[1] then $
           parinfo[ip1+i].value = parinfo[ip1+i].limits[1] - $
                                  (parinfo[ip1+i].limits[1] - $
                                   parinfo[ip1+i].limits[0])*0.1
        if parinfo[ip1+i].value le parinfo[ip1+i].limits[0] then $
           parinfo[ip1+i].value = parinfo[ip1+i].limits[0] + $
                                  (parinfo[ip1+i].limits[1] - $
                                   parinfo[ip1+i].limits[0])*0.1
     endfor
  endif

; Ha/Hb ratio
  ilratlim = 3
  iha = where(linename eq 'Halpha')
  ihb = where(linename eq 'Hbeta')
  tmp_ncomp = ncomp[iha]
  tmp_ncomp = tmp_ncomp[0]
  if tmp_ncomp gt 0 then begin
     ip1 = ppoff0 + ilratlim*maxncomp
     ip2 = ip1 + tmp_ncomp - 1
     parinfo[ip1:ip2].value = transpose(initflux[iha,0:tmp_ncomp-1]/$
                                        initflux[ihb,0:tmp_ncomp-1])
     parinfo[ip1:ip2].limited = rebin([1B,0B],2,tmp_ncomp)
     parinfo[ip1:ip2].limits  = rebin([2.86d,100d],2,tmp_ncomp)
     for i=0,tmp_ncomp-1 do begin
        if finite(parinfo[ip1+i].value,/nan) then parinfo[ip1+i].value = $
           parinfo[ip1+i].limits[0]
        if parinfo[ip1+i].value ge parinfo[ip1+i].limits[1] then $
           parinfo[ip1+i].value = parinfo[ip1+i].limits[1] - $
                                  (parinfo[ip1+i].limits[1] - $
                                   parinfo[ip1+i].limits[0])*0.1
        if parinfo[ip1+i].value le parinfo[ip1+i].limits[0] then $
           parinfo[ip1+i].value = parinfo[ip1+i].limits[0] + $
                                  (parinfo[ip1+i].limits[1] - $
                                   parinfo[ip1+i].limits[0])*0.1
     endfor
  endif

; cycle through velocity components
  for i=0,maxncomp-1 do begin

;    index offsets for this component
     foff = ppoff+i*nline*3
     woff = foff+1
     soff = foff+2

;    cycle through lines
     for iline=0, nline-1 do begin


;       indices
        ifoff = foff + iline*3
        iwoff = woff + iline*3
        isoff = soff + iline*3

        if i+1 gt ncomp[iline] then begin 
           parinfo[ifoff].value = 0d
           parinfo[iwoff].value = 0d
           parinfo[isoff].value = 0d
           parinfo[ifoff].fixed = 1b
           parinfo[iwoff].fixed = 1b
           parinfo[isoff].fixed = 1b
        endif else begin

;          initial values
           parinfo[ifoff].value = initflux[iline,i]
           parinfo[iwoff].value = linelambdaz[iline,i]
           parinfo[isoff].value = initsig[iline,i]
;          limits
           parinfo[ifoff].limited[0] = 1B
           parinfo[ifoff].limits[0]  = 0d
           parinfo[iwoff].limited = [1B,1B]
           parinfo[iwoff].limits[0] = linelambdaz[iline,i]*0.997d
           parinfo[iwoff].limits[1] = linelambdaz[iline,i]*1.003d
           parinfo[isoff].limited = [1B,1B]
           parinfo[isoff].limits = siglim
;          ties
           if (linename[iline] eq linetie[iline]) then begin
              parinfo[iwoff].tied = ''
              parinfo[isoff].tied = ''
           endif else begin
              indtie = where(linename eq linetie[iline])
              parinfo[iwoff].tied = $
                 string(linelambda[iline],'/',linelambda[indtie],$
                        '* P[',woff+indtie*3,']',$
                        format='(D0.2,A0,D0.2,A0,I0,A0)') 
              parinfo[isoff].tied = $
                 string('P[',soff+indtie*3,']',format='(A0,I0,A0)')   
           endelse
;          fixed/free
           if keyword_set(sigfix) then if sigfix[i] then $
              parinfo[isoff].fixed=1B
           if keyword_set(zfix) then if zfix[i] then parinfo[iwoff].fixed=1B

        endelse

     endfor

     ilratlim = 0
     linea = where(linename eq '[SII]6716')
     lineb = where(linename eq '[SII]6731')
     parinfo[foff+lineb*3].tied = 'P['+$
                                  string(foff+linea*3,$
                                         format='(I0)')+']/P['+$
                                  string(ppoff0+maxncomp*ilratlim+i,$
                                         format='(I0)')+']'
     ilratlim = 1
     linea = where(linename eq '[NI]5198')
     lineb = where(linename eq '[NI]5200')
     parinfo[foff+lineb*3].tied = 'P['+$
                                  string(foff+linea*3,$
                                         format='(I0)')+']*P['+$
                                  string(ppoff0+maxncomp*ilratlim+i,$
                                         format='(I0)')+']'
     ilratlim = 2
     linea = where(linename eq 'Halpha')
     lineb = where(linename eq '[NII]6583')
     parinfo[foff+lineb*3].tied = 'P['+$
                                  string(foff+linea*3,$
                                         format='(I0)')+']*P['+$
                                  string(ppoff0+maxncomp*ilratlim+i,$
                                         format='(I0)')+']'
     ilratlim = 3
     linea = where(linename eq 'Halpha')
     lineb = where(linename eq 'Hbeta')
     parinfo[foff+lineb*3].tied = 'P['+$
                                  string(foff+linea*3,$
                                         format='(I0)')+']/P['+$
                                  string(ppoff0+maxncomp*ilratlim+i,$
                                         format='(I0)')+']'


     linea = where(linename eq '[OIII]4959')
     lineb = where(linename eq '[OIII]5007')
     parinfo[foff+linea*3].tied = 'P['+$
                                  string(foff+lineb*3,$
                                         format='(I0)')+']/3.0d'
     linea = where(linename eq '[OI]6364')
     lineb = where(linename eq '[OI]6300')
     parinfo[foff+linea*3].tied = 'P['+$
                                  string(foff+lineb*3,$
                                         format='(I0)')+']/3.0d'
     linea = where(linename eq '[NII]6548')
     lineb = where(linename eq '[NII]6583')
     parinfo[foff+linea*3].tied = 'P['+ $
                                  string(foff+lineb*3,$
                                         format='(I0)')+']/3.0d'

  endfor

; Check parinit initial values vs. limits
  badpar = where((parinfo.limited[0] AND $
                  parinfo.value lt parinfo.limits[0]) OR $
                 (parinfo.limited[1] AND $
                  parinfo.value gt parinfo.limits[1]),ct)
  if ct gt 0 then begin
     print,'UHSF_GM_INITPAR: Initial values are outside limits.'
     print,'UHSF_GM_INITPAR: Offending parameters:'
     print,'Index','Value','Lower limit','Upper limit',$
           format='(A7,3A15)'
     for i=0,ct-1 do begin
        j = badpar[i]
        print,j,parinfo[j].value,parinfo[j].limits[0],$
              parinfo[j].limits[1],format='(I7,3E15.6)'
     endfor
     return,0
  endif else begin
     return,parinfo
  endelse

end
