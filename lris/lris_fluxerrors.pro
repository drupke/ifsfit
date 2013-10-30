;
; History
;  10may17  DSNR  created
;

function lris_fluxerrors,instr

;  continuum = instr.spec_nocnt
  continuum = instr.spec_nocnt - instr.specfit

  ppoff = instr.param[0]
  ncomp = instr.param[1]

  nlines = fix((n_elements(instr.param)-ppoff)/3/ncomp)
  fluxind = indgen(nlines)*3 + ppoff
  waveind = fluxind+1
  sigind = fluxind+2
  
  fluxpkerr = dblarr(nlines,ncomp)

  for i=0,nlines-1 do begin
     for j=0,ncomp-1 do begin
        dind = j*nlines*3
        sigmult = 2d
        wavelo = instr.param[waveind[i]+dind]-$
                 sigmult*instr.param[sigind[i]+dind]
        wavehi = instr.param[waveind[i]+dind]+$
                 sigmult*instr.param[sigind[i]+dind]
        iran = where(instr.wave ge wavelo AND instr.wave lt wavehi,ct)
        if ct gt 0 then $
           fluxpkerr[i,j] = sqrt(total(continuum[iran]^2d)/n_elements(iran))
;        print,wavelo,wavehi,instr.param[fluxind[i]+dind],fluxpkerr[i,j]
     endfor
  endfor

  fluxpkerr = reform(fluxpkerr,nlines*ncomp)

  return,fluxpkerr

end
