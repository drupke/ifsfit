;
; History
;  12jan02  DSNR  created
;
pro uhsf_printweq,instr,outfile,row,col,init=init
  
  if keyword_set(init) then begin

     openw,unit,outfile,/get_lun
     printf,unit,'#Row','Col','Weq(A)',format='(A-4,A4,A10)'

  endif else begin
     
; Model stellar continuum
     wave = instr.wave
     nlam = n_elements(wave)
;   polynomial near Ha/[NII]
     ypoly = dblarr(n_elements(wave))
     npoly = instr.param[4]
     wavelo = instr.param[5]
     wavehi = instr.param[6]
     n2hawave = where(wave ge wavelo AND wave le wavehi)
     ypoly[n2hawave] = $
        poly(wave[n2hawave]-mean(wave[n2hawave]),instr.param[7:7+npoly])
     specstars = instr.spec - instr.specfit + ypoly
     modstars = instr.spec - instr.spec_nocnt + ypoly

; Model Halpha emission line spectrum
     modhalines = dblarr(nlam)
     ncomp = instr.param[1]
     for i=1,ncomp do modhalines += uhsf_cmplin(instr,'Halpha',i)

; Equivalent width
     dwave = wave[1:nlam-1] - wave[0:nlam-2]
     weq = total(-modhalines[1:nlam-1]/modstars[1:nlam-1]*dwave)

     openu,unit,outfile,/get_lun,/append
     printf,unit,row,col,weq,format='(2I4,D10.4)'
     
  endelse

  free_lun,unit

end
