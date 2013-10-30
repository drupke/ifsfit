;
; History
;  12jan03  DSNR  created separate program

function gmos_componeline,instr,line,comp,center=center
  
  iline = where(instr.linelabel eq line,ct)
  ppoff = instr.param[0]
  ncomp = instr.param[1]
  ppoff0 = ppoff - (ncomp-1)

  nline = n_elements(instr.linelabel)
  indices = ppoff+(comp-1)*nline*3+iline*3
  indices = indices[0] + indgen(3)
;  print,comp,indices[0],instr.param[indices]
  flux = gaussian(instr.wave,instr.param[indices],/double)

; central wavelength
  tmp = instr.param[indices]
  center = tmp[1]

  return,flux

end
