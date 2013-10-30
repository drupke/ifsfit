function gmos_initlinelist,addnad=addnad
;
; History
;   09jul08  DSNR  created
;   10nov04  DSNR  wavelengths corrected to NIST values
;
; Wavelength sources:
;   1. Mappings III
;   2. NIST ASD ([SIII] 6312)
; 

  label = ['[OI]6300','[OI]6364',$
           '[NII]6548','Halpha','[NII]6583',$
           '[SII]6716','[SII]6731']

  wave = double([6300.30, 6363.78, $
                 6548.05, 6562.80, 6583.45, $
                 6716.44, 6730.82])

  if keyword_set(addnad) then begin
     label = ['NaD2','NaD1',label]
     wave = [5889.95d,5895.92d,wave]
  endif

  linelist = {label:label, wave:wave}

  return,linelist

end
