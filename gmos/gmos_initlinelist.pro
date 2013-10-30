function gmos_initlinelist
;
; History
;   09jul08  DSNR  created
;
; Wavelength sources:
;   1. Mappings III
;   2. NIST ASD ([SIII] 6312)
; 

  nlines = 8
  linelist = { label:strarr(nlines), wave:dblarr(nlines) }

  linelist.label = ['HeI5876','[OI]6300','[OI]6364',$
                    '[NII]6548','Halpha','[NII]6583',$
                    '[SII]6716','[SII]6731']

;    GMOS
;;   linelist.wave = double([5875.60, 6300.20, 6363.67, $
;;                           6547.96, 6562.80, 6583.34, $
;;                           6716.31, 6730.68])
;    KPNO
  linelist.wave = double([5877.289, 6302.049, 6365.536, $
                          6549.86, 6564.60, 6585.27, $
                          6718.29, 6732.67])



  return,linelist
end
