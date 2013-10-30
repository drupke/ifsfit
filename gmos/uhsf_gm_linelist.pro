;
; History
;   09jul08  DSNR  created
;   10nov04  DSNR  wavelengths corrected to NIST values
;   13sep    DSNR  added more lines for 1-slit GMOS config.
;
; Wavelength sources:
;   1. Mappings III
;   2. HeI 4686 -- found random wavelength compilation. In NIST, shows
;                  up as a complex of many fine structure lines with
;                  range of several tenths of an A.
; 
function uhsf_gm_linelist,addnad=addnad,addo1_5577=addo1_5577,$
                          twoslit=twoslit,felines=felines

  if not keyword_set(twoslit) then begin

     label = ['HeII4686', $
              'HeI6678', 'HeI7065', $
              '[OIII]4959', '[OIII]5007', $
              'Hbeta', '[NI]5198', '[NI]5200', $
              '[OI]6300', '[OI]6364', '[NII]6548', $
              'Halpha', '[NII]6583', '[SII]6716', $
              '[SII]6731']
     
     wave = double([4686.7, $
                    6678.15, 7065.19, $
                    4958.91, 5006.84, $
                    4861.32, 5197.90, 5200.26, $
                    6300.30, 6363.78, 6548.05, $
                    6562.80, 6583.45, 6716.44, $
                    6730.82])

  endif else begin

     label = ['[OI]6300', '[OI]6364', '[NII]6548', $
              'Halpha', '[NII]6583', '[SII]6716', $
              '[SII]6731']
     
     wave = double([6300.30, 6363.78, 6548.05, $
                    6562.80, 6583.45, 6716.44, $
                    6730.82])

  endelse

  if keyword_set(addnad) then begin
     label = ['NaD2','NaD1',label]
     wave = [5889.95d,5895.92d,wave]
  endif

  if keyword_set(addo1_5577) then begin
     label = ['[OI]5577',label]
     wave = [5577.34,wave]
  endif

  if keyword_set(felines) then begin
     label = ['[FeVII]5159','[FeVII]5721','[FeVII]6087',$
              '[FeX]6375',label]
     wave = [5158.89,5720.7,6087.0,6374.51,wave]
  endif
     
  linelist = {label:label, wave:wave}

  return,linelist

end
