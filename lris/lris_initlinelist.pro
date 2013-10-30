function lris_initlinelist,strong=strong
;
; History
;   09may26  DSNR  created
;   10nov04  DSNR  wavelengths corrected to NIST values
;
; Wavelength sources:
;   1. Mappings III
;   2. NIST ASD ([SIII] 6312)
; 

  labels = ['[OII]3726', '[OII]3729', 'H12', 'H11', 'H10', 'H9', $
            '[NeIII]3869', 'HeI3889', 'H8', 'HeI3965', '[NeIII]3967', $
            'Hepsilon', 'HeI4026', 'Hdelta', 'Hgamma', '[OIII]4363', $
            'HeI4472', 'Hbeta', 'HeI4922', '[OIII]4959', '[OIII]5007', $
            '[NI]5198', '[NI]5200', '[NII]5755', 'HeI5876', $
            '[OI]6300', '[SIII]6312', '[OI]6364', '[NII]6548', $
            'Halpha', '[NII]6583', 'HeI6678', '[SII]6716', $
            '[SII]6731', 'HeI7065', '[ArIII]7136', '[OII]7319', $
            '[OII]7320', '[OII]7330', '[OII]7331']

  waves = double([3726.04, 3728.80, 3750.15, 3770.63, 3797.90, $
                  3835.38, 3868.76, 3888.65, 3889.05, 3964.73, $
                  3967.47, 3970.07, 4026.19, 4101.73, 4340.46, $
                  4363.21, 4471.48, 4861.32, 4921.93, 4958.91, $
                  5006.84, 5197.90, 5200.26, 5754.59, $
                  5875.62, 6300.30, 6312.06, 6363.78, 6548.05, $
                  6562.80, 6583.45, 6678.15, 6716.44, 6730.82, $
                  7065.19, 7135.79, 7318.92, 7319.99, 7329.67, $
                  7330.73])

  if keyword_set(strong) then begin

     labels = ['[OII]3726','[OII]3729','Hgamma','Hbeta','[OIII]4959',$
               '[OIII]5007',$
               '[OI]6300','[OI]6364','[NII]6548','Halpha', '[NII]6583',$
               '[SII]6716','[SII]6731']

     waves = double([3726.04,3728.80,4340.46,4861.32,4958.91,5006.84,$
                     6300.30,6363.78,6548.05,6562.80,6583.45,$
                     6716.44, 6730.82])

  endif

  linelist = { label:labels, wave:waves }

  return,linelist

end
