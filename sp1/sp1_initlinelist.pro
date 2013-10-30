function sp1_initlinelist,strong=strong,vacuum=vacuum,quiet=quiet
;
; History
;   09nov24  DSNR  created
;   10nov04  DSNR  revised to allow for vacuum wavelengths in all lines
;                  vacuum and air wavelengths corrected to NIST values
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

  if keyword_set(vacuum) then begin

     if (~ quiet) then print,'SP1_INITLINELIST: ',$
                             'Using vacuum wavelengths.'

     waves = double([3727.10, 3729.86, 3751.22, 3771.70, 3798.98, $
                     3836.47, 3869.87, 3889.75, 3890.15, 3965.85, $
                     3968.59, 3971.20, 4027.33, 4102.89, 4341.68, $
                     4364.435, 4472.74, 4862.68, 4923.31, 4960.295, $
                     5008.240, 5199.35, 5201.71, 5756.19, $
                     5877.25, 6302.046, 6313.81, 6365.535, 6549.86, $
                     6564.61, 6585.27, 6680.00, 6718.29, 6732.67, $
                     7067.14, 7137.76, 7320.94, 7322.01, 7331.69, $
                     7332.75])

  endif else begin

     waves = double([3726.04, 3728.80, 3750.15, 3770.63, 3797.90, $
                     3835.38, 3868.76, 3888.65, 3889.05, 3964.73, $
                     3967.47, 3970.07, 4026.19, 4101.73, 4340.46, $
                     4363.21, 4471.48, 4861.32, 4921.93, 4958.91, $
                     5006.84, 5197.90, 5200.26, 5754.59, $
                     5875.62, 6300.30, 6312.06, 6363.78, 6548.05, $
                     6562.80, 6583.45, 6678.15, 6716.44, 6730.82, $
                     7065.19, 7135.79, 7318.92, 7319.99, 7329.67, $
                     7330.73])

  endelse

  if keyword_set(strong) then begin

     labels = ['[OII]3726','[OII]3729','Hbeta','[OIII]4959','[OIII]5007',$
               '[OI]6300','[OI]6364','[NII]6548','Halpha', '[NII]6583',$
               '[SII]6716','[SII]6731']

     if keyword_set(vacuum) then begin

        ; if (~ quiet) then print,'SP1_INITLINELIST: ',$
        ;                         'Using vacuum wavelengths.'

        waves = double([3727.092,3729.875,4862.69,4960.294,5008.238,$
                        6302.049,6365.536,6549.86,6564.60,6585.27,$
                        6718.29,6732.67])

     endif else begin
        
        waves = double([3726.03,3728.73,4861.32,4958.83,5006.77,$
                        6300.20,6363.67,6547.96,6562.80, 6583.34,$
                        6716.31, 6730.68])

     endelse
     

  endif

  linelist = { label:labels, wave:waves }

  return,linelist

end
