;
; History
;  09xxxYY  Jabran Zahid  created
;  10may04  DSNR          added rest-frame template option
;

function interpol_template, spec_lam, temp_lam, template, $
                            temp_lam_rest=temp_lam_rest

;spec_lam - lambda of spectra  [n_spec]
;temp_lam - lambda of template [n_temp]
;template - model flux from templates [n_temp, nn]
;new_temp - interpolated template [n_spec, nn]

  ss = size(template)

  if keyword_set(temp_lam_rest) then $
     new_temp = fltarr(n_elements(spec_lam), ss[2]*2) $
  else $
     new_temp = fltarr(n_elements(spec_lam), ss[2])

  if min(temp_lam) gt min(spec_lam) then $
     print, '  INTERPOL_TEMPLATE: WARNING: Extrapolating template from ',$
            min(temp_lam),' to ',min(spec_lam),'.',format='(A,I0,A,I0,A)'
  if max(temp_lam) lt max(spec_lam) then $
     print, '  INTERPOL_TEMPLATE: WARNING: Extrapolating template from ',$
            max(temp_lam),' to ',max(spec_lam),'.',format='(A,I0,A,I0,A)'

  for i = 0, ss[2] - 1 do begin
     new_temp[*, i] = interpol(template[*, i],temp_lam,spec_lam)
     if keyword_set(temp_lam_rest) then $
        new_temp[*,i+ss[2]] = interpol(template[*,i],temp_lam_rest,$
                                       spec_lam)
  endfor

  return, new_temp

end
