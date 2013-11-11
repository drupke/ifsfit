; docformat = 'rst'
;
;+
;
; Interpolate templates from template wavelength grid to data
; wavelength grid.
;
; :Categories:
;    UHSPECFIT
;
; :Returns:
;    The interpolated templates, of type dblarr(nwave_spec, ntemplates).
;
; :Params:
;    spec_lam: in, required, type=dblarr(nwave_spec)
;      Wavelengths of data arrays.
;    temp_lam: in, required, type=dblarr(nwave_temp)
;      Wavelengths of template arrays.
;    template: in, required, type=dblarr(nwave_temp\, ntemplates)
;      Model fluxes from templates.
;
; :Keywords:
;
; :Author:
;    Jabran Zahid and David Rupke
;
; :History:
;    Change History::
;      2009, HJZ, created
;      2013oct17, DSNR, documented
;
;-
function uhsf_interptemp, spec_lam, temp_lam, template, $
                          temp_lam_rest=temp_lam_rest

  ss = size(template)
  new_temp = fltarr(n_elements(spec_lam), ss[2])

  if min(temp_lam) gt min(spec_lam) then $
     print, '  INTERPOL_TEMPLATE: WARNING: Extrapolating template from ',$
            min(temp_lam),' to ',min(spec_lam),'.',format='(A,I0,A,I0,A)'
  if max(temp_lam) lt max(spec_lam) then $
     print, '  INTERPOL_TEMPLATE: WARNING: Extrapolating template from ',$
            max(temp_lam),' to ',max(spec_lam),'.',format='(A,I0,A,I0,A)'

  for i = 0, ss[2] - 1 do new_temp[*, i] = $
     interpol(template[*, i],temp_lam,spec_lam)

  return, new_temp

end
