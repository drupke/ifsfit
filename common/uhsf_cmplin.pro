; docformat = 'rst'
;
;+
;
; Compute emission line profile of a given line and velocity
; component.
;
; :Categories:
;    UHSPECFIT
;
; :Returns:
;    Array of fluxes representing emission line profile.
;
; :Params:
;    instr: in, required, type=structure
;      Contains parameters of line profile. Required tags are PARAM,
;      which is an array of best fit line parameters output from MPFIT;
;      LINELABEL, which is an array of line labels; and WAVE, which is
;      the wavelength array of the full spectrum.
;    line: in, required, type=string
;      Name of line for which to compute profile.
;    comp: in, required, type=integer
;      Number of velocity component for which to compute line profile.
;
; :Keywords:
;    velsig: in, optional, type=byte, default=0
;      Set if line sigma in PARAM array is in velocity space (km/s).
;      
; :Author:
;    David Rupke
;
; :History:
;    ChangeHistory::
;      2013sep12  DSNR  made into stand-alone routine
;      2013oct09, DSNR, added documentation
;
;-
function uhsf_cmplin,instr,line,comp,velsig=velsig

  c = 299792.458d
  
  iline = where(instr.linelabel eq line,ct)
  ppoff = instr.param[0]
  ncomp = instr.param[1]
  ppoff0 = ppoff - (ncomp-1)

  nline = n_elements(instr.linelabel)
  indices = ppoff+(comp-1)*nline*3+iline*3
  indices = indices[0] + indgen(3)
  gausspar = instr.param[indices]
  if keyword_set(velsig) then gausspar[2] *= gausspar[1]/c
  flux = gaussian(instr.wave,gausspar,/double)

  return,flux

end
