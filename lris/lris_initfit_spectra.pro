function lris_initfit_spectra,lab,mask

; Directory where raw data is located
  rootindir = '/Users/drupke/metalgrads/fullspec/'
; Direcotory to search for initialization files and in which to put
; output files
  rootoutdir = '/Users/drupke/metalgrads/specfits/'
; IDL save file with stellar template structure
  startempfile = '/Users/drupke/src/idl/uhspecfit/stellar_models/'+$
                 'gonzalezdelgado/SSPGeneva_z020.sav'
;  startempfile = '/Users/drupke/src/idl/uhspecfit/stellar_models/'+$
;                 'gonzalezdelgado/SSPGeneva_z008+020+040.sav'

; Defaults
;
; Velocity dispersion used to compute initial guess for emission-line
; Gaussian sigmas; applied at mean wavelength of data.
  vdisp = 100d
; (observed) Wavelength where red and blue LRIS data are stitched
; together.
  stitchwave = 5600d
; Set blue z = red z?
  zfix = 0
; Fit strong lines only?
  stronglines = 0
; Add exponential terms to continuum fit?
  addexp = 0

  labmask = lab+mask
  infile = rootindir+labmask
  errfile = rootindir+labmask+'_err'
  outfile = rootoutdir+lab+'/'+labmask
  compfile = rootoutdir+lab+'/'+labmask+'.comp.dat'
  initfile = rootoutdir+lab+'/'+labmask+'.initfit.dat'
  initparinfo = 'lris_initparinfo'

  if labmask eq 'u12914+5b' then $
     zinit = 0.015d
  if labmask eq 'u12914+5f' then $
     zinit = 0.015d
  if labmask eq 'u813+6b1' then $
     zinit = 0.0175d
  if labmask eq 'u813+6b2' then $
     zinit = 0.0175d
  if labmask eq 'u813+6f1' then begin
     zinit = [0.017d,0.018d,0.018d,0.018d,0.017d,0.017d,0.019d]
     addexp = intarr(8)
     addexp[2] = 1
  endif
  if labmask eq 'u813+6f2' then $
     zinit = 0.0175d
  if labmask eq 'u12545+6b1' then $
     zinit = 0.02d
  if labmask eq 'u12545+6b2' then $
     zinit = 0.02d
  if labmask eq 'u12545+6f1' then $
     zinit = 0.02d
  if labmask eq 'u312b' then $
     zinit = 0.015d
  if labmask eq 'u312f' then $
     zinit = 0.015d
  if labmask eq 'n3994+5b1' then $
     zinit = 0.01d
  if labmask eq 'n3994+5b2' then $
     zinit = 0.01d
  if labmask eq 'n3994+5f1' then begin
     zinit = [dblarr(11)+0.01d,0.207]
     zfix = [intarr(11),1]
     stronglines = [intarr(11),1]
  endif
  if labmask eq 'n3994+5f2' then begin
     zinit = [dblarr(3)+0.01d,0.327,dblarr(5)+0.01d,0.075,0.01d]
     zfix = [intarr(3),1,intarr(7)]
     stronglines = [intarr(3),1,intarr(5),1,0]
  endif
  if labmask eq 'tadpoleb1' then begin
     zinit = [0.0295d,dblarr(6)+0.031d]
     zfix = [1,intarr(6)]
     stronglines = [1,intarr(6)]
  endif

; Make sure data exists  
  if infile eq '' then begin
     print,"INITFIT_SPECTRA: Mask not found."
     init = -1
  endif

  init = {zinit            : zinit, $
          zfix             : zfix, $
          stronglines      : stronglines, $
          infile           : infile, $
          errfile          : errfile, $
          outfile          : outfile, $
          compfile         : compfile, $
          initfile         : initfile, $
          initparinfo      : initparinfo, $
          vdisp            : vdisp, $
          stitchwave       : stitchwave, $
          addexp           : addexp, $
          startempfile     : startempfile }

  return,init

end
