function lris_feb09_initfit_spectra,lab,mask

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
; Order for fitting red part of data for continuum fit
  redord_def = 5
  redord = redord_def

  labmask = lab+mask
  infile = rootindir+labmask
  errfile = rootindir+labmask+'_err'
  outfile = rootoutdir+lab+'/'+labmask
  compfile = rootoutdir+lab+'/'+labmask+'.comp.dat'
  initfile = rootoutdir+lab+'/'+labmask+'.initfit.dat'
  initparinfo = 'lris_feb09_initparinfo'
  fcnlinefit = 'lris_feb09_manygauss'
  fcncontfit = 'lris_fit_continuum'
  if labmask eq 'n2207m1' then begin
     zinit = dblarr(19)+0.009d
     zfix = intarr(19)
     zfix[13]=1
  endif
  if labmask eq 'n2207m2' then begin
     zinit = 0.009d
  endif
  if labmask eq 'n3994+5m1' then begin
     zinit = dblarr(14)+0.01d
  endif
  if labmask eq 'n3994+5m2' then begin
     zinit = dblarr(12)+0.01d
     redord = dblarr(12)+redord_def
     redord[10] = 11
     redord[11] = 11
  endif
  if labmask eq 'n4485+90m1' then begin
     zinit = dblarr(17)+0.0016d
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
          fcnlinefit       : fcnlinefit, $
          fcncontfit       : fcncontfit, $
          redord           : redord, $
          vdisp            : vdisp, $
          stitchwave       : stitchwave, $
          startempfile     : startempfile }

  return,init

end
