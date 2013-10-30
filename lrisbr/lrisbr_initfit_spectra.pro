function lrisbr_initfit_spectra,lab,mask

; Directory where raw data is located
  rootindir = '/Users/drupke/metalgrads/fullspec/'
; Direcotory to search for initialization files and in which to put
; output files
  rootoutdir = '/Users/drupke/metalgrads/specfits/'
; IDL save file with stellar template structure
  startempfile = '/Users/drupke/src/idl/uhspecfit/stellar_models/'+$
                 'gonzalezdelgado/SSPGeneva_z020.sav'

; Defaults
;
; Velocity dispersion used to compute initial guess for emission-line
; Gaussian sigmas; applied at mean wavelength of data.
  disperseB = 1
  disperseR = 1
  vdispB = 150d
  vdispR = 100d
; Fit strong lines only?
  stronglines = 0
; Set tilt = 0?
  zfixB = 0
  zfixR = 1
; 
  dzstel = 0d

  labmask = lab+mask
  infile = rootindir+labmask
  errfile = rootindir+labmask+'_err'
  outfile = rootoutdir+lab+'/'+labmask
  compfile = rootoutdir+lab+'/'+labmask+'.comp.dat'
  initfile = rootoutdir+lab+'/'+labmask+'.initfit.dat'
  initparinfo = 'lrisbr_initparinfo'
  linefit = 'lrisbr_manygauss'
  if labmask eq 'arp248m1' then begin
     zinit = 0.018d
     vdispB = 250d
     vdispR = 150d
  endif
  if labmask eq 'arp248m2' then begin
     zinit = 0.018d
     vdispB = 250d
     vdispR = 150d
  endif
  if labmask eq 'arp256m1' then begin
     zinit = 0.027d
     dzstel = 0.0003d
     stronglines = intarr(10)
     stronglines[6] = 1
  endif
  if labmask eq 'arp256m2' then begin
     zinit = 0.027d
     dzstel = 0.0003d
  endif
  if labmask eq 'arp256m3' then begin
     zinit = 0.027d
     dzstel = 0.0003d
  endif
  if labmask eq 'arp298m1' then begin
     zinit = 0.016d
     dzstel = 0.0003d
  endif
  if labmask eq 'arp298m2' then begin
     zinit = 0.016d
     dzstel = 0.0003d
     zfixR = indgen(10)+1
     zfixR[1] = 0
  endif

; Make sure data exists  
  if infile eq '' then begin
     print,"LRISBR_INITFIT_SPECTRA: Mask not found."
     init = -1
  endif

  init = {zinit            : zinit, $
          zfixB            : zfixB, $
          zfixR            : zfixR, $
          dzstel           : dzstel, $
          stronglines      : stronglines, $
          infile           : infile, $
          errfile          : errfile, $
          outfile          : outfile, $
          compfile         : compfile, $
          initfile         : initfile, $
          initparinfo      : initparinfo, $ 
          linefit          : linefit, $
          disperseB        : disperseB, $
          disperseR        : disperseR, $
          vdispB           : vdispB, $
          vdispR           : vdispR, $
          startempfile     : startempfile }

  return,init

end
