;Written by JTF in Spring 2011 for his awesome Honors project. This
;file sets up and writes the initilization file for SP1. It uses data
;calculated by zfit.pro to write the file and determine the initial
;guesses.  

function sp1init,dataout,imloc,file,zuse,caredshift,eqwidth,haratio,masksig=masksig

  OpenW,lun,dataout + file + '.init',/get_lun

; FITS file containing data spectrum
  m1 = imloc + file + 'd.fits'

; FITS file containing error spectrum
  m2 = imloc + file + 'e.fits'

; Prefix for output files
  m3 = dataout + file

; Lower and upper limits for fit, in A
  m4 = '3700 6800'

; Redshift estimate
  m5 = string(zuse,format='(D0.5)')

; List of component velocities
  IF haratio LT 7d THEN m6 = '-1' ELSE m6='0'

; Redshift of stellar z
  m7 = string(caredshift,format='(D0.5)')

; Estimate of galaxy velocity dispersion (sigma)
  m8 = '90'

; Spectral resolution (sigma)
  m9 = '66'

; Masksig
  m95 = '3'
  IF keyword_set(masksig) then m95 = string(masksig,format='(D0.1)')

; Spectrum is in 'vacuum' or 'air' wavelengths
  m10 = 'vacuum'

; 'fix' emission line widths to Halpha?
  m11='fix'

; 'smooth' stellar spectrum
  m12 = 'smooth'

; Fit only 'strong' or 'all' emission lines?
  m13 = 'strong'

; Using stellar z? If yes, 'z'
  m14 = 'z'


  Printf,lun,m1,m2,m3,m4,m5,m6,m7,m8,m9,m95,m10,m11,m12,m13,m14,format='(A-0)'
  Free_lun,lun


end
