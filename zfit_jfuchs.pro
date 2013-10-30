;Written by JTF in Spring 2011 for his awesome Honors project. You
;should read it if you have not already! This program determines some
;of the inital constraints for the initialization file by doing a
;rough fitting to the H-alpha and Ca-H line. It will print to the
;xterm the fits to either the H-alpha or Ca line, depending on how you
;have it set up. This program calls in the fitting program SP1 after
;the initialization file is created. After the fit is complete it
;checks that the fit was not comepletely wrong. If it was, a value is
;changed in the initialization file and the fit is repeated. This is
;what I spent a large portion of this semester working on, so enjoy! 

function zfit,file
;pro zfit,time=time,verbose=verbose

; File to analyze
;  file = '69692'

; Image Location
  imloc = '/Users/jfuchs/Documents/ellipticals/zoo/spectra/'

; Data Output Location
  dataout = '/Users/jfuchs/Documents/ellipticals/zoo/data/'

; READ SPECTRUM
  specfile = imloc + file + 'd.fits'
  errfile  = imloc + file + 'e.fits'
  spec = readspec(specfile,header=spechead)
  wave = spec[*,0]
  flux = spec[*,1]
  specerr  = readspec(errfile)
  err = specerr[*,1]

; Determine z
  zuse = SXPAR(spechead,'Z')

; Determine Ca continuum

; P = [0d,1d]
  w1 = (1d +zuse)*3934.77d -20d
  w2 = (1d +zuse)*3934.77d -12d
  w3 = (1d +zuse)*3934.77d +12d
  w4 = (1d +zuse)*3934.77d +20d
  caindex = where((wave GT w1 AND wave LT w2) OR (wave GT w3 AND wave LT w4))
  P = POLY_FIT(wave[caindex],flux[caindex],1,measure_errors= 1d /err[caindex]^2d )


; Fit Gaussian to CaII line, determine centroid and width

  caab = where(wave GT w2 AND wave LT w3)
  cacontflux = POLY(wave[caab],P)
  castart = flux[caab] - cacontflux
  cafit = MPFITPEAK(wave[caab],castart,A,error=err[caab],/NEGATIVE,nterms=3)
  capeakflux = A[0]
  cacentroid = A[1]
  casigma = A[2]

; Calculate redshift based on Ca fitting

  caredshift = (cacentroid / 3934.77d) -1d

; Determine H-alpha continuum 

;  R = [0d,1d]
  w5 = (1+zuse)*6564.52d -50d
  w6 = (1+zuse)*6564.52d -40d
  w7 = (1+zuse)*6564.52d -10d
  w8 = (1+zuse)*6564.52d 
  w9 = (1+zuse)*6564.52d +10d
  w10= (1+zuse)*6564.52d +40d
  w11= (1+zuse)*6564.52d +50d
  haindex = where((wave GT w5 AND wave LT w6) OR (wave GT w10 AND wave LT w11))
  R = POLY_FIT(wave[haindex],flux[haindex],1,measure_errors= 1d /err[haindex]^2d )

; Fit Gaussian to H-alpha

  haem = where(wave GT w7 AND wave LT w9)
  hacontflux = POLY(wave[haem],R)
  hastart = flux[haem] - hacontflux  
  hafit = MPFITPEAK(wave[haem],hastart,B,error=err[haem],nterms=3)
  hapeakflux = B[0]
  hasigma = B[2]

; Determine Equivalent Width

  haflux = B[0]*B[2]*SQRT(2*!DPI)
  hacon = R[0] + R[1]*w8
  eqwidth = (-1d * haflux + hacon) / hacon

; Determine ~S/N for Ha line
  hadev = stddev(flux[haindex])
  haratio = hapeakflux / hadev

; Ca range for plots
  carange = where(wave GT w1 AND wave LT w4)
  capoly = POLY(wave[carange],P)
  caranged = flux[carange] - capoly

; Ha Range for plots
  harange = where(wave GT w5 AND wave LT w11)
  hapoly = POLY(wave[harange],R)
  haranged = flux[harange] - hapoly 

; Plot Calcium line fits

;  cleanplot,/silent
;  set_plot,'x'
;  device,decomposed=0
;  loadct,0,/silent
;  plot,[0],xrange=[w1,w4],yrange=[min(caranged),max(caranged)]
;  loadct,13,/silent
;  oplot,wave[carange],caranged,color=125
;  oplot,wave[caab],cafit,color=255

; Plot H-alpha line fits
  cleanplot,/silent
  set_plot,'x'
  device,decomposed=0
  loadct,0,/silent
  plot,[0],xrange=[w5,w11],yrange=[min(haranged),max(haranged)]
  loadct,13,/silent
  oplot,wave[harange],haranged,color=125
  oplot,wave[haem],hafit,color=255

;  print,zuse,caredshift,eqwidth
  initfile = sp1init(dataout,imloc,file,zuse,caredshift,eqwidth,haratio)
  sp1_fit_spectra, dataout + file + '.init'
  sp1_anl_spectra, dataout + file + '.init'


  file2 = dataout + file + '.fit.dat'
  OPENR,lun,file2, /GET_LUN
  dummy=''
  array = dblarr(4)
  READF,lun,dummy
  READF,lun,array
  FREE_LUN, lun
  fiterr = array[0]
  IF fiterr GT 15d THEN BEGIN
     initfile = sp1init(dataout,imloc,file,zuse,caredshift,eqwidth,haratio,masksig=2d)
     sp1_fit_spectra, dataout + file + '.init'
     sp1_anl_spectra, dataout + file + '.init'
  ENDIF

     
end
