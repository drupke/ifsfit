PRO ifsf_fituvabssimp
  readcol, 'testscript', galaxyname, galaxyredshift, width0, width1, $
    format='(A,D,D,D,D,D,D,D)'

  galaxy_num = N_ELEMENTS(galaxyname)

  FOR I = 0, galaxy_num-1 DO BEGIN
    ;*****************************************
    ;This establishes variables and arrays that will be used later
    ;*****************************************
    filename= galaxyname[I]
    ext='.txt'
    linecount=FILE_LINES(filename+ext)
    c = 299792.458d
    
    readcol,filename+'n'+ext, wavelength, relativeflux
    readcol,filename+'f'+ext, useless, fit
    readcol,filename+ext, useless, useless, error
    error = error/fit

END