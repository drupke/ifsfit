pro sp1_fit_spectra,initfile,time=time,verbose=verbose

;
; History
;  09nov24  DSNR  copied from GMOS routine
;
; NOTE: Limit to speed comes in line fitting function (fcnlinefit);
; change to 'manygauss' over 'manygauss_slow' if constant dispersion
; is available.
;

  starttime = systime(1)
  if ~ keyword_set(time) then time=0
  if keyword_set(verbose) then quiet=0 else quiet=1
  
; Threshold for delineating linear/loglinear dispersion.  Greater than
; threshold assumes spectral resolution in km/s, loglinear dispersion,
; and loglinear resolution.  Lower than/equal to threshold assumes
; spectral resolution in A, non-loglinear dispersion, and constant
; resolution in A.
  specresthresh = 10d

; Upper limit for emission line masking width, in A.
  masklim = 25d

; Read input file
  specfile=''
  errfile=''
  outfile=''
  fitranstr=''
  mcompstr=''
  vacair=''
  sigstr=''
  dispstr=''
  strongstr=''
  zstelstr=''

  openu,ilun,initfile,/get_lun
  readf,ilun,specfile,format='(A0)'
  readf,ilun,errfile,format='(A0)'
  readf,ilun,outfile,format='(A0)'
  readf,ilun,fitranstr,format='(A0)'
  readf,ilun,zguess,format='(D)'
  readf,ilun,mcompstr,format='(A0)'
  readf,ilun,dzstel,format='(D)'
  readf,ilun,vdisp,format='(D)'
  readf,ilun,specres,format='(D)'
  readf,ilun,masksig,format='(D)'
  readf,ilun,vacair,format='(A0)'
  readf,ilun,sigstr,format='(A0)'
  readf,ilun,dispstr,format='(A0)'
  readf,ilun,strongstr,format='(A0)'
  if ~ eof(ilun) then readf,ilun,zstelstr,format='(A0)'
  free_lun,ilun

  fitran_rest=double(strsplit(fitranstr,/extract))
  mcomp=double(strsplit(mcompstr,/extract))

  if specres gt specresthresh then loglam=1 else loglam=0
  if vacair eq 'vacuum' then vacuum=1 else vacuum=0
  if dispstr eq 'smooth' then disperse = 1 else disperse=0
  if strongstr eq 'strong' then strong = 1 else strong = 0
  if zstelstr eq 'z' then zstel = dzstel else zstel = -1

; Fitting routines
  fcninitpar='sp1_initparinfo'
  fcnz = 'sp1_redshift_spec'
  fcnlinefit='manygauss_slow'
  fcncontfit='fit_continuum'

; Stellar model library
  startempfile = '/Users/drupke/src/idl/uhspecfit/stellar_models/'+$
                 'gonzalezdelgado/SSPGeneva_z020.sav'

; Get linelist
  linelist = sp1_initlinelist(strong=strong,vacuum=vacuum,quiet=quiet)
  nlines = n_elements(linelist.wave)

; Read spectrum
  spec = readspec(specfile)
  wave = spec[*,0]
  flux = spec[*,1]
  specerr  = readspec(errfile)
  err = specerr[*,1]

; Initial set of redshifts
  if mcomp[0] eq -1 then ncomp = 0 else ncomp = n_elements(mcomp)
  if ncomp gt 0 then begin
     zguessold = zguess
     zguess = refinez(zguess,wave/(1+zguess),flux)
     if ~ quiet then $
        print,'Gas redshift adjusted from ',zguessold,' to ',zguess,$
              format='(A0,D0.5,A0,D0.5)'
     dzcomp = mcomp/299792d
     zgas = zguess + dzcomp
  endif else $
     zgas = zguess
  if zstel eq -1 then zstar = zguess+dzstel else zstar = dzstel
  zstarold = zstar
  zstar = refinestz(zstar,wave,flux,err)
  if ~ quiet then $
     print,'Stellar redshift adjusted from ',zstarold,' to ',zstar,$
           format='(A0,D0.5,A0,D0.5)'
  z = {star : zstar,$
       gas  : zgas}

; Fix emission lines to 0?  This has to be hardwired presently.
  if ncomp gt 0 then $
     zerolines = dblarr(nlines,ncomp) $
  else zerolines = 0

; First fit
  if fcnlinefit eq 'manygauss_slow' AND loglam then begin
     if ~ quiet then print,'SP1_FIT_SPECTRA: Assuming loglinear resolution.'
     argslinefit = {velsig:specres}
  endif

  argsinitpar = {sigma: sigstr,$
                 specres: specres,$
                 specrthresh: specresthresh,$
                 zerolines: zerolines}
  fitran = sp1_redshift_spec(fitran_rest,z)

  structinit = fit_spectrum(wave,flux,err,startempfile,z,ncomp,/subtract,$
                            vdisp=vdisp,linelist=linelist,quiet=quiet,$
                            time=time,fitran=fitran,obj_id=outfile,$
                            fcninitpar=fcninitpar,argsinitpar=argsinitpar,$
                            fcnlinefit=fcnlinefit,argslinefit=argslinefit,$
                            fcnz=fcnz,fcncontfit=fcncontfit,$
                            loglam=loglam,disperse=disperse,vacuum=vacuum)

; Repeat fit with new z estimate and better estimate for line widths,
; if we're fitting emission lines.

  if ncomp gt 0 then begin

;    Update de-redshifted wavelengths and z array
     sp1_updatez,structinit.param,z
     if zstel eq -1 then z.star = z.gas[0] + dzstel

;    Update linelist so that maskwidths array works OK
     newlinelist = {wave:structinit.linewave,label:structinit.linelabel}

;    Estimate linewidths.  Mask at MASKSIG # of sigma away from line
;    in either direction
;    Make emission-line masks based on best-fit parameters
     linepars = sepfitpars(structinit.param,structinit.perror)
     masklines = reform(linepars.wave[*,0:ncomp-1],$
                        n_elements(linepars.wave[*,0:ncomp-1]))
     maskwidths = masksig*reform(linepars.sigma[*,0:ncomp-1],$
                                 n_elements(linepars.sigma[*,0:ncomp-1]))
     if loglam then maskwidths *= masklines / 299792d
;    Limit masking to a certain width
     toowide = where(maskwidths ge masklim,ct)
     if ct gt 0 then maskwidths[toowide] = masklim
;    Don't mask certain lines/components
;     maskwidths[indgen(7)+nlines] = 0d
;    Redshift fit range
     fitran = sp1_redshift_spec(fitran_rest,z)
 
     struct = fit_spectrum(wave,flux,err,startempfile,z,ncomp,$
                           vdisp=vdisp,/subtract,linelist=newlinelist,$
                           quiet=quiet,obj_id=lab,time=time,fitran=fitran,$
                           fcninitpar=fcninitpar,argsinitpar=argsinitpar,$
                           maskwidths=maskwidths,masklines=masklines,$
                           fcnlinefit=fcnlinefit,argslinefit=argslinefit,$
                           fcncontfit=fcncontfit,sigguess=linepars.sigma,$
                           peakguess=linepars.fluxpk,loglam=loglam,$
                           fcnz=fcnz,disperse=disperse,vacuum=vacuum)
     
;    Update z array and fitted line wavelengths, and add to structure
     sp1_updatez,struct.param,z

  endif else begin

     struct = structinit

  endelse


; Save result to binary file
  struct = jjadd_tag(struct,'z',z)
  struct = jjadd_tag(struct,'fitrange',fitran,/array)
  save,struct,file=outfile+'.xdr'
  
  if (~ quiet) then print,'Total time for calculation: ',$
                          systime(1)-starttime,' s.',$
                          format='(/,A0,I0,A0,/)'
  
end
