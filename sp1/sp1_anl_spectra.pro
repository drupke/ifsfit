pro sp1_anl_spectra,initfile
;
; History
;  09nov24  DSNR  created
;

  fwhmtosig = 2d*sqrt(2d*alog(2d))
  
; Threshold for delineating linear/loglinear dispersion.  Greater than
; threshold assumes spectral resolution in km/s, loglinear dispersion,
; and loglinear resolution.  Lower than/equal to threshold assumes
; spectral resolution in A, non-loglinear dispersion, and constant
; resolution in A.
  specresthresh = 10d

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
  readf,ilun,vacair,format='(A0)'
  readf,ilun,sigstr,format='(A0)'
  readf,ilun,dispstr,format='(A0)'
  readf,ilun,strongstr,format='(A0)'
  free_lun,ilun

  infile = outfile+'.xdr'

; Get linelist
  if specres gt specresthresh then velsig=1 else velsig=0
  if vacair eq 'vacuum' then vacuum=1 else vacuum=0
  if strongstr eq 'strong' then strong = 1 else strong = 0
  linelist = sp1_initlinelist(strong=strong,vacuum=vacuum,/quiet)

; Restore fit from binary file
  restore,file=infile

  if n_elements(struct.param) gt 1 then $
     ncomp = struct.param[1] $
  else $
     ncomp = 0

; Separate best-fit parameters
  if ncomp gt 0 then begin
     linepars = sepfitpars(struct.param,struct.perror)
     nlines = n_elements(linepars.flux)
  endif

; Plot results
  if ~ keyword_set(noplot) then begin
     plotstelfit,struct,outfile+'_stel',/zbuf
     if ncomp gt 0 then $
        plotstronglines,struct,outfile+'_strg',/zbuf,/ploti,$
                        velsig=velsig
  endif

; Print line fluxes to a text file
  if ncomp gt 0 then begin
     outlines=['[OII]3726','[OII]3729','Hbeta','[OIII]5007',$
               '[OI]6300','Halpha','[NII]6583','[SII]6716','[SII]6731']
     sp1_printlinepars,[''],[0],[0],outfile+'.lines.dat',/init,$
                       whichlines=outlines
     sp1_printlinepars,struct.linelabel,linepars.flux,$
                       linepars.fluxerr,$
                       outfile+'.lines.dat',/append,$
                       whichlines=outlines
  endif

; Print fit parameters to a text file
  openw,fitunit,outfile+'.fit.dat',/get_lun
  printf,fitunit,'#Rchi2','Niter','fwhm','z',$
         format='(A-6,A6,A8,A9)'
  if ncomp gt 0 then begin
     printf,fitunit,$
            struct.redchisq,$
            struct.niter,$
            format='(D6.2,I6.0,$)'
     iha = where(linelist.label eq 'Halpha')
     for i=0,ncomp-1 do begin
        printf,fitunit,fwhmtosig*linepars.sigma[iha,i],$
               struct.z.gas[i],format='(D8.2,D9.5,$)'
     endfor
  endif else printf,fitunit,1,1,-1,-1,format='(I8.0,I8.0,D8.2,D9.5,$)'

  printf,fitunit,''
  free_lun,fitunit

end
