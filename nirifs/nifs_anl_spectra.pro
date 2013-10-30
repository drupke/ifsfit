pro nifs_anl_spectra,gal,bin,cols=cols,rows=rows,$
                     noplots=noplots,verbose=verbose
;
; History
;  13mar05  DSNR  created
;

  fwhmtosig = 2d*sqrt(2d*alog(2d))

; Get fit initialization
  initdat = nifs_initfit_spectra(gal,bin)
  ncompinit = initdat.ncomp
  infile = initdat.infile
  ctoutfile = initdat.ctoutfile
  outdir = initdat.outdir
  plotstelfit = initdat.plotstelfit
  qsocntargs = initdat.qsocntargs
  fcnstrlines = initdat.fcnstrlines
  argstrlines = initdat.argstrlines
  zbuf = 1
  argslinelist = initdat.argslinelist
  outlines = initdat.outlines
  doubleline = initdat.doubleline
  
; which lines to output
  if ~keyword_set(outlines) then outlines=['']

  if gal eq 'mrk231' then begin
     if bin eq 999 then nocube=1 else nocube=0
     cube = nifs_readcube(infile,header=header,nocube=nocube)
  endif
  if gal eq 'f08572nw' then cube = osiris_readcube(infile,header=header)
  data = cube.dat
  wave = cube.wave
  if gal eq 'f08572nw' then wave = wave*10d
  nz = cube.nz
  ncols = cube.ncols
  nrows = cube.nrows

; Initialize output files
; Output text files
  nifs_printlinepars,[''],[0],[0],outdir+gal+'.lines.dat',-1,-1,/init,$
                     whichlines=outlines
  openw,fitunit,outdir+gal+'.fit.dat',/get_lun
  printf,fitunit,'#Col','Row','Cmp','Rchi2 ','Niter','FWHM','z',$
         format='(A-4,2A4,A7,A5,A8,A10)'

; Output continuum data cubes
  if gal eq 'f08572nw' then ctoutcube = dblarr(nz,ncols,nrows)
  if gal eq 'mrk231' then ctoutcube = dblarr(ncols,nrows,nz)

; Cycle through spectra
  if ~ keyword_set(cols) then cols=[1,ncols] $
  else if n_elements(cols) eq 1 then cols = [cols,cols]
  cols = fix(cols)
  for i=cols[0]-1,cols[1]-1 do begin

     if keyword_set(verbose) then $
        print,'Column ',i+1,' of ',ncols,format='(A,I0,A,I0)'
     
     if ~ keyword_set(rows) then rows=[1,nrows] $
     else if n_elements(rows) eq 1 then rows = [rows,rows]
     rows = fix(rows)
     for j=rows[0]-1,rows[1]-1 do begin

        if keyword_set(verbose) then $
           print,'  Row ',j+1,' of ',nrows,format='(A,I0,A,I0)'
        lab = string(i+1,'_',j+1,format='(I04,A,I04)')

;       Load raw data, to check if there was no data
        if gal eq 'mrk231' then $
           if bin eq 999 then flux=data else flux = data[i,j,*]
        if gal eq 'f08572nw' then flux = data[*,i,j]
        nodata = where(flux ne 0d,ct)

        if ct ne 0 then begin

           if (i eq 1 AND j eq 1) then append=0 else append=1

; Get linelist
           if keyword_set(argslinelist) then $
              linelist = call_function('nifs_initlinelist',$
                                       _extra=argslinelist[i,j]) $
           else linelist = nifs_initlinelist()

;          Restore fit
           infile = outdir+gal+'_'+lab+'.xdr'
           outfile = outdir+gal+'_'+lab
           if file_test(infile) then begin
              restore,file=infile
           endif else begin
              print,'NIFS_ANL_SPECTRA: Spectrum ',i+1,', ',j+1,$
                    ' does not exist.',$
                    format='(A0,I4,A0,I4,A0)'
              goto,nofit
           endelse

           izstart = value_locate(wave,struct.wave[0])
           nzfit = n_elements(struct.wave)

;          Write pure continuum spectrum to file
           
           if gal eq 'f08572nw' then $
              ctoutcube[izstart:izstart+nzfit-1,i,j] = $
              struct.spec - struct.specfit
           if gal eq 'mrk231' then ctoutcube[i,j,izstart:izstart+nzfit-1] = $
              struct.spec - struct.specfit
           
;          # of components
           if n_elements(struct.param) gt 1 then $
              ncomp = struct.param[1] $
           else $
              ncomp = 0

;          Separate best-fit parameters
           if ncomp gt 0 then begin
              linepars = sepfitpars(struct.param,struct.perror)
              nlines = n_elements(linepars.flux)
              fwhm_c1 = fwhmtosig*linepars.sigma[0,0]
              z_c1 = struct.z.gas[0]
              c1 = 1
           endif else begin 
              fwhm_c1 = -1
              z_c1 = struct.z.star
              c1 = 0
           endelse
             
           if ~ keyword_set(noplots) then begin
              if keyword_set(qsocntargs) then $
                 call_procedure,plotstelfit,struct,outfile+'_stel',$
                                zbuf=zbuf,qsocntargs=qsocntargs $
              else $
                 call_procedure,plotstelfit,struct,outfile+'_stel',$
                                zbuf=zbuf
              if ncomp gt 0 then $
                 nifs_plotstronglines,struct,outfile+'_strg',zbuf=1,$
                                      /plot 
           endif

;          Print fit parameters to a text file
           if struct.redchisq gt 1000 then redchisq = '>1000' $
           else redchisq = string(struct.redchisq,'(D7.2)')
           printf,fitunit,i+1,j+1,c1,redchisq,struct.niter,$
                  fwhm_c1,z_c1,$
                  format='(I4,I4,I4,A7,I5,D8.2,D10.6)'
           for k = 1,ncomp-1 do begin
              cmp = k+1
              if doubleline[i,j] eq 'paa' AND k eq 2 then cmp=4
              printf,fitunit,i+1,j+1,cmp,-1,-1,$
                     fwhmtosig*linepars.sigma[0,k],struct.z.gas[k],$
                     format='(I4,I4,I4,I7,I5,D8.2,D10.6)'
           endfor

;         Print line fluxes to a text file
           if ncomp then begin

              nifs_printlinepars,struct.linelabel,linepars.flux,$
                                 linepars.fluxerr,$
                                 outdir+gal+'.lines.dat',i+1,j+1,/append,$
                                 whichlines=outlines,$
                                 doubleline=doubleline[i,j]
              
           endif

        endif

nofit:

     endfor

  endfor

; Write output cubes
  writefits,ctoutfile,ctoutcube,header.dat

  free_lun,fitunit

end
