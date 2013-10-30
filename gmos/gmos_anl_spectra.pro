pro gmos_anl_spectra,gal,bin,rows=rows,cols=cols,$
                     fibers=fibers,noplots=noplots,$
                     verbose=verbose
;
; History
;  09may13  DSNR  created
;  09jun07  DSNR  added multiple components
;  10may27  DSNR  started re-write for new data
;

  fwhmtosig = 2d*sqrt(2d*alog(2d))
  if keyword_set(fibers) then fibers=1 else fibers=0

; Get fit initialization
  initdat = gmos_initfit_spectra(gal,bin=bin)
  ncompinit = initdat.ncomp
  infile = initdat.infile
  ctoutfile = initdat.ctoutfile
  qsoutfile = initdat.qsoutfile
  qsoutfile_ha = initdat.qsoutfile_ha
  outdir = initdat.outdir
  plotstelfit = initdat.plotstelfit
  qsocntargs = initdat.qsocntargs
  fcnstrlines = initdat.fcnstrlines
  argstrlines = initdat.argstrlines
  nadfitord = initdat.nadfitord
  restcomp = initdat.nadrestcomp
  zbuf = 1
  argslinelist = initdat.argslinelist
  outlines = initdat.outlines
  
; Get linelist
  if keyword_set(argslinelist) then $
     linelist = call_function('gmos_initlinelist',_extra=argslinelist) $
  else linelist = gmos_initlinelist()
; which lines to output
  if ~keyword_set(outlines) then outlines=['[OI]6300','Halpha','[NII]6583']
     

  if fibers then begin
;    Read data
     data = readfits(infile,header,ext=2,/silent)
     var = readfits(infile,ext=3,/silent)
     dq = readfits(infile,ext=4,/silent)
     datasize = size(data)
     nz    = datasize[1]
     ncols = datasize[2]
     nrows = 1
  endif else begin
     data = readfits(infile,header,ext=1,/silent)
     var = readfits(infile,ext=2,/silent)
     dq = readfits(infile,ext=3,/silent)
     datasize = size(data)
     ncols = datasize[1]
     nrows = datasize[2]
     nz    = datasize[3]
  endelse

; Initialize output files
; Output text files
  gmos_printlinepars,[''],[0],[0],outdir+gal+'.lines.dat',-1,-1,/init,$
                     whichlines=outlines,/weq
  gmos_printweqha,0,outdir+gal+'.weqha.dat',-1,-1,/init
  openw,fitunit,outdir+gal+'.fit.dat',/get_lun
  printf,fitunit,'#Col','Row','Cmp','Rchi2','Niter','FWHM','z',$
         format='(A-4,2A4,2A6,A8,A10)'
; Output data cubes
  if ~ fibers then begin
     ctoutcube = dblarr(ncols,nrows,nz)
     qsoutcube = dblarr(ncols,nrows,nz)
     qsoutcube_ha = dblarr(ncols,nrows,nz)
  endif

; Cycle through spectra
  if ~ keyword_set(cols) then cols=[1,ncols] $
  else if n_elements(cols) eq 1 then cols = [cols,cols]
  cols = fix(cols)
  firstspectrum=0
  for i=cols[0]-1,cols[1]-1 do begin

     if keyword_set(verbose) then $
        print,'Column ',i+1,' of ',ncols,format='(A,I0,A,I0)'

     if ~ keyword_set(rows) then rows=[1,nrows] $
     else if n_elements(rows) eq 1 then rows = [rows,rows]
     rows = fix(rows)
     for j=rows[0]-1,rows[1]-1 do begin

;       Load raw data
        if fibers then begin
           flux = data[*,i]
           err = sqrt(abs(var[*,i]))
           lab = string(i+1,format='(I04)')
        endif else begin
           if keyword_set(verbose) then $
              print,'  Row ',j+1,' of ',nrows,format='(A,I0,A,I0)'
           lab = string(i+1,'_',j+1,format='(I04,A,I04)')
           flux = data[i,j,*]
           err = sqrt(abs(var[i,j,*]))
        endelse
        nodata = where(flux ne 0d,ct)
        if ct ne 0 then begin

           if (i eq 1 AND j eq 1) then append=0 else append=1

;          Restore fit
           infile = outdir+gal+'_'+lab+'.xdr'
           outfile = outdir+gal+'_'+lab
           if file_test(infile) then restore,file=infile $
           else begin
              print,'GMOS_ANL_SPECTRA: Spectrum ',i+1,', ',j+1,$
                    ' does not exist.',$
                    format='(A0,I4,A0,I4,A0)'
              goto,nofit
           endelse

;          Initialize output data cube         
           if ~ firstspectrum AND keyword_set(qsocntargs) then begin

              waveout = struct.wave[50:n_elements(struct.wave)-51]
              istart = where(struct.wave eq waveout[0])
              nzout = n_elements(struct.wave) - 100

              haz = gmos_redshift_spec([6563d],struct.z)
              istart_haz = value_locate(struct.wave,haz - 30d)
              iend_haz = value_locate(struct.wave,haz + 30d)
              nzout_haz = iend_haz-istart_haz+1

              ctoutcube = dblarr(ncols,nrows,nzout)
              qsoutcube = dblarr(ncols,nrows,nzout)
              qsoutcube_ha = dblarr(ncols,nrows,nzout_haz)

              firstspectrum = 1

           endif

;          Restore original error
           struct.spec_err = err[struct.gd_indx]
           
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

;          Case with QSO:
           if keyword_set(qsocntargs) then begin

              if ~ keyword_set(noplots) then begin
                 call_procedure,plotstelfit,struct,outfile+'_stel',$
                                zbuf=zbuf,qsocntargs=qsocntargs
                 if ncomp gt 0 then $
                    gmos_plotstronglines,struct,outfile+'_strg',zbuf=1,$
                                         /plot 
              endif

;             Compute pure starburst spectrum
              ctspec = gmos_subtractqso(struct,qsocntargs=qsocntargs,$
                                        qsomod=qsomod)
              istart = where(struct.wave eq waveout[0])
              haz = gmos_redshift_spec([6563d],struct.z)
              istart_haz = value_locate(struct.wave,haz - 30d)
              ctoutcube[i,j,*] = ctspec[istart:istart+nzout-1]
              qsoutcube[i,j,*] = qsomod[istart:istart+nzout-1]
              qsoutcube_ha[i,j,*] = qsomod[istart_haz:istart_haz+nzout_haz-1]

;             Print pure starburst spectrum to a text file
              gmos_printnadspec,struct,outfile+'_nad_spec.dat',$
                                nadfitord,wavenorm,fluxnorm,parnorm,weq,$
                                qsocntargs=qsocntargs,restcomp=restcomp
              gmos_plotnad,struct,outfile+'_nad',$
                           wavenorm,fluxnorm,parnorm,zbuf=zbuf,$
                           qsocntargs=qsocntargs

;          Case w/o QSO:
           endif else begin
              if ~ fibers then begin
                 if ~ keyword_set(noplots) then begin
                    plotstelfit,struct,outfile+'_stel',zbuf=zbuf,$
                                outstelfit=outstelfit
                 endif
                 gmos_printnadspec,struct,outfile+'_nad_spec.dat',$
                                   nadfitord,wavenorm,fluxnorm,parnorm,weq,$
                                   restcomp=restcomp
                 if ~ keyword_set(noplots) then $
                    gmos_plotnad,struct,outfile+'_nad',$
                                 wavenorm,fluxnorm,parnorm,zbuf=zbuf
              endif
              if ncomp gt 0 AND ~ keyword_set(noplots) then begin
                 if keyword_set(argstrlines) then $
                    call_procedure,fcnstrlines,struct,outfile+'_strg',$
                                   zbuf=zbuf,$
                                   /plot,_extra=argstrlines $
                 else $
                    call_procedure,fcnstrlines,struct,outfile+'_strg',$
                                   zbuf=zbuf,/plot
              endif
           endelse

;          Print fit parameters to a text file
           printf,fitunit,i+1,j+1,c1,struct.redchisq,struct.niter,$
                  fwhm_c1,z_c1,$
                  format='(I4,I4,I4,D6.2,I6,D8.2,D10.6)'
           for k = 1,ncomp-1 do $
              printf,fitunit,i+1,j+1,k+1,-1,-1,$
                     fwhmtosig*linepars.sigma[0,k],struct.z.gas[k],$
                     format='(I4,I4,I4,2I6,D8.2,D10.6)'

;         Print line fluxes and Halpha Weq to a text file
           if ncomp then begin

              gmos_printlinepars,struct.linelabel,linepars.flux,$
                                 linepars.fluxerr,$
                                 outdir+gal+'.lines.dat',i+1,j+1,/append,$
                                 whichlines=outlines,weq=weq
              gmos_printweqha,struct,outdir+gal+'.weqha.dat',i+1,j+1

           endif

        endif

nofit:

     endfor

  endfor

; Write output cube
  if (~ fibers AND keyword_set(qsocntargs)) then begin
     header_out = header
     sxaddpar,header_out,'CRPIX3',1
     sxaddpar,header_out,'CRVAL3',waveout[0]
     writefits,ctoutfile,ctoutcube,header_out
     writefits,qsoutfile,qsoutcube,header_out
     writefits,qsoutfile_ha,qsoutcube_ha,header_out
  endif

  free_lun,fitunit

end
