; docformat = 'rst'
;
;+
;
; This procedure is the core routine to plot the continuum and emission
; lines fits to a spectrum.
;
; As input, it requires a structure of initializaiton parameters. The
; tags for this structure can be found in ...
;
;
; :Categories:
;    UHSPECFIT/GMOS
;
; :Returns:
;    IDL save file (.xdr)
;
; :Params:
;    gal: in, required, type=string
;    bin: in, required, type=double
;
; :Keywords:
;    rows: in, optional, type=dblarr
;    cols: in, optional, type=dblarr
;    noplots: in, optional, type=byte
;    sky: in, optional, type=byte
;    verbose: in, optional, type=byte
; 
; :Author:
;    David Rupke
;
; :History:
;    Change History::
;      2009may13, DSNR, created
;      2013oct04, DSNR, started re-write for new data
;      2013oct09, DSNR, documented
;
;
;-
pro uhsf_gm_anlspec,gal,bin,cols=cols,rows=rows,$
                    fibers=fibers,noplots=noplots,$
                    sky=sky,verbose=verbose

  fwhmtosig = 2d*sqrt(2d*alog(2d))
  if keyword_set(fibers) then fibers=1 else fibers=0

; Get fit initialization
  if keyword_set(sky) then initdat = $
     call_function('uhsf_gm_init_'+gal,bin,/sky) $
  else initdat = $
     call_function('uhsf_gm_init_'+gal,bin)
  
; Get linelist
  if tag_exist(initdat,'argslinelist') then linelist = $
     call_function('uhsf_gm_linelist',_extra=initdat.argslinelist) $
  else linelist = uhsf_gm_linelist()
  nlines = n_elements(linelist.label)

  if fibers then begin
;    Read data
     data = readfits(initdat.infile,header,ext=2,/silent)
     var = readfits(initdat.infile,ext=3,/silent)
     dq = readfits(initdat.infile,ext=4,/silent)
     datasize = size(data)
     nz    = datasize[1]
     ncols = datasize[2]
     nrows = 1
  endif else begin
     data = readfits(initdat.infile,header,ext=1,/silent)
     var = readfits(initdat.infile,ext=2,/silent)
     dq = readfits(initdat.infile,ext=3,/silent)
     datasize = size(data)
     ncols = datasize[1]
     nrows = datasize[2]
     nz    = datasize[3]
  endelse

; Initialize output files
  if not tag_exist(initdat,'outlines') then outlines = linelist.label $
  else outlines = initdat.outlines
  uhsf_printlinpar,[''],[0],[0],initdat.outdir+gal+'.lines.dat',-1,-1,/init,$
                   whichlines=outlines
  openw,fitunit,initdat.outdir+gal+'.fit.dat',/get_lun
  printf,fitunit,'#Col','Row','Cmp','Rchi2','Niter','FWHM','z',$
         format='(A-4,2A4,2A6,A8,A10)'

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
           flux = reform(data[i,j,*],nz)
           err = reform(sqrt(abs(var[i,j,*])),nz)
        endelse
        nodata = where(flux ne 0d,ct)
        if ct ne 0 then begin

           if (i eq 1 AND j eq 1) then append=0 else append=1

;          Restore fit
           infile = initdat.outdir+gal+'_'+lab+'.xdr'
           outfile = initdat.outdir+gal+'_'+lab
           if file_test(infile) then restore,file=infile $
           else begin
              print,'UHSF_GM_ANLSPEC: Spectrum ',i+1,', ',j+1,$
                    ' does not exist.',$
                    format='(A0,I4,A0,I4,A0)'
              goto,nofit
           endelse

;          Restore original error
           struct.spec_err = err[struct.gd_indx]
           
;; ;          # of components
;;            if n_elements(struct.param) gt 1 then begin
;;               ncomp = reform(initdat.ncomp[i,j,*],nlines)
;;               linepars = uhsf_sepfitpars(struct.param,struct.perror)
;;            endif

           if not fibers then begin
              if ~ keyword_set(noplots) then begin
                 if tag_exist(initdat,'fcnpltcont') then $
                    fcnpltcont=initdat.fcnpltcont else $
                       fcnpltcont='uhsf_pltcont'
                 call_procedure,fcnpltcont,struct,outfile+'_cnt'
              endif
              if n_elements(struct.param) gt 1 AND $
                 ~ keyword_set(noplots) then begin
                 if tag_exist(initdat,'fcnpltlin') then $
                    fcnpltlin=initdat.fcnpltlin else $
                       fcnpltlin='uhsf_pltlin'
                 if tag_exist(initdat,'argspltlin1') then $
                    call_procedure,fcnpltlin,struct,initdat.argspltlin1,$
                                   outfile+'_lin1',/velsig
                 if tag_exist(initdat,'argspltlin2') then $
                    call_procedure,fcnpltlin,struct,initdat.argspltlin2,$
                                   outfile+'_lin2',/velsig
              endif
           endif
              
;; ;          Print fit parameters to a text file
;;            printf,fitunit,i+1,j+1,c1,struct.redchisq,struct.niter,$
;;                   fwhm_c1,z_c1,$
;;                   format='(I4,I4,I4,D6.2,I6,D8.2,D10.6)'
;;            for k = 1,ncomp-1 do $
;;               printf,fitunit,i+1,j+1,k+1,-1,-1,$
;;                      fwhmtosig*linepars.sigma[0,k],struct.z.gas[k],$
;;                      format='(I4,I4,I4,2I6,D8.2,D10.6)'

;; ;         Print line fluxes and Halpha Weq to a text file
;;            if ncomp then begin
              
;;               gmos_printlinepars,struct.linelabel,linepars.flux,$
;;                                  linepars.fluxerr,$
;;                                  outdir+gal+'.lines.dat',i+1,j+1,/append,$
;;                                  whichlines=outlines,weq=weq

;;            endif

        endif

nofit:

     endfor

  endfor

  free_lun,fitunit

end
