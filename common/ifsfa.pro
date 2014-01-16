; docformat = 'rst'
;
;+
;
; This procedure is the core routine to plot the continuum and emission
; lines fits to a spectrum.
;
; As input, it requires a structure of initialization parameters and the output 
; XDR file from IFSF. The tags for the initialization structure can be found in 
; INITTAGS.txt.
;
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Various plots and data files.
;
; :Params:
;    initproc: in, required, type=string
;      Name of procedure to initialize the fit.
;
; :Keywords:
;    cols: in, optional, type=intarr, default=all
;      Columns to fit, in 1-offset format. Either a scalar or a
;      two-element vector listing the first and last columns to fit.
;    rows: in, optional, type=intarr, default=all
;      Rows to fit, in 1-offset format. Either a scalar or a
;      two-element vector listing the first and last rows to fit.
;    noplots: in, optional, type=byte
;      Disable plotting.
;    oned: in, optional, type=byte
;      Data is assumed to be in a 2d array; choose this switch to
;      input data as a 1d array.
;    verbose: in, optional, type=byte
;      Print error and progress messages. Propagates to most/all
;      subroutines.
; 
; :Author:
;    David Rupke
;
; :History:
;    Change History::
;      2009may13, DSNR, created
;      2013oct04, DSNR, started re-write for new data
;      2013oct09, DSNR, documented
;      2013nov25, DSNR, renamed, added copyright and license; changed
;                       required parameters from 'gal' and 'bin' to
;                       'initproc', and optional parameter 'fibers' to
;                       'oned', to make it more general
;      2014jan13, DSNR, propagated use of hashes; 
;                       updated for new linelist routine; 
;                       re-wrote IFSF_PRINTLINPAR and created IFSF_PRINTFITPAR
;      2014jan16, DSNR, bugfixes
;
;-
pro ifsfa,initproc,cols=cols,rows=rows,noplots=noplots,oned=oned,verbose=verbose

  fwhmtosig = 2d*sqrt(2d*alog(2d))
  if keyword_set(oned) then oned=1 else oned=0

; Get fit initialization
  initdat=call_function(initproc)
  
; Get linelist
; Get linelist
  linelist = ifsf_linelist(initdat.lines)
  nlines = linelist.count()

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  cube = ifsf_readcube(initdat.infile,quiet=quiet,oned=oned)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Initialize output files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if not tag_exist(initdat,'outlines') then outlines = linelist->keys() $
  else outlines = initdat.outlines
  ifsf_printlinpar,outlines,linlun,$
                   outfile=initdat.outdir+initdat.label+'.lin.dat'
  ifsf_printfitpar,fitlun,$
                   outfile=initdat.outdir+initdat.label+'.fit.dat'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Loop through spaxels
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if not keyword_set(cols) then cols=[1,cube.ncols] $
  else if n_elements(cols) eq 1 then cols = [cols,cols]
  cols = fix(cols)
  for i=cols[0]-1,cols[1]-1 do begin

     if keyword_set(verbose) then $
        print,'Column ',i+1,' of ',cube.ncols,format='(A,I0,A,I0)'

     if not keyword_set(rows) then rows=[1,cube.nrows] $
     else if n_elements(rows) eq 1 then rows = [rows,rows]
     rows = fix(rows)
     for j=rows[0]-1,rows[1]-1 do begin

        if oned then begin
           flux = cube.dat[*,i]
           err = sqrt(abs(cube.var[*,i]))
           lab = string(i+1,format='(I04)')
        endif else begin
           if keyword_set(verbose) then $
              print,'  Row ',j+1,' of ',cube.nrows,format='(A,I0,A,I0)'
           flux = reform(cube.dat[i,j,*],cube.nz)
           err = reform(sqrt(abs(cube.var[i,j,*])),cube.nz)
           lab = string(i+1,'_',j+1,format='(I04,A,I04)')
        endelse

;       Restore fit, after a couple of sanity checks
;       TODO: Mark zero-ed spectra as bad earlier!
        infile = initdat.outdir+initdat.label+'_'+lab+'.xdr'
        outfile = initdat.outdir+initdat.label+'_'+lab
        nodata = where(flux ne 0d,ct)
        if file_test(infile) OR ct gt 0 then restore,file=infile $
        else begin
           print,'IFSFA: No XDR file for ',i+1,', ',j+1,$
                 ', or data is everywhere 0.',$
                 format='(A0,I4,A0,I4,A0)'
           goto,nofit
        endelse

;       Restore original error.
;       TODO: Is this necessary?
        struct.spec_err = err[struct.gd_indx]

;       Get line fit parameters
        linepars = ifsf_sepfitpars(linelist,struct.param,struct.perror,$
                                   struct.parinfo)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plot
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        if not keyword_set(noplots) then begin

;          Plot continuum
           if tag_exist(initdat,'fcnpltcont') then $
              fcnpltcont=initdat.fcnpltcont $
           else fcnpltcont='ifsf_pltcont'
           call_procedure,fcnpltcont,struct,outfile+'_cnt'
;          Plot emission lines
           if not linepars.nolines then begin
              if tag_exist(initdat,'fcnpltlin') then $
                 fcnpltlin=initdat.fcnpltlin else fcnpltlin='ifsf_pltlin'
              if tag_exist(initdat,'argspltlin1') then $
                 call_procedure,fcnpltlin,struct,initdat.argspltlin1,$
                                outfile+'_lin1'
              if tag_exist(initdat,'argspltlin2') then $
                 call_procedure,fcnpltlin,struct,initdat.argspltlin2,$
                                outfile+'_lin2'
           endif
              
;          Print fit parameters to a text file
           ifsf_printfitpar,fitlun,i+1,j+1,struct

;          Print line fluxes and Halpha Weq to a text file
           if not linepars.nolines then $ 
              ifsf_printlinpar,outlines,linlun,i+1,j+1,initdat.maxncomp,linepars


        endif

nofit:

     endfor

  endfor

  free_lun,fitlun
  free_lun,linlun

end
