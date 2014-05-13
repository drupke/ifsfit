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
;      2014jan29, DSNR, added _extra parameter to permit passing parameters
;                       to initialization routine; added some lines to deal
;                       properly with case of 1d data "cube"
;      2014feb26, DSNR, replaced ordered hashes with hashes
;      2014apr15, DSNR, bugfix in file sanity check
;
; :Copyright:
;    Copyright (C) 2013-2014 David S. N. Rupke
;
;    This program is free software: you can redistribute it and/or
;    modify it under the terms of the GNU General Public License as
;    published by the Free Software Foundation, either version 3 of
;    the License or any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;    General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program.  If not, see
;    http://www.gnu.org/licenses/.
;
;-
pro ifsfa,initproc,cols=cols,rows=rows,noplots=noplots,oned=oned,$
          verbose=verbose,_extra=ex

  bad = 1d99
  fwhmtosig = 2d*sqrt(2d*alog(2d))
  if keyword_set(verbose) then quiet=0 else quiet=1
  if keyword_set(oned) then oned=1 else oned=0

; Get fit initialization
  initdat = call_function(initproc,_extra=ex)
  
; Get linelist
  linelist = ifsf_linelist(initdat.lines)
  nlines = linelist.count()

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if not tag_exist(initdat,'datext') then datext=1 else datext=initdat.datext
  if not tag_exist(initdat,'varext') then varext=2 else varext=initdat.varext
  if not tag_exist(initdat,'dqext') then dqext=3 else dqext=initdat.dqext

  cube = ifsf_readcube(initdat.infile,quiet=quiet,oned=oned,$
                       datext=datext,varext=varext,dqext=dqext)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Initialize output files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if not tag_exist(initdat,'outlines') then $
     outlines = (linelist->keys())->toarray() $
  else outlines = initdat.outlines
  ifsf_printlinpar,outlines,linlun,$
                   outfile=initdat.outdir+initdat.label+'.lin.dat'
  ifsf_printfitpar,fitlun,$
                   outfile=initdat.outdir+initdat.label+'.fit.dat'
                   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Initialize line hash
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  linmaps = hash()
  foreach line,outlines do $
    linmaps[line] = dblarr(cube.ncols,cube.nrows,initdat.maxncomp,4) + bad

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Loop through spaxels
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Switch to track when first NaD normalization done
  firstnadnorm=1

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
        if ~ file_test(infile) then begin
           print,'IFSFA: No XDR file for ',i+1,', ',j+1,'.',$
           format='(A0,I4,A0,I4,A0)'
           goto,nofit
        endif
        nodata = where(flux ne 0d,ct)
        if ct le 0 then begin
           print,'IFSFA: Data is everywhere 0 for ',i+1,', ',j+1,'.',$
           format='(A0,I4,A0,I4,A0)'
           goto,nofit
        endif
        restore,file=infile

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

        endif 
              
;       Print fit parameters to a text file
        ifsf_printfitpar,fitlun,i+1,j+1,struct

        foreach line,outlines do $
           linmaps[line,i,j,*,*] = [[linepars.flux[line,*]],$
                                   [linepars.fluxerr[line,*]],$
                                   [linepars.wave[line,*]],$
                                   [linepars.sigma[line,*]]]

;       Print line fluxes and Halpha Weq to a text file
        if not linepars.nolines then $ 
           ifsf_printlinpar,outlines,linlun,i+1,j+1,initdat.maxncomp,linepars

;       Normalize data around NaD and plot
        if tag_exist(initdat,'donad') then begin
           if tag_exist(initdat,'argsnormnad') then $
              normnad = ifsf_normnad(struct.wave,$
                                     struct.cont_dat/struct.cont_fit,$
                                     struct.spec_err,struct.zstar,fitpars,$
                                     _extra=initdat.argsnormnad) else $
              normnad = ifsf_normnad(struct.wave,$
                                     struct.cont_dat/struct.cont_fit,$
                                     struct.spec_err,struct.zstar,fitpars)
           if not keyword_set(noplots) then $
              if tag_exist(initdat,'argspltnormnad') then $
                 ifsf_pltnaddat,normnad,fitpars,struct.zstar,$
                    outfile+'_nad_norm',$
                    _extra=initdat.argspltnormnad else $
                 ifsf_pltnaddat,normnad,fitpars,struct.zstar,$
                    outfile+'_nad_norm'
           if firstnadnorm then begin
              nadcube = $
                 {wave: dblarr(cube.ncols,cube.nrows,n_elements(normnad[*,0])),$
                  dat: dblarr(cube.ncols,cube.nrows,n_elements(normnad[*,0])),$
                  err: dblarr(cube.ncols,cube.nrows,n_elements(normnad[*,0]))}
              firstnadnorm = 0
           endif
           nadcube.wave[i,j,*] = normnad[*,0]
           nadcube.dat[i,j,*] = normnad[*,2]
           nadcube.err[i,j,*] = normnad[*,3]
        endif

nofit:

     endfor

  endfor

  free_lun,fitlun
  free_lun,linlun

  save,linmaps,file=initdat.outdir+initdat.label+'.lin.xdr'
  save,nadcube,file=initdat.outdir+initdat.label+'.nadspec.xdr'

end
