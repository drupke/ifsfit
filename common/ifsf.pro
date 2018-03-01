; docformat = 'rst'
;
;+
;
; This procedure is the core routine to fit the continuum and emission
; lines of a spectrum.
;
; As input, it requires a structure of initialization parameters. The
; tags for this structure can be found in INITTAGS.txt.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    IDL save file (.xdr)
;
; :Params:
;    initproc: in, required, type=string
;      Name of procedure to initialize the fit.
;
; :Keywords:
;    cols: in, optional, type=intarr, default=all
;      Columns to fit, in 1-offset format. Either a scalar or a
;      two-element vector listing the first and last columns to fit.
;    ncores: in, optional, type=int, default=1
;      Number of cores to split processing over.
;    rows: in, optional, type=intarr, default=all
;      Rows to fit, in 1-offset format. Either a scalar or a
;      two-element vector listing the first and last rows to fit.
;    oned: in, optional, type=byte
;      Data is assumed to be in a 2d array; choose this switch to
;      input data as a 1d array.
;    onefit: in, optional, type=byte
;      Option to skip second fit; primarily for testing.
;    verbose: in, optional, type=byte
;      Print error and progress messages. Propagates to most/all
;      subroutines.
; 
; :Author:
;    David S. N. Rupke::
;      Rhodes College
;      Department of Physics
;      2000 N. Parkway
;      Memphis, TN 38104
;      drupke@gmail.com
;
; :History:
;    ChangeHistory::
;      2009jul08, DSNR, copied from LRIS routine to GMOS
;      2010may27, DSNR, started re-write for new data
;      2013oct04, DSNR, started re-write for new data
;      2013nov25, DSNR, renamed, added copyright and license; changed
;                       required parameters from 'gal' and 'bin' to
;                       'initproc', and optional parameter 'fibers' to
;                       'oned', to make it more general
;      2013nov26, DSNR, added code to turn input hashes into arrays
;                       for each spaxel
;      2013dec10, DSNR, testing with PPXF and bug fixes
;      2013dec17, DSNR, inserted new IFSF_READCUBE function in place of code 
;                       block to read data cube; started propagation of hashes
;                       through code and implementation of new calling sequence
;                       rubric for IFSF_FITSPEC
;      2013dec19, DSNR, more progress on propagating use of hashes
;                       through code, through first fit
;      2014jan13, DSNR, finished propagating use of hashes
;      2014jan16, DSNR, updated treatment of redshifts; bugfixes
;      2014jan17, DSNR, bugfixes; implemented SIGINIT_GAS, TWEAKCNTFIT keywords
;      2014jan29, DSNR, added _extra parameter to permit passing parameters
;                       to initialization routine; added some lines to deal
;                       properly with case of 1d data "cube"
;      2014feb26, DSNR, replaced ordered hashes with hashes
;      2014may08, DSNR, added ability to check components automagically.
;      2016jan06, DSNR, allow no emission line fit with initdat.noemlinfit
;      2016feb12, DSNR, changed treatment of sigma limits for emission lines
;                       so that they can be specified on a pixel-by-pixel basis
;      2016sep13, DSNR, added internal logic to check if emission-line fit present
;      2016sep16, DSNR, changed logic of ONEFIT so component checking doesn't
;                       happen
;      2016sep17, DSNR, removed GOTO statements; moved separate row/col loops
;                       to a single loop; moved loop commands to IFSF_FITLOOP
;                       to enable multicore processing
;      2016sep20, DSNR, implemented multicore processing
;      2018feb08, DSNR, updated call to IFSF_READCUBE
;    
; :Copyright:
;    Copyright (C) 2013--2018 David S. N. Rupke
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
pro ifsf,initproc,cols=cols,rows=rows,oned=oned,onefit=onefit,ncores=ncores,$
         verbose=verbose,_extra=_extra
  
  starttime = systime(1)
  time = 0
  if ~ keyword_set(ncores) then ncores=1
  if ncores gt !CPU.HW_NCPU -1 then $
     message,'Number of processes to spawn greater than number of CPUs available.'
  if keyword_set(oned) then oned=1b else oned=0b
  if keyword_set(onefit) then onefit=1b else onefit=0b
  if keyword_set(verbose) then quiet=0b else quiet=1b
; Get fit initialization
  if keyword_set(_extra) then initdat = call_function(initproc,_extra=_extra) $
  else initdat = call_function(initproc)
  
; Get linelist
  if tag_exist(initdat,'lines') then begin
     if tag_exist(initdat,'waveunit') then $
        linelist = ifsf_linelist(initdat.lines,waveunit=initdat.waveunit) $
     else linelist = ifsf_linelist(initdat.lines)
     nlines = linelist.count()
  endif else begin
     linelist = hash()
  endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if not tag_exist(initdat,'datext') then datext=1 else datext=initdat.datext
  if not tag_exist(initdat,'varext') then varext=2 else varext=initdat.varext
  if not tag_exist(initdat,'dqext') then dqext=3 else dqext=initdat.dqext
  
  if tag_exist(initdat,'vormap') then vormap=initdat.vormap else vormap=0b
  if tag_exist(initdat,'argsreadcube') then $
     cube = ifsf_readcube(initdat.infile,quiet=quiet,oned=oned,$
                          datext=datext,varext=varext,dqext=dqext,$
                          vormap=vormap,_extra=initdat.argsreadcube) $
  else $
     cube = ifsf_readcube(initdat.infile,quiet=quiet,oned=oned,$
                          datext=datext,varext=varext,dqext=dqext,$
                          vormap=vormap)
                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Loop through spaxels
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; case of specifying a column and row in voronoi format
  if keyword_set(cols) AND keyword_set(rows) $
     AND tag_exist(initdat,'vormap') then begin
     
     if n_elements(cols) eq 1 AND n_elements(rows) eq 1 then begin
        cols = vormap[cols[0]-1,rows[0]-1]
        rows = 1        
     endif else begin
        message,'Can only specify 1 spaxel at a time, or all spaxels, in Voronoi mode.'
     endelse
  endif

  if ~ keyword_set(cols) then begin
    cols=[1,cube.ncols]
    ncols = cube.ncols
  endif else if n_elements(cols) eq 1 then begin
    cols = [cols,cols]
    ncols = 1
  endif else begin
    ncols = cols[1]-cols[0]+1
  endelse
  cols = fix(cols)
  if ~ keyword_set(rows) then begin
    rows=[1,cube.nrows]
    nrows = cube.nrows
  endif else if n_elements(rows) eq 1 then begin
    rows = [rows,rows]
    nrows = 1
  endif else begin
    nrows = rows[1]-rows[0]+1
  endelse
  rows = fix(rows)

  colarr = rebin(dindgen(ncols)+cols[0]-1d,ncols,nrows)
  rowarr = rebin(transpose(dindgen(nrows)+rows[0]-1d),ncols,nrows)
  nspax = ncols*nrows
  
  dolog=0b
  if ncores eq 1 then begin
     if tag_exist(initdat,'logfile') then begin
        dolog=1b
        logloop=replicate(initdat.logfile,nspax)
     endif
     for ispax=0,nspax-1 do begin
        if ispax eq 0 AND dolog then file_delete,initdat.logfile,/allow
        ifsf_fitloop,ispax,colarr,rowarr,cube,initdat,linelist,$
                     oned,onefit,quiet,logfile=logloop
     endfor
  endif else begin
     if nspax lt ncores then begin
        message,'nspax < ncores; setting ncores = nspax',/cont
        ncores=nspax
     endif
     if tag_exist(initdat,'logfile') then $
        logfiles=replicate(initdat.logfile+'_',ncores)+$
                 string(indgen(ncores)+1,format='(I0)') $
     else message,'Logfile must be specified.'
     logloop = strarr(nspax)
     min_per_core=floor(double(nspax)/double(ncores))
     num_per_core=intarr(ncores)+min_per_core
     remainder=round((double(nspax)/double(ncores)-double(min_per_core))*$
               double(ncores))
     for i=0,ncores-1 do begin
        if i+1 le remainder then num_per_core[i]++
        logloop[i+indgen(num_per_core[i])*ncores] = $
           replicate(logfiles[i],num_per_core[i])
     endfor
     drt_bridgeloop,nspax,ncores,initdat.batchfile,initdat.batchdir,$
                    invar=['colarr','rowarr','cube','initdat',$
                           'linelist','oned','onefit','quiet','logloop'],$
                    loopvar='ispax',logfiles=logfiles
  endelse
  
  print,'Total time for calculation: ',systime(1)-starttime,' s.',$
        format='(/,A0,I0,A0,/)'

end
