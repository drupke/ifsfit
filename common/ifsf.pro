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
pro ifsf,initproc,cols=cols,rows=rows,oned=oned,onefit=onefit,$
         verbose=verbose,_extra=ex
  
  starttime = systime(1)
  time = 0
  if keyword_set(verbose) then quiet=0 else quiet=1
  if keyword_set(oned) then oned=1 else oned=0

; Get fit initialization
  initdat = call_function(initproc,_extra=ex)
  
; Get linelist
  linelist = ifsf_linelist(initdat.lines)
  nlines = linelist.count()

  masksig_secondfit_def = 2d

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if not tag_exist(initdat,'datext') then datext=1 else datext=initdat.datext
  if not tag_exist(initdat,'varext') then varext=2 else varext=initdat.varext
  if not tag_exist(initdat,'dqext') then dqext=3 else dqext=initdat.dqext
  
  cube = ifsf_readcube(initdat.infile,quiet=quiet,oned=oned,$
                       datext=datext,varext=varext,dqext=dqext)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Loop through spaxels
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if ~ keyword_set(cols) then cols=[1,cube.ncols] $
  else if n_elements(cols) eq 1 then cols = [cols,cols]
  cols = fix(cols)
  for i=cols[0]-1,cols[1]-1 do begin

     print,'Column ',i+1,' of ',cube.ncols,format='(A,I0,A,I0)'

     if ~ keyword_set(rows) then rows=[1,cube.nrows] $
     else if n_elements(rows) eq 1 then rows = [rows,rows]
     rows = fix(rows)
     for j=rows[0]-1,rows[1]-1 do begin

        if oned then begin
           flux = cube.dat[*,i]
;          absolute value takes care of a few deviant points
           err = sqrt(abs(cube.var[*,i]))
           bad = cube.dq[*,i]
        endif else begin
           print,'  Row ',j+1,' of ',cube.nrows,format='(A,I0,A,I0)'
           flux = reform(cube.dat[i,j,*],cube.nz)
           err = reform(sqrt(abs(cube.var[i,j,*])),cube.nz)
           bad = reform(cube.dq[i,j,*],cube.nz)
        endelse

;       Apply DQ plane
        indx_bad = where(bad gt 0,ct)
        if ct gt 0 then begin
           flux[indx_bad] = 0d
           err[indx_bad] = max(err)*100d
        endif

        nodata = where(flux ne 0d,ct)
        if ct ne 0 then begin
           
;          Extract # of components and initial redshift guesses
;          specific to this spaxel, and write as hashes.
           ncomp = hash(initdat.lines)
           foreach line,initdat.lines do $
              if oned then ncomp[line] = (initdat.ncomp)[line,i] $
              else ncomp[line] = (initdat.ncomp)[line,i,j]

             
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; First fit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
           
;          Initialize stellar redshift for this spaxel
           if oned then zstar = initdat.zinit_stars[i] $
           else zstar = initdat.zinit_stars[i,j]
           
;          Remove NaI D line for purposes of continuum fit by maximizing
;          error. 
           if not tag_exist(initdat,'keepnad') then begin
              nadran_rest = [5850d,5900d]
              nadran = (1d + zstar) * nadran_rest
              indx_nad = where(cube.wave ge nadran[0] AND $
                               cube.wave le nadran[1],ct)
              if ct gt 0 then err[indx_nad]=max(err)
           endif

;          Option to tweak cont. fit
           if tag_exist(initdat,'tweakcntfit') then $
              tweakcntfit = reform(initdat.tweakcntfit[i,j,*,*],3,$
                                   n_elements(initdat.tweakcntfit[i,j,0,*])) $
           else tweakcntfit = 0

;          Initialize starting wavelengths
           linelistz = hash(initdat.lines)
           foreach line,initdat.lines do $
              if keyword_set(oned) then linelistz[line] = $
                 reform(linelist[line]*(1d + (initdat.zinit_gas)[line,i,*]),$
                        initdat.maxncomp) $
              else linelistz[line] = $
                 reform(linelist[line]*(1d + (initdat.zinit_gas)[line,i,j,*]),$
                        initdat.maxncomp)
                 
           structinit = ifsf_fitspec(cube.wave,flux,err,zstar,linelist,$
                                     linelistz,ncomp,initdat,quiet=quiet,$
                                     tweakcntfit=tweakcntfit)
           
           testsize = size(structinit)
           if testsize[0] eq 0 then begin
              print,'IFSF: Aborting.'
              goto,nofit
           endif
           if not quiet then print,'FIT STATUS: ',structinit.fitstatus
           if structinit.fitstatus eq -16 then begin
              print,'IFSF: Aborting.'
              goto,nofit
           endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Second fit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
           
           if not keyword_set(onefit) then begin
           
;             Set emission line mask parameters
              linepars = ifsf_sepfitpars(linelist,structinit.param,$
                                         structinit.perror,structinit.parinfo)
              linelistz = linepars.wave
              if tag_exist(initdat,'masksig_secondfit') then $
                 masksig_secondfit = initdat.masksig_secondfit $
              else masksig_secondfit = masksig_secondfit_def           
              maskwidths = hash(initdat.lines)
              foreach line,initdat.lines do $
                 maskwidths[line] = masksig_secondfit*linepars.sigma[line]

              struct = ifsf_fitspec(cube.wave,flux,err,structinit.zstar,$
                                    linelist,$
                                    linelistz,ncomp,initdat,quiet=quiet,$
                                    maskwidths=maskwidths,$
                                    peakinit=linepars.fluxpk,$
                                    siginit_gas=linepars.sigma,$
                                    tweakcntfit=tweakcntfit)           
              testsize = size(struct)
              if testsize[0] eq 0 then begin
                 print,'IFSF: Aborting.'
                 goto,nofit
              endif
              if not quiet then print,'FIT STATUS: ',struct.fitstatus
              if struct.fitstatus eq -16 then begin
                 print,'IFSF: Aborting.'
                 goto,nofit
              endif

           endif else struct = structinit

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Save result to a file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

           if oned then $
              save,struct,file=$
                   string(initdat.outdir,initdat.label,'_',i+1,'.xdr',$
                          format='(A,A,A,I04,A,I04,A)') $
           else $
              save,struct,file=$
                   string(initdat.outdir,initdat.label,'_',i+1,'_',j+1,$
                          '.xdr',format='(A,A,A,I04,A,I04,A)')

nofit:
        
        endif

     endfor

  endfor

  print,'Total time for calculation: ',systime(1)-starttime,' s.',$
        format='(/,A0,I0,A0,/)'

end
