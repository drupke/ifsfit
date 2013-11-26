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
;    
; :Copyright:
;    Copyright (C) 2013 David S. N. Rupke
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
pro ifsf,initproc,cols=cols,rows=rows,oned=oned,$
         verbose=verbose
  
  starttime = systime(1)
  time = 0
  if keyword_set(verbose) then quiet=0 else quiet=1
  if keyword_set(oned) then oned=1 else oned=0

; Get fit initialization
  initdat = call_function(initproc)
  
; Get linelist
  if tag_exist(initdat,'argslinelist') then linelist = $
     call_function('ifsf_linelist',initdat.lines,_extra=initdat.argslinelist) $
  else linelist = ifsf_linelist(initdat.lines)
  nlines = linelist.count()

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  if oned then begin
     data = readfits(initdat.infile,header,ext=2,/silent)
     var = readfits(initdat.infile,ext=3,/silent)
     dq = readfits(initdat.infile,ext=4,/silent)
     datasize = size(data)
     nz    = datasize[1]
     ncols = datasize[2]
     nrows = 1
;    Wavelength solution
     wave = dindgen(nz)
     crval = double(sxpar(header,'CRVAL1',/silent))
     crpix = double(sxpar(header,'CRPIX1',/silent))
     cdelt = double(sxpar(header,'CDELT1',/silent))
     wave = crval + cdelt*(wave-crpix+1) 
  endif else begin
     data = readfits(initdat.infile,header,ext=1,/silent)
     var = readfits(initdat.infile,ext=2,/silent)
     dq = readfits(initdat.infile,ext=3,/silent)
     datasize = size(data)
     ncols = datasize[1]
     nrows = datasize[2]
     nz    = datasize[3]
     wave = dindgen(nz)
     crval = double(sxpar(header,'CRVAL3',/silent))
     crpix = double(sxpar(header,'CRPIX3',/silent))
;     cdelt = double(sxpar(header,'CDELT3',/silent))
     cdelt = double(sxpar(header,'CD3_3',/silent))
     wave = crval + cdelt*(wave-crpix+1) 
  endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Loop through spaxels
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if ~ keyword_set(cols) then cols=[1,ncols] $
  else if n_elements(cols) eq 1 then cols = [cols,cols]
  cols = fix(cols)
  for i=cols[0]-1,cols[1]-1 do begin

     print,'Column ',i+1,' of ',ncols,format='(A,I0,A,I0)'

     if ~ keyword_set(rows) then rows=[1,nrows] $
     else if n_elements(rows) eq 1 then rows = [rows,rows]
     rows = fix(rows)
     for j=rows[0]-1,rows[1]-1 do begin

        if oned then begin
           flux = data[*,i]
;          absolute value takes care of a few deviant points
           err = sqrt(abs(var[*,i]))
           bad = dq[*,i]
        endif else begin
           print,'  Row ',j+1,' of ',nrows,format='(A,I0,A,I0)'
           flux = reform(data[i,j,*],nz)
           err = reform(sqrt(abs(var[i,j,*])),nz)
           bad = reform(dq[i,j,*],nz)
        endelse

;       Apply DQ plane
        indx_bad = where(bad gt 0,ct)
        if ct gt 0 then begin
           flux[indx_bad] = 0d
           err[indx_bad] = max(err)*100d
        endif

        nodata = where(flux ne 0d,ct)
        if ct ne 0 then begin

; Turn hashes into spaxel-specific arrays

           ncomp = dblarr(nlines)
           linewave = dblarr(nlines)
           linetie = dblarr(nlines)
           zinit_gas = dblarr(nlines,initdat.maxncomp)
           for k=0,n_elements(initdat.lines)-1 do begin
              ncomp[k] = initdat.ncomp[initdat.lines[k],i,j]
              linewave[k] = linelist[initdat.lines[k]]
              linetie[k] = initdat.linetie[initdat.lines[k]]
              zinit_gas[k,*] = initdat.zinit_gas[lines[k],i,j,*]
           endfor


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; First fit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
           
;          Initialize redshift structure
           z = {star: initdat.zinit_stars[i,j],$
                gas: zinit_gas}
           
;          Remove NaI D line for purposes of continuum fit by maximizing
;          error. 
           if not tag_exist(initdat,'keepnad') then begin
              nadran_rest = [5850d,5900d]
              nadran = (1d + z.star) * nadran_rest
              indx_nad = where(wave ge nadran[0] AND wave le nadran[1],ct)
              if ct gt 0 then err[indx_nad]=max(err)
           endif

;          Initialize starting wavelengths
           linewavez = rebin(linewave,nlines,initdat.maxncomp) * (1d + z.gas)
           
           structinit = ifsf_fitspec(wave,flux,err,z,lines,linewave,linewavez,$
                                     linetie,ncomp,initdat,quiet=quiet)

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
           
;          Set emission line mask parameters
           linepars = ifsf_sepfitpars(structinit.param,structinit.perror)
           nlines_tot = nlines*initdat.maxncomp
           linewavez = linepars.wave
           masksig = 2d
           newmaskwidths = masksig*reform(linepars.sigma,nlines_tot)

;          Put fit initialization parameters into structure
           initdat_use = initdat
           if tag_exist(initdat,'maskwidths') then initdat_use.maskwidths = $
              newmaskwidths $
           else initdat_use = $
              jjadd_tag(initdat_use,'maskwidths',newmaskwidths,/array_tag)
           if tag_exist(initdat,'peakinit') then initdat_use.peakinit = $
              linepars.fluxpk $
           else initdat_use = $
              jjadd_tag(initdat_use,'peakinit',linepars.fluxpk,/array_tag)
           if structinit.sig_stars gt 0 then $
              initdat_use.siginit_stars = structinit.sig_stars
           if tag_exist(initdat_use,'fcnoptstelz') then initdat_use = $
              rem_tag(initdat_use,'fcnoptstelz')
           if tag_exist(initdat_use,'fcnoptstelsig') then initdat_use = $
              rem_tag(initdat_use,'fcnoptstelsig')
           if tag_exist(initdat_use,'sigfitvals') then initdat_use = $
              rem_tag(initdat_use,'sigfitvals')

           struct = ifsf_fitspec(wave,flux,err,z,lines,linewave,linewavez,$
                                 linetie,ncomp,initdat_use,quiet=quiet)
           
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Save result to a file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

           struct = jjadd_tag(struct,'zstar',z.star)

           if oned then $
              save,struct,file=$
                   string(initdat.outdir,initdat.gal,'_',i+1,'.xdr',$
                          format='(A,A,A,I04,A,I04,A)') $
           else $
              save,struct,file=$
                   string(initdat.outdir,initdat.gal,'_',i+1,'_',j+1,$
                          '.xdr',format='(A,A,A,I04,A,I04,A)')

nofit:
        
        endif

     endfor

  endfor

  print,'Total time for calculation: ',systime(1)-starttime,' s.',$
        format='(/,A0,I0,A0,/)'

end
