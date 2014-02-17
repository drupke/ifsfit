; docformat = 'rst'
;
;+
;
; Extract data and information from an IFS data cube FITS file.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    IDL structure
;
; :Params:
;    infile: in, required, type=string
;      Name and path of input FITS file.
;
; :Keywords:
;    header: out, optional, type=structure
;      Headers for each extension.
;    quiet: in, optional, type=byte
;      Suppress progress messages.
;    oned: in, optional, type=byte
;      Input cube has only one non-wavelength dimension.
;    datext: in, optional, type=integer, default=1
;      Extension # of data plane.
;    varext: in, optional, type=integer, default=2
;      Extension # of variance plane.
;    dqext: in, optional, type=integer, default=3
;      Extension # of DQ plane.
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
;      2010jun08, DSNR, created as GMOS_READCUBE
;      2013dec17, DSNR, ported to IFSF_READCUBE
;      2014jan29, DSNR, added ability to change default extensions
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
function ifsf_readcube,infile,header=header,quiet=quiet,oned=oned,$
                       datext=datext,varext=varext,dqext=dqext

; Set header keywords that describe wavelength solution.
  if ~ keyword_set(oned) then begin
     crvalstr = 'CRVAL3'
     crpixstr = 'CRPIX3'
     cdeltstr = 'CD3_3'
  endif else begin
     crvalstr = 'CRVAL1'
     crpixstr = 'CRPIX1'
     cdeltstr = 'CDELT1'
  endelse

  if ~ keyword_set(quiet) then print,'IFSF_READCUBE: Loading data.'
  if ~ keyword_set(datext) then datext=1
  if ~ keyword_set(varext) then varext=2
  if ~ keyword_set(dqext) then dqext=3

; Read fits file.
  phu = readfits(infile,header_phu,ext=0,silent=quiet)
  dat = readfits(infile,header_dat,ext=datext,silent=quiet)
  var = readfits(infile,header_var,ext=varext,silent=quiet)
  dq = readfits(infile,header_dq,ext=dqext,silent=quiet)

; Get #s of rows, columns, and wavelength pixels.     
  datasize = size(dat)
  if ~ keyword_set(oned) then begin
     ncols = datasize[1]
     nrows = datasize[2]
     nz    = datasize[3]
  endif else begin
     nz    = datasize[1]
     ncols = datasize[2]
     nrows = 1
  endelse

; Create wavelength array.
  pix = dindgen(nz)
  crval = double(sxpar(header_dat,crvalstr,silent=quiet))
  crpix = double(sxpar(header_dat,crpixstr,silent=quiet))
  cdelt = double(sxpar(header_dat,cdeltstr,silent=quiet))
  wave = crval + cdelt*(pix-crpix+1) 

  cube = { $
         phu:  phu,$
         dat:  dat,$
         var:  var,$
         dq:   dq,$
         wave: wave,$
         nrows: nrows,$
         ncols: ncols,$
         nz: nz,$
         crval: crval,$
         cdelt: cdelt,$
         crpix: crpix $
         }
  
  if keyword_set(header) then $
     header = { $
              phu: header_phu,$
              dat: header_dat,$
              var: header_var,$
              dq: header_dq $
              }

  return,cube

end
