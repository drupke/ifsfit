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
;      [Deprecated.] Input cube has only one non-wavelength dimension.
;    datext: in, optional, type=integer, default=1
;      Extension # of data plane. Set to a negative number if the correct 
;      extension is 0, since an extension of 0 ignores the keyword.
;    varext: in, optional, type=integer, default=2
;      Extension # of variance plane.
;    dqext: in, optional, type=integer, default=3
;      Extension # of DQ plane. Set to a negative number if there is no DQ; DQ 
;      plane is then set to 0.
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
;      2014aug05, DSNR, small tweak to allow single spectra and no DQ
;      2015may20, DSNR, updated logic in reading # of rows, cols, and wavelength
;                       points to be more flexible; added wavedim to output
;                       structure
;      2016sep12, DSNR, fixed DATEXT logic so can specify an extension of 0;
;                       added warning when correct wavelength keywords not found;
;                       added second dispersion keyword option.
;    
; :Copyright:
;    Copyright (C) 2013--2016 David S. N. Rupke
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

  if ~ keyword_set(quiet) then print,'IFSF_READCUBE: Loading data.'

  if ~ keyword_set(datext) then begin
     datext=1
     print,'IFSF_READCUBE: Setting data extension to 1.'
     print,'IFSF_READCUBE: Set DATEXT to a negative # if it needs to be 0.'
  endif
  if datext lt 0 then datext=0
  if ~ keyword_set(varext) then varext=2
  if ~ keyword_set(dqext) then dqext=3

; Read fits file.
  if datext ne 0 then begin
     phu = readfits(infile,header_phu,ext=0,silent=quiet)
  endif else begin
     phu = 0b
     header_phu = ''
  endelse
  dat = readfits(infile,header_dat,ext=datext,silent=quiet)
  var = readfits(infile,header_var,ext=varext,silent=quiet)
  if dqext ge 0 then $
     dq = readfits(infile,header_dq,ext=dqext,silent=quiet) $
  else dq = dat*0d

; Get #s of rows, columns, and wavelength pixels. Logic for 2-d case assumes 
; that # wavelength pts > # cols.
  datasize = size(dat)
  IF datasize[0] eq 3 then begin
     ncols = datasize[1]
     nrows = datasize[2]
     nz = datasize[3]
     wavedim = 3
     crvalstr = 'CRVAL3'
     crpixstr = 'CRPIX3'
     cdeltstr = 'CD3_3'
     cdeltstropt = 'CDELT3'
  endif else if datasize[0] eq 2 then begin
     ncols = (datasize[1] gt datasize[2]) ? datasize[2] : datasize[1]
     nz = (datasize[1] gt datasize[2]) ? datasize[1] : datasize[2]
     wavedim = (datasize[1] gt datasize[2]) ? 1 : 2
     nrows = 1
     crvalstr = 'CRVAL1'
     crpixstr = 'CRPIX1'
     cdeltstr = 'CDELT1'
  endif else if datasize[0] eq 1 then begin
     ncols = 1
     nrows = 1
     nz = datasize[1]
     wavedim = 1
     crvalstr = 'CRVAL1'
     crpixstr = 'CRPIX1'
     cdeltstr = 'CDELT1'
  endif

; Create wavelength array.
  pix = dindgen(nz)
  crval = double(sxpar(header_dat,crvalstr,silent=quiet,count=countval))
  crpix = double(sxpar(header_dat,crpixstr,silent=quiet,count=countpix))
  cdelt = double(sxpar(header_dat,cdeltstr,silent=quiet,count=countdel))
  if countdel eq 0 then $
     cdelt = double(sxpar(header_dat,cdeltstropt,silent=quiet,count=countdel))
  if countval eq 0 OR countpix eq 0 OR countdel eq 0 then begin
     print,'IFSF_READCUBE: WARNING -- missing a header parameter. Wavelength'
     print,'               solution will be incorrect.'
  endif
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
         wavedim: wavedim,$
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
