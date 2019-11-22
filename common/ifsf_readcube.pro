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
;    datext: in, optional, type=integer, default=1
;      Extension # of data plane. Set to a negative number if the correct
;      extension is 0, since an extension of 0 ignores the keyword.
;    dqext: in, optional, type=integer, default=3
;      Extension # of DQ plane. Set to a negative number if there is no DQ; DQ
;      plane is then set to 0.
;    error: in, optional, type=byte, default=0
;      If the data cube contains errors rather than variance, the routine will
;      convert to variance.
;    header: out, optional, type=structure
;      Headers for each extension.
;    invvar: in, optional, type=byte
;      Set if the data cube holds inverse variance instead of variance. The
;      output structure will still contain the variance.
;    linearize: in, optional, type=byte
;      If set, resample the input wavelength scale so it is linearized.
;    oned: in, optional, type=byte
;      [Deprecated.] Input cube has only one non-wavelength dimension.
;    quiet: in, optional, type=byte
;      Suppress progress messages.
;    varext: in, optional, type=integer, default=2
;      Extension # of variance plane.
;    waveext: in, optional, type=integer
;      The extention number of a wavelength array.
;    zerodq: in, optional, type=byte
;      Zero out the DQ array.
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
;      2018feb08, DSNR, added WAVEEXT, INVVAR, and ZERODQ keywords
;      2018feb23, DSNR, added LINEARIZE keyword
;      2018aug12, DSNR, ensure values of DATEXT, VAREXT, DQEXT don't get changed
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
function ifsf_readcube,infile,header=header,quiet=quiet,oned=oned,$
                       datext=datext,varext=varext,dqext=dqext,$
                       vormap=vormap,error=error,waveext=waveext,$
                       invvar=invvar,zerodq=zerodq,linearize=linearize,$
                       gooddq=gooddq

  if ~ keyword_set(quiet) then print,'IFSF_READCUBE: Loading data.'

  if ~ keyword_set(datext) then begin
     datext=1
     print,'IFSF_READCUBE: Setting data extension to 1.'
     print,'IFSF_READCUBE: Set DATEXT to a negative # if it needs to be 0.'
  endif
  if datext lt 0 then datext_use=0 else datext_use=datext
  if keyword_set(varext) then begin
     if varext lt 0 then varext_use=0 else varext_use = varext
  endif else varext_use=2
  if ~ keyword_set(dqext) then dqext_use=3 else dqext_use = dqext

; Read fits file.
  if datext_use ne 0 then begin
     phu = readfits(infile,header_phu,ext=0,silent=quiet)
  endif else begin
     phu = 0b
     header_phu = ''
  endelse
  dat = readfits(infile,header_dat,ext=datext_use,silent=quiet)
  if varext_use gt 0 then begin
     var = readfits(infile,header_var,ext=varext_use,silent=quiet)
     if keyword_set(invvar) then var = 1d/var
  endif else begin
     var = dat*0d
     header_var = ''
  endelse
  if dqext_use ge 0 then begin
     dq = readfits(infile,header_dq,ext=dqext_use,silent=quiet)
     if keyword_set(gooddq) then begin
        igddq = []
        for i=0,n_elements(gooddq)-1 do begin
           igddqtmp = where(dq eq gooddq[i],ctgddq)
           if ctgddq gt 0 then igddq = [igddq,igddqtmp]
        endfor
        dq = dq*0d + 1d
        if n_elements(igddq) gt 0 then dq[igddq] = 0d
     endif
     if keyword_set(zerodq) then dq*=0d
  endif else begin
     dq = dat*0d
     header_dq = ''
  endelse
  if keyword_set(waveext) then begin
     wave = readfits(infile,header_wave,ext=waveext,silent=quiet)
;    these were chosen to match MANGA defaults ...
     crval = 0
     crpix = 1
     cdelt = 1
  endif else begin
     header_wave = ''
  endelse

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
     cdeltstropt = 'CD1_1'
  endif else if datasize[0] eq 1 then begin
     ncols = 1
     nrows = 1
     nz = datasize[1]
     wavedim = 1
     crvalstr = 'CRVAL1'
     crpixstr = 'CRPIX1'
     cdeltstr = 'CDELT1'
     cdeltstropt = 'CD1_1'
  endif

; Create wavelength array.
  if ~ keyword_set(waveext) then begin
     pix = dindgen(nz)+1d
     crval = double(sxpar(header_dat,crvalstr,silent=quiet,count=countval))
     crpix = double(sxpar(header_dat,crpixstr,silent=quiet,count=countpix))
     cdelt = double(sxpar(header_dat,cdeltstr,silent=quiet,count=countdel))
     if countdel eq 0 then $
        cdelt = double(sxpar(header_dat,cdeltstropt,silent=quiet,count=countdel))
;    Assume absence of CRPIX quantity means it is 1
     if countpix eq 0 then crpix = 1d
     if countval eq 0 OR countdel eq 0 then begin
        print,'IFSF_READCUBE: WARNING -- missing a header parameter. Wavelength'
        print,'               solution will be incorrect.'
     endif
     wave = crval + cdelt*(pix-crpix)
  endif

  if keyword_set(vormap) then begin

    ncols = max(vormap)
    nrows = 1
    vordat = dblarr(ncols,nrows,nz)
    vorvar = dblarr(ncols,nrows,nz)
    vordq = dblarr(ncols,nrows,nz)
    vorcoords = intarr(ncols,2)
    nvor = intarr(ncols)
    for i=1,ncols do begin
      ivor = where(vormap eq i,ctivor)
      xyvor = array_indices(vormap,ivor[0])
      vordat[i-1,0,*] = dat[xyvor[0],xyvor[1],*]
      vorvar[i-1,0,*] = var[xyvor[0],xyvor[1],*]
      vordq[i-1,0,*] = dq[xyvor[0],xyvor[1],*]
      vorcoords[i-1,*] = xyvor
      nvor[i-1] = ctivor
    endfor
    dat = vordat
    var = vorvar
    dq = vordq
    
  endif

  if keyword_set(error) then var=var^2d

  if keyword_set(linearize) then begin
     waveold = wave
     datold = dat
     varold = var
     dqold = dq
     crpix = 1d
     cdelt = double(wave[nz-1]-wave[0])/double(nz-1)
     wave = double(wave[0]) + dindgen(nz)*cdelt
     crval = wave[0]
     IF datasize[0] eq 3 then begin
        for i=0,ncols-1 do begin
           for j=0,nrows-1 do begin
               dat[i,j,*] = interpol(datold[i,j,*],waveold,wave,/spline)
               var[i,j,*] = interpol(varold[i,j,*],waveold,wave,/spline)
               dq[i,j,*] = interpol(dqold[i,j,*],waveold,wave)
           endfor
        endfor
     ENDIF
     IF datasize[0] eq 2 then begin
        for i=0,ncols-1 do begin
              dat[i,*] = interpol(datold[i,*],waveold,wave,/spline)
              var[i,*] = interpol(varold[i,*],waveold,wave,/spline)
              dq[i,*] = interpol(dqold[i,*],waveold,wave)
        endfor
     ENDIF
     IF datasize[0] eq 1 then begin
        dat = interpol(datold,waveold,wave,/spline)
        var = interpol(varold,waveold,wave,/spline)
        dq = interpol(dqold,waveold,wave)
     endif
     print,'IFSF_READCUBE: Interpolating DQ; values > 0.01 set to 1.'
     ibd = where(dq gt 0.01,ctbd)
     if ctbd gt 0 then dq[ibd] = 1
  endif
  
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
  if keyword_set(vormap) then $
     cube = create_struct(cube,'vorcoords',vorcoords,'nvor',nvor)
  
  if keyword_set(header) then $
     header = { $
              phu: header_phu,$
              dat: header_dat,$
              var: header_var,$
              dq: header_dq,$
              wave: header_wave $
              }

  return,cube

end
