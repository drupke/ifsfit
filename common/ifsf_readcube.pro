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
;    error: in, optional, type=bool, default=0
;      If the data cube contains errors rather than variance, the routine will
;      convert to variance.
;    fluxnorm: in, optional, type=double
;      Divide data and variance by this value.
;    header: out, optional, type=structure
;      Headers for each extension.
;    invvar: in, optional, type=bool
;      Set if the data cube holds inverse variance instead of variance. The
;      output structure will still contain the variance.
;    linearize: in, optional, type=bool
;      If set, resample the input wavelength scale so it is linearized.
;    oned: in, optional, type=bool
;      [Deprecated.] Input cube has only one non-wavelength dimension.
;    quiet: in, optional, type=bool
;      Suppress progress messages.
;    varext: in, optional, type=integer, default=2
;      Extension # of variance plane.
;    waveext: in, optional, type=integer
;      The extention number of a wavelength array.
;    zerodq: in, optional, type=bool
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
;      2020may05, DSNR, new treatment of default axes in 2D images; added CUNIT
;                       and BUNIT to output
;      2021jan04, DSNR, added NDIM to output; fixed bug in linearization for 
;        ndim=2 case (reversed indices)
;      2022jan27, DSNR, added BUNIT_VAR to output
;      2022jul13, DSNR, added FLUXNORM keyword; if CUNIT is 'nm', convert to A
;      2022jul28, DSNR, if var or dq is absent, set either to -1 (output of
;         readfits)
;    
; :Copyright:
;    Copyright (C) 2013--2022 David S. N. Rupke
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
                       gooddq=gooddq,fluxnorm=fluxnorm

  if ~ keyword_set(quiet) then print,'IFSF_READCUBE: Loading data.'

  if ~ keyword_set(datext) then begin
     datext=1
     if ~ keyword_set(quiet) then begin
        print,'IFSF_READCUBE: Setting data extension to 1.'
        print,'IFSF_READCUBE: Set DATEXT to a negative # if it needs to be 0.'
     endif
  endif
  if datext lt 0 then datext_use=0 else datext_use=datext
  if keyword_set(varext) then begin
     if varext lt 0 then varext_use=0 else varext_use = varext
  endif else varext_use=2
  if ~ keyword_set(dqext) then dqext_use=3 else dqext_use = dqext

; Read fits file.
  if datext_use ne 0 then begin
     phu = mrdfits(infile,0,header_phu,silent=quiet)
  endif else begin
     phu = 0b
     header_phu = ''
  endelse
  dat = mrdfits(infile,datext_use,header_dat,silent=quiet)
  if varext_use gt 0 then begin
     novar=0b
     var = mrdfits(infile,varext_use,header_var,silent=quiet)
     if n_elements(var) eq 1 then begin
        if var[0] eq -1 then begin
           novar=1b
           header_var = ''
        endif
     endif
     if keyword_set(invvar) and ~novar then var = 1d/var
  endif else begin
     var = dat*0d
     ; in case any of dat has a nan
     inan = where(finite(dat,/nan),ctnan)
     if ctnan gt 0 then var[inan] = 0d
     header_var = header_dat
     sxaddpar,header_var,'EXTNAME','VAR',silent=quiet
     if ~ keyword_set(quiet) then begin
        print,'IFSR_READCUBE: No variance extension.'
        print,'   Setting var=0*dat and variance header to data header.'
     endif
  endelse
  if dqext_use ge 0 then begin
     dq = mrdfits(infile,dqext_use,header_dq,silent=quiet)
     if n_elements(dq) eq 1 then begin
        if dq[0] eq -1 then begin
           nodq=1b
           header_dq=''
        endif
     endif
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
     dq = dat*0b
     ; in case any of dat has a nan
     inan = where(finite(dat,/nan),ctnan)
     if ctnan gt 0 then dq[inan] = 0b
     header_dq = header_dat
     sxaddpar,header_dq,'EXTNAME','DQ'
     if ~ keyword_set(quiet) then begin
        print,'IFSR_READCUBE: No dq extension.'
        print,'   Setting dq=0*dat and dq header to data header.'
     endif
  endelse
  if keyword_set(waveext) then begin
     wave = mrdfits(infile,waveext,header_wave,silent=quiet)
;    these were chosen to match MANGA defaults ...
     crval = 0
     crpix = 1
     cdelt = 1
  endif else begin
     header_wave = ''
  endelse

; Get #s of rows, columns, and wavelength pixels.
  datasize = size(dat)
  IF datasize[0] eq 3 then begin
     ncols = datasize[1]
     nrows = datasize[2]
     nz = datasize[3]
     wavedim = 3
     ndim = 3
     crvalstr = 'CRVAL3'
     crpixstr = 'CRPIX3'
     cdeltstr = 'CD3_3'
     cdeltstropt = 'CDELT3'
     cunitstr = 'CUNIT3'
  endif else if datasize[0] eq 2 then begin
     if ~ keyword_set(quiet) then $
        print,'IFSF_READCUBE: Reading 2D image. Assuming dispersion direction is'+$
           '   along first dimension.'
; Old logic: # wavelength pts > # cols.           
;     ncols = (datasize[1] gt datasize[2]) ? datasize[2] : datasize[1]
;     nz = (datasize[1] gt datasize[2]) ? datasize[1] : datasize[2]
;     wavedim = (datasize[1] gt datasize[2]) ? 1 : 2
; New logic: Assume dispersion along rows. Note that "ncols" then really
;     means "nrows"; keep this language for now.
     wavedim = 1
     nz = datasize[1]
     ncols = datasize[2]
     nrows = 1
     ndim = 2
     crvalstr = 'CRVAL1'
     crpixstr = 'CRPIX1'
     cdeltstr = 'CDELT1'
     cdeltstropt = 'CD1_1'
     cunitstr = 'CUNIT1'
  endif else if datasize[0] eq 1 then begin
     ncols = 1
     nrows = 1
     nz = datasize[1]
     wavedim = 1
     ndim = 1
     crvalstr = 'CRVAL1'
     crpixstr = 'CRPIX1'
     cdeltstr = 'CDELT1'
     cdeltstropt = 'CD1_1'
     cunitstr = 'CUNIT1'
  endif

  bunit = sxpar(header_dat,'BUNIT',silent=quiet,count=countbunit)
  if countbunit eq 0 then bunit = ''
  if ~novar then begin
     bunit_var = sxpar(header_var,'BUNIT',silent=quiet,count=countbunit_var)
     if countbunit_var eq 0 then bunit_var = ''
  endif else begin
     bunit_var = ''
  endelse

; Create wavelength array.
  if ~ keyword_set(waveext) then begin
     pix = dindgen(nz)+1d
     crval = double(sxpar(header_dat,crvalstr,silent=quiet,count=countval))
     crpix = double(sxpar(header_dat,crpixstr,silent=quiet,count=countpix))
     cdelt = double(sxpar(header_dat,cdeltstr,silent=quiet,count=countdel))
     cunit = sxpar(header_dat,cunitstr,silent=quiet,count=countunit)
     if countdel eq 0 then $
        cdelt = double(sxpar(header_dat,cdeltstropt,silent=quiet,count=countdel))
     if countunit eq 0 then cunit = ''
;    Assume absence of CRPIX quantity means it is 1
     if countpix eq 0 then crpix = 1d
     if countval eq 0 OR countdel eq 0 then begin
        print,'IFSF_READCUBE: WARNING -- missing a header parameter. Wavelength'
        print,'   solution will be incorrect.'
     endif
     wave = crval + cdelt*(pix-crpix)
     ; convert to A
     if strtrim(cunit) eq 'nm' then wave *= 10.

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
              dat[*,i] = interpol(datold[*,i],waveold,wave,/spline)
              var[*,i] = interpol(varold[*,i],waveold,wave,/spline)
              dq[*,i] = interpol(dqold[*,i],waveold,wave)
        endfor
     ENDIF
     IF datasize[0] eq 1 then begin
        dat = interpol(datold,waveold,wave,/spline)
        var = interpol(varold,waveold,wave,/spline)
        dq = interpol(dqold,waveold,wave)
     endif
     if ~ keyword_set(quiet) then $
        print,'IFSF_READCUBE: Interpolating DQ; values > 0.01 set to 1.'
     ibd = where(dq gt 0.01,ctbd)
     if ctbd gt 0 then dq[ibd] = 1
  endif
  
  if keyword_set(fluxnorm) then begin
     dat /= fluxnorm
     if ~novar then var /= fluxnorm*fluxnorm
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
         ndim: ndim,$
         wavedim: wavedim,$
         crval: crval,$
         cdelt: cdelt,$
         crpix: crpix,$
         cunit: cunit,$
         bunit: bunit,$
         bunit_var: bunit_var $
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
