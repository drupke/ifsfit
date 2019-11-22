; docformat = 'rst'
;
;+
;
; Load Valdes et al. 2004 (ApJS, 152, 221) stellar spectra into a
; structure of wavelength and flux arrays.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    IDL save file. The file consists of a structure (named
;    'template') containing three tags: lambda, flux, and temps. Lambda
;    is an n-element array, while flux is and nxm array of the format
;    [n_wavelengths, n_ages]. Ages is the array of age values.
;
; :Params:
;    infile: in, required, type=string
;      Filename and path of input model table.
;    outfile: in, required, type=string
;      Filename and path of output file.
;
; :Keywords:
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
;      2018apr14, DSNR, copied from IFSF_GDTEMP
;    
; :Copyright:
;    Copyright (C) 2018--2019 David S. N. Rupke
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
pro ifsf_valdestemp,infile,outfile,wavelo=wavelo,wavehi=wavehi

   bad = 1d99
   valdesbad = -2147483647

   disp = 0.4d
   if ~ keyword_set(wavelo) then wavelo = 3465.0d
   if ~ keyword_set(wavehi) then wavehi = 7000.0d
   
;  get filenames
   readcol,infile,files,/silent,format='(A)'
  
;  process each file
   nuse = 0
   foreach file,files do begin
      openr,lun,file,/get_lun
      usespec = 1b
      isheadline = 1b
      while isheadline AND usespec do begin
         headline = strarr(1)
         readf,lun,headline
         if strmatch(headline[0],'#K*') then begin
            if strmatch(headline[0],'*OBJECT*') then begin
               split = strsplit(headline[0]," '",/extract)
               thisobj = split[2]
            endif
            if strmatch(headline[0],'*TEFF*') then begin
               split = strsplit(headline[0],/extract)
               thisteff = double(split[2])
            endif
            if strmatch(headline[0],'*LOGG*') then begin
               split = strsplit(headline[0],/extract)
               thislogg = double(split[2])
            endif
            if strmatch(headline[0],'*FEH*') then begin
               split = strsplit(headline[0],/extract)
               thisfeh = double(split[2])
            endif
            if strmatch(headline[0],'*COVERAGE*') then begin
               split = strsplit(headline[0],' -',/extract)
               wave0 = double(split[2])
               wave1 = double(split[3])
               if wave0 gt wavelo OR wave1 lt wavehi then usespec=0b
            endif
            if strmatch(headline[0],'*GAPS*') then begin
               split = strsplit(headline[0],' -,',/extract)
               ngaps = (n_elements(split)-1)/2
               for i=0,ngaps-1 do begin
                  wgap0 = double(split[2+2*i])
                  wgap1 = double(split[3+2*i])
                  dgap = wgap1-wgap0
                  if dgap gt 50d then usespec=0b
               endfor
            endif
         endif else begin
            isheadline = 0b
         endelse
      endwhile

      if usespec then begin
         nuse++
         if nuse eq 1 then begin
            readcol,file,thiswave,thisflux,/silent,comment='#',format='(D,D)'
            ilo = value_locate(thiswave,wavelo)
            ihi = value_locate(thiswave,wavehi)
            wave = thiswave[ilo:ihi]
            thisflux /= max(thisflux[ilo:ihi])
            flux = thisflux[ilo:ihi]
            obj = thisobj
            if thisteff ne valdesbad then teff = thisteff else teff = bad
            if thislogg ne valdesbad then logg = thislogg else logg = bad
            if thisfeh ne valdesbad then feh = thisfeh else feh = bad
         endif else begin
            readcol,file,thiswave,thisflux,/silent,comment='#',format='(D,D)'
            ilo = value_locate(thiswave,wavelo)
            ihi = value_locate(thiswave,wavehi)
            thisflux /= max(thisflux[ilo:ihi])
            flux = [[flux],[thisflux[ilo:ihi]]]
            obj = [obj,thisobj]
            if thisteff ne valdesbad then teff = [teff,thisteff] $
               else teff = [teff,bad]
            if thislogg ne valdesbad then logg = [logg,thislogg] $
               else logg = [logg,bad]
            if thisfeh ne valdesbad then feh = [feh,thisfeh] $
               else feh = [feh,bad]         
         endelse
      endif

      free_lun,lun
   endforeach
  
   template = {lambda: wave, flux: flux, obj: obj, $
               teff: teff, logg: logg, feh: feh}

   save,template,filename=outfile
  
end
