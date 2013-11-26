; docformat = 'rst'
;
;+
;
; A lookup table for wavelengths.
;
; Wavelength sources:
;   1. NIST ASD
;   2. HeI 4686 -- found random wavelength compilation. In NIST, shows
;                  up as a complex of many fine structure lines with
;                  range of several tenths of an A.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Hash of wavelengths for specified input lines.
;
; :Params:
;    lines: in, required, type=strarr
;      Line labels for which to retrieve wavelengths.
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
;      2009jul08, DSNR, created
;      2010nov04, DSNR, wavelengths corrected from Mappings III to
;                       NIST values
;      2013sep, DSNR, added more lines for 1-slit GMOS config.
;      2013nov22, DSNR, documented, renamed, added copyright and
;                       license, and turned output into a hash
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
function ifsf_linelist,inlines

; Associated line labels:
  lines = hash()
  lines['Halpha'] = 6562.80d
  lines['Hbeta'] = 4861.32d
  lines['HeI6678'] = 6678.15d
  lines['HeI7065'] = 7065.19d
  lines['HeII4686'] = 4686.7d
  lines['[NI]5198'] = 5197.90d
  lines['[NI]5200'] = 5200.26d
  lines['[NII]6548'] = 6548.05d
  lines['[NII]6583'] = 6583.45d
  lines['[OI]5577'] = 5577.34d
  lines['[OI]6300'] = 6300.30d
  lines['[OI]6364'] = 6363.78d
  lines['[OIII]4959'] = 4958.91d
  lines['[OIII]5006'] = 5006.84d
  lines['[SII]6716'] = 6716.44d
  lines['[SII]6731'] = 6730.82d
  lines['NaD2'] = 5889.95d
  lines['NaD1'] = 5895.92d
  lines['[FeVII]5159'] = 5158.89d
  lines['[FeVII]5721'] = 5720.7d
  lines['[FeVII]6087'] = 6087.0d
  lines['[FeX]6375'] = 6374.51d
  
  outlines = hash()
  for i=0, inlines-1 do begin
     imatch = where(inlines[i] eq lines->keys(),ctmatch)
     if ctmatch eq 1 then outlines[inlines[i]] = lines[inlines[i]] $
     else print,'IFSF_LINELIST: ERROR: ',inlines[i],$
                ' not found in wavelength list.'
  endfor

  return,outlines

end
