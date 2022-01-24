; docformat = 'rst'
;
;+
;
; Read single spectrum and output to structure QSOTEMPLATE in an IDL save file.
;
; :Categories:
;    IFSFIT
;
; :Returns:
; 
;    structure QSOTEMPLATE in an IDL save file
;
; :Params:
;    infits: in, required, type=string
;    outxdr: in, required, type=string
;
; :Keywords:
;    waveext: in, optional, type=integer
;      The extension number of a wavelength array.
;    datext: in, optional, type=integer, default=1
;      The extension number of the data array. Set to 0 or -1 for extension 0.
;    varext: in, optional, type=integer, default=2
;      The extension number of the variance array.
;    dqext: in, optional, type=integer, default=3
;      The extension number of the dq array.
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
;      2016dec14, DSNR, created
;      2018feb08, DSNR, added WAVEEXT keyword
;      2018may24, DSNR, added DATEXT keyword
;      2022jan05, DSNR, added VAREXT keyword; cleaned up documenation
;
; :Copyright:
;    Copyright (C) 2016--2021 David S. N. Rupke
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
pro ifsf_makeqsotemplate,infits,outxdr,datext=datext,dqext=dqext,varext=varext,$
   waveext=waveext

   if ~ keyword_set(waveext) then waveext=0
   if ~ keyword_set(datext) then datext=1
   if datext eq -1 then datext=0
   if ~ keyword_set(varext) then varext=2
   if ~ keyword_set(dqext) then dqext=3

   spec = ifsr_readspec(infits,ext=datext,waveext=waveext)
   var = ifsr_readspec(infits,ext=varext)
   dq = ifsr_readspec(infits,ext=dqext)
   qsotemplate = {wave: double(spec[*,0]), flux: double(spec[*,1]), $
      var: double(var[*,1]), dq: dq[*,1]}
   save,qsotemplate,file=outxdr

end
