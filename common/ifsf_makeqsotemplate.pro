; docformat = 'rst'
;
;+
;
;
; :Categories:
;    IFSFIT
;
; :Returns:
;
; :Params:
;
; :Keywords:
;    waveext: in, optional, type=integer
;      The extention number of a wavelength array.
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
;
; :Copyright:
;    Copyright (C) 2016--2018 David S. N. Rupke
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
pro ifsf_makeqsotemplate,infits,outxdr,dqext=dqext,waveext=waveext

   if ~ keyword_set(dqext) then dqext=3
   if ~ keyword_set(waveext) then waveext=0

   spec = ifsr_readspec(infits,ext=1,waveext=waveext)
   dq = ifsr_readspec(infits,ext=dqext)
   qsotemplate = {wave: double(spec[*,0]), flux: double(spec[*,1]), dq: dq[*,1]}
   save,qsotemplate,file=outxdr

end
