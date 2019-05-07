; docformat = 'rst'
;
;+
;
; Pull previous stellar continuum fit from Binary FITS tables from Christy 
; Tremonti.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;
; :Params:
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
;      2019jan24, DSNR, created
;
; :Copyright:
;    Copyright (C) 2019 David S. N. Rupke
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
pro ifsf_makehizeacont,infile,cols,outxdr,fitran

   fxbopen,tablun,infile,1
   fxbreadm,tablun,cols,wave,cont

   cont /= max(cont)
   iran = where(wave ge fitran[0] AND wave le fitran[1])

   template = {lambda: wave[iran],$
               flux: [[cont[iran]]],$
               ages: [0d]}
   save,template,file=outxdr

end
