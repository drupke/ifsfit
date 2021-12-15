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
;      2018may11, DSNR, added option for non-GD05 templates
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
pro ifsf_makestellartemplate,inxdr,outxdr,refit=refit,$
                             stellartemplates=stellartemplates

   if keyword_set(stellartemplates) then $
      startempfile = stellartemplates $
   else $
     startempfile = '/Users/drupke/Documents/stellar_models/'+$
                    'gonzalezdelgado/SSPGeneva_z020.sav'
   restore,startempfile
   restore,inxdr
   if keyword_set(refit) then stelspec = template.flux # struct.ct_coeff.stel $
   else stelspec = template.flux # struct.ct_coeff
   stelspec /= median(stelspec)
   template = {lambda: template.lambda,$
               flux: [[stelspec]],$
               ages: [0d]}
   save,template,file=outxdr

end
