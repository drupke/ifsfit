; docformat = 'rst'
;
;+
;
; Find stellar continuum for LZIFU data from LZIFU output file AND data cube.
; Normalization is measured from data and 
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
;      2017mar07, DSNR, created
;
; :Copyright:
;    Copyright (C) 2017 David S. N. Rupke
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
pro ifsf_makelzifucont,infile,incube,outxdr,fitran

   rcontnorm_cube = ifsf_readcube(infile,quiet=quiet,datext=2,varext=-1,dqext=-1)
   rcontnorm = rcontnorm_cube.dat
   rwave = rcontnorm_cube.wave
   rcontmask = mrdfits(infile,'R_CONT_MASK',/silent)
   rdat = mrdfits(incube,0,/silent)

;  Note that normalization assumed to be measured only from positive, finite,
;  unmasked points. Seems to work fine except for extremely low S/N spaxels.
   cubesize = size(rdat)
   for i=0,cubesize[1]-1 do begin
      for j=0,cubesize[2]-1 do begin
         igd = where(rcontmask[i,j,*] gt 0 AND $
                     finite(rcontnorm[i,j,*]) AND finite(rdat[i,j,*]) AND $
                     rcontnorm[i,j,*] gt 0 AND rdat[i,j,*] gt 0,ctgd)
         if ctgd gt 0 then begin
            norm = median(rcontnorm[i,j,igd]/rdat[i,j,igd],/double)
            rcontnorm[i,j,*] /= norm
         endif
      endfor
   endfor

   iran = where(rwave ge fitran[0] AND rwave le fitran[1])
      
   template = {lambda: rwave[iran],$
               flux: rcontnorm[*,*,iran]}
   save,template,file=outxdr

end
