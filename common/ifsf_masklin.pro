; docformat = 'rst'
;
;+
;
; Mask data near emission lines.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Array of lambda-array indices indicating non-masked wavelengths.
;
; :Params:
;    lambda: in, required, type=dblarr
;      Wavelengths of spectrum
;    linelambda: in, required, type=dblarr(nlines)
;      Central wavelengths of lines to mask
;    halfwidth: in, required, type=dblarr(nlines)
;      Half width (in km/s) of masking region around each line
;
; :Keywords:
;    nomaskran: in, optional, type=dblarr(2\,nreg)
;      Set of lower and upper wavelength limits of regions *not* to
;      mask.
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
;      2009, DSNR, copied base code from Harus Jabran Zahid
;      2013oct, DSNR, documented, added nomaskran keyword
;      2013nov21, DSNR, renamed, added license and copyright
;      2013dec11, DSNR, renamed a couple of variables for clarity
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
function ifsf_masklin, lambda, linelambda, halfwidth, nomaskran=nomaskran

  c = 299792.458d

  indgd = indgen(n_elements(lambda))

  foreach line,linelambda do begin
    for i=0,n_elements(linelambda[line])-1 do begin
      ind_indgd = $
        where(lambda[indgd] lt linelambda[line,i]*$
                               (1d - halfwidth[line,i]/c) OR $
              lambda[indgd] gt linelambda[line,i]*$
                               (1d + halfwidth[line,i]/c),ct)
      if ct gt 0 then indgd = indgd[ind_indgd]
    endfor
  endforeach

  if keyword_set(nomaskran) then begin
     for j=0,n_elements(nomaskran)/2 - 1 do begin
        indgdadd = where(lambda ge nomaskran[0,j] AND $
                        lambda le nomaskran[1,j],ctadd)
        if ctadd gt 0 then begin
           if ct gt 0 then indgd = [indgd,indgdadd] else indgd = indgdadd
           ct += ctadd
        endif
     endfor
  endif

  indgd = indgd[sort(indgd)]
  indgd = indgd[uniq(indgd)]

  return, indgd

end
