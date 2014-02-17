; docformat = 'rst'
;
;+
;
; Add smooth functions to stellar continuum templates for spectral
; fitting: polynomials of form x^i (i.e., 1, x, x^2, etc.), their
; array reversals in wavelength space (i.e., reverse(x), etc.), and
; optionally 8 exponentials (of form +/-e^-x, +/-e^-2x, +/-e^-4x,
; +/-e^-8x).
;
; :Categories:
;    IFSFIT
;
; :Returns:
;
;    Array of continuum fitting templates, of type=dblarr(N,M'), with
;    N pixels and M' templates. M' is the number of original templates
;    plus NTERMS polynomials and 8 exponentials (i.e., M' = M + NTERMS
;    + 8).
;
; :Params:
;    template: in, required, type=dblarr(N\,M)
;      Stellar continuum templates for continuum fitting, with N
;      pixels and M templates.
;
; :Keywords:
;    nterms: in, required, type=integer, default=3
;      Use to prevent detailed output to screen. Default is to print
;      detailed output.
;    addexp: in, required, type=byte, default=0
;      Set to add 8 exponential terms.
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
;      2010apr13, DSNR, added exponentials
;      2013oct, DSNR, renamed, added documentation, added
;                     reverse-ordered x^i terms
;      2013nov11, DSNR, renamed, added license and copyright 
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
function ifsf_addpoly2temp, template, nterms=nterms, addexp=addexp

  if ~ keyword_set(nterms) then nterms = 3

  nwave = double(n_elements(template[*,0]))
  new_temp = template

  new_temp = [ [new_temp], [dblarr(nwave)+1d] ]
  for i=1, nterms-1 do begin
     new_temp = [ [new_temp], [(dindgen(nwave)/nwave)^i] ]
     new_temp = [ [new_temp], [reverse((dindgen(nwave)/nwave)^i)] ]
  endfor

  if keyword_set(addexp) then begin
     new_temp = [ [new_temp], [exp(-dindgen(nwave)/nwave)] ]
     new_temp = [ [new_temp], [exp(-2d*dindgen(nwave)/nwave)] ]
     new_temp = [ [new_temp], [exp(-4d*dindgen(nwave)/nwave)] ]
     new_temp = [ [new_temp], [exp(-8d*dindgen(nwave)/nwave)] ]
     new_temp = [ [new_temp], [1d -exp(-dindgen(nwave)/nwave)] ]
     new_temp = [ [new_temp], [1d -exp(-2d*dindgen(nwave)/nwave)] ]
     new_temp = [ [new_temp], [1d -exp(-4d*dindgen(nwave)/nwave)] ]
     new_temp = [ [new_temp], [1d -exp(-8d*dindgen(nwave)/nwave)] ]
  endif

  return, new_temp

end
