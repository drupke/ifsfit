; docformat = 'rst'
;
;+
;
; Evaluate the integral of a Gaussian function over all x:
;   f(x) = norm * exp(-ax^2)
;   Area = sqrt(Pi/a)
; If a = 0, area is set to 0.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Structure with tags area and, optionally, area_err.
;
; :Params:
;    a: in, required, dblarr
;      Coefficient of exponential argument
;
; :Keywords:
;    aerr: in, optional, dblarr
;      Error in a
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
;      2014apr10, DSNR, documented, added license and copyright; added treatment
;                       of a = 0 case; moved to IFSFIT
;
; :Copyright:
;    Copyright (C) 2014 David S. N. Rupke
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
function ifsf_gaussarea,a,aerr=aerr

  ngauss = n_elements(a)
  sqrtpi = sqrt(!DPI)

  out = dblarr(ngauss)
  igda = where(a gt 0,ctgda)
  if ctgda gt 0 then begin
     sqrta = sqrt(a[igda])
     out[igda] = sqrtpi/sqrta
  endif
  outstr={area:out}
  
  if keyword_set(aerr) then begin
     outerr = dblarr(ngauss)
     if ctgda gt 0 then outerr[igda] = out[igda]*0.5d/a[igda]*aerr[igda]
     outstr = add_tag(outstr,outerr,'area_err')
  endif

  return,outstr

end
