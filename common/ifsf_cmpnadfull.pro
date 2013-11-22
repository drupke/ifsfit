; docformat = 'rst'
;
;+
;
; Take the data used for the Na D fit and the best-fit parameters and
; compute the best-fit spectrum at the same wavelengths. Optionally
; compute equivalent width from fitted profile (including absorption
; only).
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    type=dblarr(Nwaves,5); five elements of second dimension are
;    wavelengths, data, errors, model (abs+em), and model (abs only)
;
; :Params:
;    datfile: in, required, type=string
;      Full path and name of file containing data (with emission lines
;      subtracted), wavelength, and error in three columns.
;    parfile: in, required, type=
;      Parameter file output by Na D fitting routine.
;    z: in, required, type=double
;      Redshift.
;
; :Keywords:
;    weq: in, optional, type=double
;      Return equivalent width computed from fitted profile. NOTE:
;      INCLUDES ABSORPTION ONLY.
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
;      2011jul28, DSNR, created
;      2013nov21, DSNR, documented, renamed, added license and copyright 
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
function ifsf_cmpnadfull,datfile,parfile,z,weq=weq

  readcol300,datfile,specwave,specflux,specerr,/silent,/skip,$
             format='(D,D,D)'

  nadpars = ifsf_readnadpar(parfile)

  modflux = dblarr(n_elements(specflux))+1d
  modabs = dblarr(n_elements(specflux))+1d
  for i=0,nadpars.nabs-1 do modabs *= ifsf_cmpnad(specwave,nadpars.abs[*,i])
  for i=0,nadpars.nem-1 do begin
     arg = ((specwave-nadpars.em[1,i])/$
            (nadpars.em[1,i]*nadpars.em[2,i]/299792d))^2d
     modflux += nadpars.em[0,i]*exp(-arg)
  endfor
  modflux *= modabs

; Compute weq
  dlam = specwave - shift(specwave,1)
  weq = total((1-modabs)*dlam)

  return,[[specwave],[specflux],[specerr],[modflux],[modabs]]

end
