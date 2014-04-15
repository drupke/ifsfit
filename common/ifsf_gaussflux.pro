; docformat = 'rst'
;
;+
;
; Evaluate the integral of a Gaussian function over all x:
;   f(x) = norm * exp(-x/2sigma^2)
; Calls IFSF_GAUSSAREA to do the evaluation.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Structure with tags flux and flux_err.
;
; :Params:
;    norm: in, required, dblarr
;      Coefficient of exponential
;    sigma: in, required, dblarr
;      Standard deviation
;    
; :Keywords:
;    normerr: in, optional, dblarr
;      Error in norm.
;    sigerr: in, optional, dblarr
;      Error in sigma
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
;                       of sigma = 0, norm = 0 cases; moved to IFSFIT
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
function ifsf_gaussflux,norm,sigma,normerr=normerr,sigerr=sigerr

  ngauss = n_elements(norm)

  a = dblarr(ngauss)
  aerr = dblarr(ngauss)
  igda = where(sigma gt 0,ctgda)
  if ctgda gt 0 then begin
     a[igda] = 0.5d/sigma[igda]^2d
     if keyword_set(sigerr) then aerr[igda] = sigerr[igda]/sigma[igda]^3d
  endif
  gint = ifsf_gaussarea(a,aerr=aerr)
  flux = norm*gint.area

  if keyword_set(sigerr) then ginterr = gint.area_err $
  else ginterr=dblarr(ngauss)
  if ~ keyword_set(normerr) then normerr = dblarr(ngauss)

  fluxerr = dblarr(ngauss)
  igdn = where(norm gt 0 AND gint.area gt 0,ctgdn)
  if ctgdn gt 0 then $
     fluxerr[igdn] = flux[igdn]*sqrt((normerr[igdn]/norm[igdn])^2d + $
                                     (ginterr[igdn]/gint.area[igdn])^2d)
  
  outstr = {flux: flux, flux_err: fluxerr}

  return,outstr

end
