; docformat = 'rst'
;
;+
;
; Refine initial guess for redshift by finding emission-line peak.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Redshift of emission-line peak.
;
; :Params:
;    zinit: in, required, type=double
;      Initial redshift guess.
;    lambda: in, required, type=dblarr(N)
;      Wavelengths of spectrum.
;    flux: in, required, type=dblarr(N)
;      Fluxes of spectrum.
;
; :Keywords:
;    lref: in, optional, type=double, default=6562.8
;      Rest wavelength of expected emission line
;    searchwidth: in, optional, type=double, default=10
;      Half-width of search region, in A.
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
;      2010jan19, DSNR, created
;      2013nov25, DSNR, documented, renamed, added license and copyright 
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
function ifsf_refinez,zinit,lambda,flux,lref=lref,searchwidth=searchwidth

  if ~ keyword_set(searchwidth) then searchwidth=10d
  if ~ keyword_set(lref) then lref = 6562.8d
  
  searchwindow = [lref-searchwidth,lref+searchwidth]
  iref0 = value_locate(lambda,searchwindow[0])
  iref1 = value_locate(lambda,searchwindow[1])
  lambdawindow = lambda[iref0:iref1]
  fluxwindow = flux[iref0:iref1]
  fmax = max(fluxwindow,imax)
  znew = zinit + (lambdawindow[imax] / lref) - 1d
  
  return,znew
  
end
