; docformat = 'rst'
;
;+
;
; Function to fit polynomial continuum to spectrum.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    The best fit continuum spectrum (over all wavelengths, not just
;    those fit).
;
; :Params:
;    lambda: in, required, type=dblarr(N)
;    flux: in, required, type=dblarr(N)
;    weight: in, required, type=dblarr(N)
;    template_flux: in, required, type=any
;      Ignored. When stellar templates are fit, contatins templates.
;    index: in, required, type=intarr
;      Contains indices of continuum regions to fit
;    ct_coeff: out, required, type=integer, default=0
;      When stellar templates are fit, contains best-fit coefficients.
;
; :Keywords:
;    fitord: in, optional, type=integer
;      Specifies order of polynomial renormalization (default = 3).
;    quiet: in, optional, type=byte
;    refit: in, optional, type=structure
;      If set, contains structure with array of continuum regions to
;      re-fit [tag RAN, type=dblarr(2, X), where X is the number of
;      regions to refit] and array of polynomial orders to fit [tag
;      ORD, type=intarr(X)].
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
;    Change History::
;      2014jan29, DSNR, copied from IFSF_FITCONT
;      2015jan18, DSNR, Switched from POLY_FIT to MPFIT for fit; found that
;                       order much above 6 gave bizarre results otherwise.
;                       There are probably still limits to how high fit order 
;                       can go!
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
function ifsf_s7cont,lambda,flux,weight,template_lambdaz,template_flux,index,$
                     ct_coeff,zstar,quiet=quiet,colrow=colrow
                     
  continuum = template_flux[colrow[0]-1,colrow[1]-1,*]
  ct_coeff=0d
  return,continuum

end
