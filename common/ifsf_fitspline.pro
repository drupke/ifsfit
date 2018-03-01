; docformat = 'rst'
;
;+
;
; Function to do b-spline interpolation to spectrum. Uses IDLUTILS 
; implementation (BSPLINE_ITERFIT).
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    The best fit continuum spectrum (over all wavelengths, not just
;    those fit).
;
; :Params:
;    galaxypint: in, required, type=string(A)
;    weight: in, required, type=dblarr(N)
;    template_flux: in, required, type=any
;      Ignored. When stellar templates are fit, contatins templates.
;    index: in, required, type=intarr
;      Contains indices of continuum regions to fit
;    ct_coeff: out, required, type=integer, default=0
;      When stellar templates are fit, contains best-fit coefficients.
;
; :Keywords:
;    argsbkpts: in, optional, type=structure
;      If set, contains parameters that control BSPLINE_BKPTS; i.e., how
;      the continuum is chopped up for splining.
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
;      2015jan19, DSNR, copied from IFSF_FITPOLY
;    
; :Copyright:
;    Copyright (C) 2015 David S. N. Rupke
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
function ifsf_fitspline,lambda,flux,weight,template_wave,template_flux,index,$
                        ct_coeff,quiet=quiet,refit=refit,argsbkpts=argsbkpts
                        

  err = 1d/sqrt(weight)
  mask = bytarr(n_elements(lambda))
  mask[index]=1b
  if keyword_set(argsbkpts) then $
     sset = bspline_iterfit(lambda,flux,invvar=1/weight,inmask=mask,$
                            _extra=argsbkpts) $
  else $
     sset = bspline_iterfit(lambda,flux,invvar=1/weight,inmask=mask,everyn=50)
  continuum = bspline_valu(lambda,sset)
  ct_coeff = 0
  
  if keyword_set(refit) then begin  
     icontinuum = continuum[index]
     ilambda=lambda[index]
     iflux=flux[index]
     iweight=weight[index]
     ierr = 1d/sqrt(iweight)

     for i = 0,n_elements(refit.ord)-1 do begin
        tmp_ind = where(lambda ge refit.ran[0,i] AND $
                        lambda le refit.ran[1,i])
        tmp_iind = where(ilambda ge refit.ran[0,i] AND $
                         ilambda le refit.ran[1,i])
        parinfo = replicate({value:0d},refit.ord[i]+1)
        tmp_pars = mpfitfun('poly',ilambda[tmp_iind],$
                            iflux[tmp_iind]-icontinuum[tmp_iind],$
                            ierr[tmp_iind],parinfo=parinfo,/quiet)
        continuum[tmp_ind] += poly(lambda[tmp_ind],tmp_pars)
     endfor
  endif

  return,continuum

end
