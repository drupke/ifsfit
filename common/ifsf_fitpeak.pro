; docformat = 'rst'
;
;+
;
; Function to fit a peak in a spectrum. Developed to serve as a "continuum" 
; (for, e.g., quasar spectra), perhaps as an input to IFSF_FITMULTICONT.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    The best fit  spectrum (over all wavelengths, not just
;    those fit).
;
; :Params:
;    galaxypint: in, required, type=string(A)
;    weight: in, required, type=dblarr(N)
;    template_wave: in, required, type=any
;      Ignored. When stellar templates are fit, contatins template wavelengths.
;    template_flux: in, required, type=any
;      Ignored. When stellar templates are fit, contatins templates.
;    index: in, required, type=intarr
;      Contains indices of continuum regions to fit
;    ct_coeff: out, required, type=integer, default=0
;      When stellar templates are fit, contains best-fit coefficients.
;
; :Keywords:
;    baseline: in, optional, type=integer, default=0
;      Set to 1 to fit a constant baseline, 2 for a linear polynomial.
;    estimates: in, optional, type=dblarr(Npar)
;      Best-fit estimates.
;    fixed: in, optional, type=bytarr(Npar)
;      Set to 1 for each parameter that should be fixed.
;    limited:
;    limits:
;    peaktype: in, required, type=string
;      Set to 'gaussian', 'lorentzian', or 'moffat' (Moffat 
;      not properly set up yet ...)
;    quiet: in, optional, type=byte
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
;      2018feb15, DSNR, copied from IFSF_FITPOLY
;    
; :Copyright:
;    Copyright (C) 2018 David S. N. Rupke
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
function ifsf_fitpeak,lambda,flux,weight,template_wave,template_flux,index,$
                      ct_coeff,quiet=quiet,refit=refit,$
                      estimates=estimates,baseline=baseline,$
                      fixed=fixed,peaktype=peaktype,limits=limits,limited=limited
                        
  ilambda=lambda[index]
  iflux=flux[index]
  iweight=weight[index]
  ierr = 1d/sqrt(iweight)

  nterms = 3
  if keyword_set(baseline) then begin
     if baseline eq 1 then nterms++
     if baseline eq 2 then nterms += 2
  endif
  parinfo = replicate({value:0d,fixed:0,limits:[0d,0d],limited:[0d,0d]},nterms)
  if keyword_set(fixed) then parinfo.fixed = fixed
  if keyword_set(limited) then parinfo.limited = limited
  if keyword_set(limits) then parinfo.limits = limits

  if peaktype eq 'gaussian' then gaussian=1 else gaussian=0
  if peaktype eq 'lorentzian' then lorentzian=1 else lorentzian=0
  if peaktype eq 'moffat' then moffat=1 else moffat=0

  yfit = mpfitpeak(ilambda,iflux,pars,weights=iweight,parinfo=parinfo,$
                   estimates=estimates,nterms=nterms,$
                   status=status,gaussian=gaussian,lorentzian=lorentzian,$
                   moffat=moffat,quiet=quiet)
  if peaktype eq 'lorentzian' then yfit = drt_lorentz(lambda,pars)
  if peaktype eq 'gaussian' then yfit = drt_gaussian(lambda,pars)
;  if peaktype eq 'moffat' then yfit = drt_moffat(lambda,pars)
  ct_coeff = 0

  return,yfit

end
