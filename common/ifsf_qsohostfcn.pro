; docformat = 'rst'
;
;+
;
; Function to fit a continuum with a QSO template plus a smooth additive 
; template. Presently, the additive template is the sum of polynomials up to 
; order FITORD, which can also be multiplied by exponential functions. The 
; multiplicative exponentials are as follows:
; 
; EXPTERMS = 0: none
; EXPTERMS = 1: C_1 exp^(-(i-i_1)/N),
;     where i is the wavelength index
;           N is the number of wavelength points
;           and C_1 and i_1 are fitted parameters
; EXPTERMS = 2: first exp term + C_2 (1 - exp^(-(ilo-i_2)/N))
; EXPTERMS = 3: first two exp terms + C_3 (1 - exp^(-(j-i_3)/N))
;     where j is the reverse wavelength index (w.r.t. the last wavelength)
; EXPTERMS = 4: first two exp terms + C_4 (1 - exp^(-(j-i_4)/N))
; 
; The QSO is also multiplied by a sum of polynomials up to order QSOORD, plus
; whatever exponential terms are specified.
;
; Order of parameter array:
;
; Notes:
; - Fitting does not seem sensitive to how we define xc.  x - max(x)
;   gives identical results.
; - Fitting seems most sensitive to whether or not we include the
;   first exponential term for both the stellar and QSO continua.
;   It's not clear that the second exp. term matters.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    The best fit continuum spectrum (over all wavelengths, not just
;    those fit).
;
; :Params:
;    x: in, required, type=dblarr(N)
;      Wavelength array
;    p: in, required, type=dblarr(N)
;      Parameter array
;    ymod: out, required, type=dblarr(N)
;      Model flux
;
; :Keywords:
;    blrterms: in, optional, type=integer, default=0
;      Number of Gaussian terms to include (3 x number of Gaussian components).
;    expterms: in, optional, type=integer, default=0
;      Number of exponential terms by which to normalize, up to 4
;    fitord: in, optional, type=integer, default=3
;      Specifies order of additive renormalization
;    quiet: in, optional, type=byte
;    qsoord: in, optional, type=integer, default=3
;      Specifies order of multiplicative renormalization
;    qsoonly: in, optional, type=byte
;      Do not add stellar continuum.
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
;      2010jul17, DSNR, created
;      2015jan24, DSNR, copied from GMOS_QSO_CONT_FCN
;      2015sep04, DSNR, added option to include broad components as a way
;                       to either fit a BLR or to deal with scattered
;                       light issues with BLR emission (PG1411+442 et al.)
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
;-;
pro ifsf_qsohostfcn,x,p,ymod,fitord=fitord,$
                    qsoflux=qsoflux,qsoord=qsoord,$
                    qsoonly=qsoonly,expterms=expterms,$
                    blrterms=blrterms

  if ~ keyword_set(fitord) then fitord=3
  if ~ keyword_set(qsoord) then qsoord=3
  if ~ keyword_set(expterms) then expterms=0
  if ~ keyword_set(blrterms) then blrterms=0

  npar = n_elements(p)
  nx = n_elements(x)
  dx = max(x) - min(x)
  xc = x - min(x)
  xc_hi = max(x) - x

  ymod = dblarr(nx)
  qsoscl = dblarr(nx)
  itot = fitord+qsoord+2*expterms

; Additive polynomial
  if ~ keyword_set(qsoonly) then begin
     for i=0,fitord-1 do $
        ymod += p[i]        * (xc-p[itot+i])^double(i)
     if expterms ge 1 then $
        ymod += p[fitord]   * exp(-(xc-p[itot+fitord])/dx)
     if expterms ge 2 then $
        ymod += p[fitord+1] * (1d - exp(-(xc-p[itot+fitord+1])/dx))
     if expterms ge 3 then $
        ymod += p[fitord+2] * exp(-(xc_hi-p[itot+fitord+2])/dx)
     if expterms eq 4 then $
        ymod += p[fitord+3] * (1d - exp(-(xc_hi-p[itot+fitord+3])/dx))
  endif

; QSO continuum
  for i=0,qsoord-1 do $
     qsoscl += p[fitord+expterms+i] * $
                (xc-p[itot+fitord+expterms+i])^double(i)
  if expterms ge 1 then $
     qsoscl += p[fitord+expterms+qsoord] * $
               exp(-(xc-p[itot+fitord+expterms+qsoord])/dx)
  if expterms ge 2 then $
     qsoscl += p[fitord+expterms+qsoord+1] * $
               (1d - exp(-(xc-p[itot+fitord+expterms+qsoord+1])/dx)) 
  if expterms ge 3 then $
     qsoscl += p[fitord+expterms+qsoord+2] * $
               exp(-(xc_hi-p[itot+fitord+expterms+qsoord+2])/dx)
  if expterms eq 4 then $
     qsoscl += p[fitord+expterms+qsoord+3] * $
               (1d - exp(-(xc_hi-p[itot+fitord+expterms+qsoord+3])/dx))
  ymod += qsoscl*qsoflux

; BLR model
  ngauss = fix(blrterms/3)
  for i=0,ngauss-1 do $
     ymod += double(gaussian(x,p[npar-blrterms+0+3*i:npar-blrterms+3*(i+1)-1]))

end
