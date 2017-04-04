; docformat = 'rst'
;
;+
;
; Function to fit a continuum with a QSO template plus a smooth additive 
; template. The additive template is the sum of Legendre polynomials up to 
; order FITORD. The QSO is also multiplied by a sum of Legendre polynomials 
; up to order QSOORD. We use Legendre polynomials rather than a "standard"
; polynomial b/c it is more well-behaved with MPCURVEFIT as the order increases.
; In principle the two should yield the same result, but in practice Legendre
; polynomials are more faithful to the data for high orders (for unknown
; numerical reasons, though probably having something to do with Legendre
; polynomials being orthogonal over -1 < x < 1).
;
; Notes:
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
;    hostonly: in, optional, type=byte
;      Output host only.
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
;      2016sep02, DSNR, added HOSTONLY keyword; set so that BLR emission
;                       goes in QSO only spectrum
;      2016sep13, DSNR, fixed small bug in last if statement
;      2016sep21, DSNR, added abs() in front of each ymod term to ensure
;                       host spectrum cannot go negative!
;      2016nov11, DSNR, re-wrote to use only Legendre polynomials
;    
; :Copyright:
;    Copyright (C) 2015--2016 David S. N. Rupke
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
pro ifsf_qsohostfcn,x,p,ymod,nxfull=nxfull,ixfull=ixfull,$
                    hostonly=hostonly,qsoonly=qsoonly,blronly=blronly,$
                    qsoflux=qsoflux,qsoscl=qsoscl,blrterms=blrterms

  if ~ keyword_set(blrterms) then blrterms=0

  npar = n_elements(p)
  
  nx = n_elements(x)
  if ~ keyword_set(nxfull) then nxfull=nx
  if ~ keyword_set(ixfull) then ixfull=indgen(nx)
  
  xfull_norm01 = dindgen(nxfull)/double(nxfull-1) ; array of numbers from 0 to 1
  xfull_norm10 = reverse(xfull_norm01)

  x_norm01 = xfull_norm01[ixfull]
  x_norm10 = xfull_norm10[ixfull]
  
  ymod = dblarr(nx)

; Additive polynomial.
  if ~ keyword_set(qsoonly) AND ~ keyword_set(blronly) then begin
     ymod += p[0]*exp(-p[1]*x_norm01)
     ymod += p[2]*exp(-p[3]*x_norm10)
     ymod += p[4]*(1d - exp(-p[5]*x_norm01))
     ymod += p[6]*(1d - exp(-p[7]*x_norm10))
  endif

; QSO continuum
  if ~ keyword_set(hostonly) AND ~ keyword_set(blronly) then begin
     qsoscl = dblarr(nx)
     qsoscl += p[8]*exp(-p[9]*x_norm01)
     qsoscl += p[10]*exp(-p[11]*x_norm10)
     qsoscl += p[12]*(1d - exp(-p[13]*x_norm01))
     qsoscl += p[14]*(1d - exp(-p[15]*x_norm10))
     ymod += qsoscl*qsoflux
  endif

; BLR model
  if ~ keyword_set(hostonly) AND keyword_set(blrterms) then begin
     ngauss = fix(blrterms/3)
     for i=0,ngauss-1 do $
        ymod += double(gaussian(x,p[npar-blrterms+0+3*i:npar-blrterms+3*(i+1)-1]))
  endif

end
