; docformat = 'rst'
;
;+
;
; Creates model emission line spectra for a given line, with taus summed 
; over all components.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Hash with 3 keys:
;      vel, type=dblarr, corresponding velocities w.r.t. systemic
;      tau, type=hash[lines,dblarr(ncols,nrows,nvel)]
;      tauerr, type=hash[lines,dblarr(ncols,nrows,nvel)]
;      cumtaunorm, type=hash[lines,dblarr(ncols,nrows,nvel)], cumulative taus 
;        at each wavelength, summed from low to high wavelengths, normalized so that
;        area under curve is 1.
;
; :Params:
;    abswav: in, required, type=hash
;    abssig: in, required, type=hash
;    abstau: in, required, type=hash
;    abstauerr: in, required, type=hash
;    ncomp: in, required, type=integer
;      Maximum number of components.
;    linelist: in, required, type=hash
;      Rest wavelengths of emission lines.
;    zref: in, required, type=double
;      Systemic redshift to which to refer output velocities.
;
; :Keywords:      
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
;      2016sep28, DSNR, created; copied from IFSF_CMPLINSPECMAPS
;      2016oct05, DSNR, updated error treatment to include wave, sigma errors
;    
; :Copyright:
;    Copyright (C) 2016 David S. N. Rupke
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
function ifsf_cmpcvdf_abs,abswav,abswaverr,abssig,abssigerr,abstau,abstauerr,$
                          ncomp,linewave,zref

   bad = 1d99
   c = 299792.458
   min = sqrt((machar(/double)).xmin)
   minexp = alog(min)

   modvelran = [-1d4,1d4]
   modvelstep = 1d
   nmod = (round(modvelran[1]-modvelran[0])/modvelstep)+1
   modvel = dindgen(nmod)
   modvel*=modvelstep
   modvel+=modvelran[0]
   beta = modvel/c
   dz = sqrt((1d + beta)/(1d - beta)) - 1d
   modwaves = linewave*(1d + dz)*(1d + zref)

   size_cube = size(abswav)
   
   abscvdf = hash()
   abscvdf['vel'] = modvel
   abscvdf['tau']=dblarr(size_cube[1],size_cube[2],nmod)
   abscvdf['tauerr']=dblarr(size_cube[1],size_cube[2],nmod)
   abscvdf['cumtaunorm']=dblarr(size_cube[1],size_cube[2],nmod)
   abscvdf['cumtaunormerr']=dblarr(size_cube[1],size_cube[2],nmod)

   for i=0,ncomp-1 do begin
      rbpktaus = rebin(abstau[*,*,i],size_cube[1],size_cube[2],nmod)
      rbpktauerrs = rebin(abstauerr[*,*,i],size_cube[1],size_cube[2],nmod)
      rbpkwaves = rebin(abswav[*,*,i],size_cube[1],size_cube[2],nmod)
      rbpkwaveerrs = rebin(abswaverr[*,*,i],size_cube[1],size_cube[2],nmod)
      rbsigmas = rebin(abssig[*,*,i],size_cube[1],size_cube[2],nmod)
      rbsigmaerrs = rebin(abssigerr[*,*,i],size_cube[1],size_cube[2],nmod)
      rbmodwaves = rebin(reform(modwaves,1,1,nmod),$
                         size_cube[1],size_cube[2],nmod)
      inz = where(rbsigmas gt 0 AND rbsigmas ne bad AND $
                  rbsigmaerrs ne bad AND $
; If sigma fits at lower limit then errors may be 0. Should not exclude these points.
;                  rbsigmaerrs gt 0 AND rbsigmaerrs ne bad AND $
                  rbpkwaves gt 0 AND rbpkwaves ne bad AND $
                  rbpkwaveerrs gt 0 AND rbpkwaveerrs ne bad AND $
                  rbpktaus gt 0 AND rbpktaus ne bad AND $
                  rbpktauerrs gt 0 AND rbpktauerrs ne bad,ctnz)
      if ctnz gt 0 then begin
         exparg = dblarr(size_cube[1],size_cube[2],nmod) -minexp
         exparg[inz] = ((rbmodwaves[inz]/rbpkwaves[inz] - 1d) / $
                        (rbsigmas[inz]/c))^2d / 2d
         i_no_under = where(exparg LT -minexp, ct_no_under)
         if ct_no_under gt 0 then begin
            abscvdf['tau',i_no_under] += $
               rbpktaus[i_no_under]*exp(-exparg[i_no_under])
            df_norm = rbpktauerrs[i_no_under]*exp(-exparg[i_no_under])
            df_wave = rbpktaus[i_no_under]*$
                      abs(rbmodwaves[i_no_under] - rbpkwaves[i_no_under])/$
                      (rbsigmas[i_no_under]/c*rbpkwaves[i_no_under])^2d*$
                      rbpkwaveerrs[i_no_under]*$
                      exp(-exparg[i_no_under])
            df_sig = rbpktaus[i_no_under]*$
                     (rbmodwaves[i_no_under] - rbpkwaves[i_no_under])^2d/$
                     (rbsigmas[i_no_under]/c*rbpkwaves[i_no_under])^2d * $
                     rbsigmaerrs[i_no_under]/rbsigmas[i_no_under]*$
                     exp(-exparg[i_no_under])
            dfsq = dblarr(size_cube[1],size_cube[2],nmod)
            dfsq = dfsq[i_no_under]
            i_no_under_2 = $
               where(df_norm gt min AND df_wave gt min AND df_sig gt min,$
                     ct_no_under_2)
            if ct_no_under_2 gt 0 then $
               dfsq[i_no_under_2] = df_norm[i_no_under_2]^2d + $
                                    df_wave[i_no_under_2]^2d + $
                                    df_sig[i_no_under_2]^2d
            abscvdf['tauerr',i_no_under] += dfsq
         endif
      endif
   endfor

   inz = where(abscvdf['tau'] gt 0,ctnz)
   if ctnz gt 0 then $
      abscvdf['tauerr',inz] = sqrt(abscvdf['tauerr',inz])
   dmodwaves = modwaves[1:nmod-1] - modwaves[0:nmod-2]
   dmodwaves = rebin(reform([dmodwaves[0],dmodwaves],1,1,nmod),$
                     size_cube[1],size_cube[2],nmod)
   taunorm = abscvdf['tau']*dmodwaves
   taunormerr = abscvdf['tauerr']*dmodwaves
   tauint = rebin(TOTAL(taunorm,3),size_cube[1],size_cube[2],nmod)
   inz = where(tauint ne 0,ctnz)
   if ctnz gt 0 then begin
      taunorm[inz]/=tauint[inz]
      taunormerr[inz]/=tauint[inz]
   endif

   abscvdf['cumtaunorm',*,*,0] = taunorm[*,*,0]
   for i=1,nmod-1 do $
      abscvdf['cumtaunorm',*,*,i] = $
         abscvdf['cumtaunorm',*,*,i-1] + taunorm[*,*,i]
   abscvdf['cumtaunormerr'] = taunormerr
 
   return,abscvdf

end
