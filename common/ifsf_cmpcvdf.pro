; docformat = 'rst'
;
;+
;
; Creates model emission line spectra for a given line, with fluxes summed 
; over all components.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Hash with 3 keys:
;      vel, type=dblarr, corresponding velocities w.r.t. systemic
;      flux, type=hash[lines,dblarr(ncols,nrows,nvel)]
;      fluxerr, type=hash[lines,dblarr(ncols,nrows,nvel)]
;      cumfluxnorm, type=hash[lines,dblarr(ncols,nrows,nvel)], cumulative fluxes 
;        at each wavelength, summed from low to high wavelengths, normalized so that
;        area under curve is 1.
;
; :Params:
;    emlwav: in, required, type=hash
;    emlsig: in, required, type=hash
;    emlflx: in, required, type=hash
;    emlflxerr: in, required, type=hash
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
function ifsf_cmpcvdf,emlwav,emlwaverr,emlsig,emlsigerr,emlflx,emlflxerr,$
                      ncomp,linelist,zref

   bad = 1d99
   c = 299792.458
;  Code from GAUSSIAN function in ASTROLIB: Get smallest value expressible
;  on computer. SQRT because sometimes we multiply two v. small values by each
;  other.
   min = sqrt((machar(/double)).xmin)
   minexp = alog(min)

   modvelran = [-1d4,1d4]
   modvelstep = 10d
   nmod = (round(modvelran[1]-modvelran[0])/modvelstep)+1
   modvel = dindgen(nmod)
   modvel*=modvelstep
   modvel+=modvelran[0]
   
   emlcvdf = hash()
   emlcvdf['vel'] = modvel
   emlcvdf['flux'] = hash()
   emlcvdf['fluxerr'] = hash()
   emlcvdf['cumfluxnorm'] = hash()
   emlcvdf['cumfluxnormerr'] = hash()

   outlines = emlflx['ftot']->keys()
   size_cube = size(emlflx['ftot',outlines[0]])
      
   foreach line,outlines do begin
      emlcvdf['flux',line]=dblarr(size_cube[1],size_cube[2],nmod)
      emlcvdf['fluxerr',line]=dblarr(size_cube[1],size_cube[2],nmod)
      emlcvdf['cumfluxnorm',line]=dblarr(size_cube[1],size_cube[2],nmod)
      emlcvdf['cumfluxnormerr',line]=dblarr(size_cube[1],size_cube[2],nmod)
      beta = modvel/c
;     redshift w.r.t. galaxy systemic
      dz = sqrt((1d + beta)/(1d - beta)) - 1d
      modwaves = linelist[line]*(1d + dz)*(1d + zref)
;     This comes from assuming
;     1 + dz = l_gas^gal/l_gal^gal
;     where _gas is the gas
;           _gal is the galaxy (systemic)
;           ^Earth is w.r.t. Earth
;           ^gal is w.r.t. the galaxy (systemic)
;     An observer in the galaxy will measure
;        l_gal^gal = rest wavelength of transition
;        l_gas^gal = red/blue-shifted wavelength due to gas motion
;        dz = corresponding redshift
;     An observer on Earth will measure
;        l_gal^Earth = (1+z_gal)l_gal^gal
;        l_gas^Earth = (1+z_gal)l_gas^gal
;     if we assume the gas parcel assumes the same cosmological redshift as the
;     galaxy.
;     Then 
;        1+dz = l_gas^Earth / l_gal^Earth
;             = observed wavelength / [rest wavelength*(1+z_gal)]
      for i=0,ncomp-1 do begin
         cstr='c'+string(i+1,format='(I0)')
         rbpkfluxes = $
            rebin(emlflx['f'+cstr+'pk',line],size_cube[1],size_cube[2],nmod)
         rbpkfluxerrs = $
            rebin(emlflxerr['f'+cstr+'pk',line],size_cube[1],size_cube[2],nmod)
         rbpkwaves = $
            rebin(emlwav[cstr,line],size_cube[1],size_cube[2],nmod)
         rbpkwaveerrs = $
            rebin(emlwaverr[cstr,line],size_cube[1],size_cube[2],nmod)
         rbsigmas = $
            rebin(emlsig[cstr,line],size_cube[1],size_cube[2],nmod)
         rbsigmaerrs = $
            rebin(emlsigerr[cstr,line],size_cube[1],size_cube[2],nmod)
         rbmodwaves = rebin(reform(modwaves,1,1,nmod),$
                            size_cube[1],size_cube[2],nmod)
;        Not sure how much of this is necessary but better safe than sorry.
;        Don't include sigmaerrs b/c error is 0 when sigma pegs at a limit.
         inz = where(rbsigmas gt 0 AND rbsigmas ne bad AND $
;                     rbsigmaerrs gt 0 AND rbsigmaerrs ne bad AND $
                     rbpkwaves gt 0 AND rbpkwaves ne bad AND $
                     rbpkwaveerrs gt 0 AND rbpkwaveerrs ne bad AND $
                     rbpkfluxes gt 0 AND rbpkfluxes ne bad AND $
                     rbpkfluxerrs gt 0 AND rbpkfluxerrs ne bad,ctnz)
         if ctnz gt 0 then begin
;           The minexp at the end ensures that the logical WHERE below doesn't
;           add back in points already rejected in the logical WHERE above -- 
;           i.e. I_NO_UNDER should be a subset of INZ.
            exparg = dblarr(size_cube[1],size_cube[2],nmod) -minexp
            exparg[inz] = ((rbmodwaves[inz]/rbpkwaves[inz] - 1d) / $
                           (rbsigmas[inz]/c))^2d / 2d
;           Set low values to 0 to avoid floating underflow.
            i_no_under = where(exparg LT -minexp, ct_no_under)
;            if ct_no_under gt 0 then $
;               iuse = cgsetintersection(inz,i_no_under,count=ctuse) $
;            else ctuse = 0
            if ct_no_under gt 0 then begin
               emlcvdf['flux',line,i_no_under] += $
                  rbpkfluxes[i_no_under]*exp(-exparg[i_no_under])
;              Flux errors. Calculated using Gaussian formulas and assuming
;              independent errors on each quantity. Usual error propagation rules.
               df_norm = rbpkfluxerrs[i_no_under]*exp(-exparg[i_no_under])
               df_wave = rbpkfluxes[i_no_under]*$
                         abs(rbmodwaves[i_no_under] - rbpkwaves[i_no_under])/$
                         (rbsigmas[i_no_under]/c*rbpkwaves[i_no_under])^2d*$
                         rbpkwaveerrs[i_no_under]*$
                         exp(-exparg[i_no_under])
               df_sig = rbpkfluxes[i_no_under]*$
                        (rbmodwaves[i_no_under] - rbpkwaves[i_no_under])^2d/$
                        (rbsigmas[i_no_under]/c*rbpkwaves[i_no_under])^2d * $
                        rbsigmaerrs[i_no_under]/rbsigmas[i_no_under]*$
                        exp(-exparg[i_no_under])
;              Set low values to 0 to [again] avoid floating underflow.
               dfsq = dblarr(size_cube[1],size_cube[2],nmod)
               dfsq = dfsq[i_no_under]
               i_no_under_2 = $
                  where(df_norm gt min AND df_wave gt min AND df_sig gt min,$
                        ct_no_under_2)
               if ct_no_under_2 gt 0 then $
                  dfsq[i_no_under_2] = df_norm[i_no_under_2]^2d + $
                                       df_wave[i_no_under_2]^2d + $
                                       df_sig[i_no_under_2]^2d
               emlcvdf['fluxerr',line,i_no_under] += dfsq
            endif
         endif
      endfor

      inz = where(emlcvdf['flux',line] gt 0,ctnz)
      if ctnz gt 0 then $
         emlcvdf['fluxerr',line,inz] = sqrt(emlcvdf['fluxerr',line,inz])

;     size of each model bin
      dmodwaves = modwaves[1:nmod-1] - modwaves[0:nmod-2]
;     ... rebinned to full cube.
      dmodwaves = rebin(reform([dmodwaves[0],dmodwaves],1,1,nmod),$
                               size_cube[1],size_cube[2],nmod)
;     Multiply flux density by bin size so that smaller bin sizes are 
;     downweighted. Otherwise cumulative distribution doesn't work.
      fluxnorm = emlcvdf['flux',line]*dmodwaves
      fluxnormerr = emlcvdf['fluxerr',line]*dmodwaves
;     Sum (renormalized) flux in wavelength, and then rebin to whole grid.
      fluxint = rebin(TOTAL(fluxnorm,3),size_cube[1],size_cube[2],nmod)
      inz = where(fluxint ne 0,ctnz)
      if ctnz gt 0 then begin
         fluxnorm[inz]/=fluxint[inz]
         fluxnormerr[inz]/=fluxint[inz]
      endif

      emlcvdf['cumfluxnorm',line,*,*,0] = fluxnorm[*,*,0]
      for i=1,nmod-1 do $
         emlcvdf['cumfluxnorm',line,*,*,i] = $
            emlcvdf['cumfluxnorm',line,*,*,i-1] + fluxnorm[*,*,i]
      emlcvdf['cumfluxnormerr',line,*,*,*] = fluxnormerr[*,*,*]

   endforeach
   
   return,emlcvdf

end
