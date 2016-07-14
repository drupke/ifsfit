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
;    Structure with 5 tags:
;      wave, type=dblarr(nwave), wavelengths of model
;      vel, type=dblarr(nwave), corresponding velocities w.r.t. systemic
;      flux, type=dblarr(ncols, nrows, nwave)
;      fluxnorm, type=dblarr(ncols, nrows, nwave), flux divided by total flux
;      cumfluxnorm, type=dblarr(ncols, nrows, nwave), cumulative fluxes at each 
;        wavelength, summed from low to high wavelengths, normalized so that
;        area under curve is 1.
;      cumfluxnorm, type=dblarr(ncols, nrows, nwave), cumulative fluxes at each
;        wavelength, summed from low to high wavelengths
;
; :Params:
;    pkfluxes: in, required, type=dblarr(ncols,nrows,maxncomp)
;      Peak fluxes of Gaussian emission-line fits.
;    pkwaves: in, required, type=dblarr(ncols,nrows,maxncomp)
;      Peak wavelengths of Gaussian emission-line fits.
;    sigmas: in, required, type=dblarr(ncols,nrows,maxncomp)
;      Sigmas of Gaussian emission-line fits.
;    maxncomp: in, required, type=integer
;      Maximum number of components.
;    waveref: in, required, type=double
;      Rest wavelength of emission line.
;    zref: in, required, type=double
;      Systemic redshift to which to refer output velocities.
;
; :Keywords:
;    ignore: in, optional, type=dblarr(maxncomp)
;      Set each element to 0 (don't ignore) or 1 (ignore) to force routine to
;      ignore a component.
;    ebv: in, optional, type=dblarr(ncols,nrows,maxncomp)
;      E(B-V) for each component, to de-redden Balmer line fluxes for that
;      component.
;      
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
;      2014jun03, DSNR, created
;      2016jul12, DSNR, added un-normalized, cumulative fluxes
;    
; :Copyright:
;    Copyright (C) 2014--2016 David S. N. Rupke
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
function ifsf_cmplinspecmaps,pkfluxes,pkwaves,sigmas,maxncomp,waveref,zref,$
                             ignore=ignore,ebv=ebv

   bad = 1d99
   c = 299792.458
   
   if ~ keyword_set(ignore) then ignore = dblarr(maxncomp)

   igd = where(pkwaves ne 0d AND pkwaves ne bad AND $
               pkfluxes ne 0d AND pkfluxes ne bad AND $
               sigmas ne 0d AND sigmas ne bad)
   modwaveran = [min(pkwaves[igd]*(1d - 5d*sigmas[igd]/c)),$
                 max(pkwaves[igd]*(1d + 5d*sigmas[igd]/c))]
   nwave = 1001
   modwaves = modwaveran[0]+dindgen(nwave)/double(nwave-1)*$
                            (modwaveran[1]-modwaveran[0])
   
   size_fluxes = size(pkfluxes)
   
   modfluxes = dblarr(size_fluxes[1],size_fluxes[2],nwave)
   for i=0,maxncomp-1 do begin
      if ~ ignore[i] then begin
         rbpkfluxes = rebin(pkfluxes[*,*,i],size_fluxes[1],size_fluxes[2],nwave)
         rbpkwaves = rebin(pkwaves[*,*,i],size_fluxes[1],size_fluxes[2],nwave)
         rbsigmas = rebin(sigmas[*,*,i],size_fluxes[1],size_fluxes[2],nwave)
         rbmodwaves = rebin(reform(modwaves,1,1,nwave),$
                            size_fluxes[1],size_fluxes[2],nwave)
         inz = where(rbsigmas gt 0 AND rbsigmas ne bad AND $
                     rbpkfluxes gt 0 AND rbpkfluxes ne bad,ctnz)
         if ctnz gt 0 then begin
            exparg = dblarr(size_fluxes[1],size_fluxes[2],nwave)
            exparg[inz] = ((rbmodwaves[inz]/rbpkwaves[inz] - 1d) / $
                           (rbsigmas[inz]/c))^2d
;           Code from GAUSSIAN function in ASTROLIB: Get smallest value expressible 
;           on computer. Set lower values to 0 to avoid floating underflow.
            min = (machar(/double)).xmin
            minexp = alog(min)
            inz = where(exparg LT -2d*minexp AND $
                        rbpkfluxes GT min AND $
                        rbpkfluxes ne bad, ctnz)
            if ctnz gt 0 then begin
               cmodfluxes = dblarr(size_fluxes[1],size_fluxes[2],nwave)
               cmodfluxes[inz] = rbpkfluxes[inz]*exp(-exparg[inz]/2d)
            endif
            if keyword_set(ebv) then begin
               cebv = ebv[*,*,i]
               iebv = where(cebv ne bad,ctebv)
               for j=0,nwave-1 do begin
                  cmodfluxes_tmp = cmodfluxes[*,*,j]
                  cmodfluxes_dr = dblarr(size_fluxes[1],size_fluxes[2])
                  cmodfluxes_dr[iebv] = $
                     ifsf_dustcor_ccm(waveref,cmodfluxes_tmp[iebv],cebv[iebv])
                  cmodfluxes[*,*,j] = cmodfluxes_dr
               endfor
            endif
            if ctnz gt 0 then modfluxes[inz] += cmodfluxes[inz]
         endif
      endif
   endfor
   
;  TOTAL(modfluxes,3) gives model flux summed over wavelength; REBIN then grids
;  this number back over wavelength for dividing into MODFLUXES. MODFLUXES_NORM
;  is a somewhat ambiguous quantity, since it then depends on the binning of the 
;  model. CUMFLUXES_NORM, however, doesn't depend on the binning; it
;  regardless represents the percent area that has been covered as one moves 
;  through the distribution.
   totalfluxes = rebin(total(modfluxes,3),size_fluxes[1],size_fluxes[2],nwave)
   modfluxes_norm = dblarr(size_fluxes[1],size_fluxes[2],nwave)
   inz = where(totalfluxes ne 0,ctnz)
   if ctnz gt 0 then modfluxes_norm[inz] = modfluxes[inz] / totalfluxes[inz]

   cumfluxes = dblarr(size_fluxes[1],size_fluxes[2],nwave)
   cumfluxes_norm = dblarr(size_fluxes[1],size_fluxes[2],nwave)
   cumfluxes_norm[*,*,0] = modfluxes_norm[*,*,0]
   for i=1,nwave-1 do $
      cumfluxes_norm[*,*,i] = cumfluxes_norm[*,*,i-1] + modfluxes_norm[*,*,i]

;  relativistic velocity shift;
;  see http://hyperphysics.phy-astr.gsu.edu/hbase/relativ/reldop2.html
   zdiff = modwaves/waveref-1d - zref
   modvels = c * ((zdiff+1d)^2d - 1d) / ((zdiff+1d)^2d + 1d)

   return,{wave: modwaves,$
           vel: modvels,$
           flux: modfluxes,$           
           cumfluxnorm: cumfluxes_norm}

end
