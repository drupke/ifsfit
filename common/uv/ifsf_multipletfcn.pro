; docformat = 'rst'
;
;+
; This function is input to MPFIT, which uses it to compute the
; multiplet spectrum.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    An absorption specturm.
;
; :Params:
;    wave: in, required, type=dblarr(N)
;      Wavelengths
;    param: in, required, type=dblarr
;      Best-fit parameter array output by MPFIT.
;
; :Keywords:
;    cont: in, optional, type=dblarr(N)
;      Continuum used to originally normalize data for fits.
;    modabs: in, optional, type=dblarr(4,ndoubletabs)
;      Returns doublet absorption spectrum.
;    vels: in, optional, type=dblarr(N)
;      Velocity array w.r.t. the red, low-tau line of the doublet.
;    vwtabs: out, optional, type=dblarr(2)
;      The depth-weighted average velocity and velocity RMS.
;    weq: out, optional, type=structure
;      Returns equivalent widths. Structure tags is 'abs,' 
;      an array with the first element being the total equivalent width and 
;      subsequent elements being the equivalent widths of separate components.
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
;      2014may09, DSNR, created
;      2014may14, DSNR, removed floating underflow
;      2014jun10, DSNR, cleaned up treatment of optional individual component 
;                       outputs; added optional computation of Weq and emission
;                       line flux
;      2015xxxYY, AT, adapted for UV doublets
;      2016aug15, DSNR, code clean up for general doublet case
;      2018jun05, DSNR, added measurements of depth-weighted velocity and 
;                       RMS
;      2019nov20, DSNR, adapted from IFSF_DOUBLETFCN
;    
; :Copyright:
;    Copyright (C) 2014--2019 David S. N. Rupke, Anthony To
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
function ifsf_multipletfcn, wave,param,refloglf=refloglf,$
                            refmultwave=refmultwave,loglf=loglf,multwave=multwave,$
                            cont=cont,weq=weq,modabs=modabs, $
                            vels=vels,vwtabs=vwtabs,specres=specres,$
                            upsample=upsample

   tratio = dblarr(n_elements(loglf))
   lratio = dblarr(n_elements(loglf))
   for i=0,n_elements(loglf)-1 do begin
      dloglf = loglf[i] - refloglf
      tratio[i] = 10d^dloglf
      lratio[i] = refmultwave/multwave[i]
   endfor

   fwhmtosig = 2d*sqrt(2d*alog(2d))
   c = 299792.458d

;  Numbers of components
   nabs = param[0]

   nwave = n_elements(wave)
   dwave = wave[1:nwave-1] - wave[0:nwave-2]

;  Upsample wave
;  [Upsample must be odd.] This is the factor by which to "upsample" the
;  spectra (opposite of bin!) for the model computation. Problem for cases
;  where sigma << spectral resolution.
;  Factor here should be high enough that model contains say ~10 points within
;  mean +/- 2.5 sigma (99.8% of area). So if minimum possible sigma = 5 km/s,
;  upsample = 7, native dispersion = 0.44 (as for WiFeS), and a line at 6000 A,
;  then there is a point every 3 km/s, vs. 5-sigma = 25 km/s for the line.
   if keyword_set(upsample) then begin
      waveuse = rebin(reform(wave,nwave),nwave*upsample)
      nwaveuse = n_elements(waveuse)
;   waveuse = waveuse[0:nwaveuse-fac_upsample]
;   nwaveuse -= fac_upsample-1
      dwaveuse = waveuse[1:nwaveuse-1] - waveuse[0:nwaveuse-2]
;  Indices for downsampling back to original resolution
      dslo = ceil(double(upsample)/2d)
      dshi = nwaveuse-(1+floor(double(upsample)/2d))
   endif else begin
      waveuse=wave
      nwaveuse=nwave
      dwaveuse=dwave
      dslo=1
      dshi=nwave-1
   endelse


   modflux = dblarr(nwaveuse)+1d
;  absorption
   ilo = 1
   if nabs gt 0 then begin
      modabsuse = dblarr(nwaveuse,nabs)+1d
      modabs = dblarr(nwave,nabs)+1d
      modtmpuse = dblarr(nwaveuse)+1d
      modtmp = dblarr(nwave)+1d
      for i=0,nabs-1 do begin
         modabsuse[*,i]*=ifsf_cmpmultiplet(waveuse,param[ilo+i*4:ilo+(i+1)*4-1],$
                                        tratio,lratio)
         modtmpuse*=ifsf_cmpmultiplet(waveuse,param[ilo+i*4:ilo+(i+1)*4-1],$
                                      tratio,lratio)
         modtmpabs = rebin(modabsuse[dslo:dshi,i],nwave-1)
         modabs[*,i] = [modtmpabs[0],modtmpabs]
      endfor
      modflux *= modtmpuse
      modtmptmp = rebin(modtmpuse[dslo:dshi],nwave-1)
      modtmp = [modtmptmp[0],modtmptmp]
   endif else modabs=0

;  Optionally, compute equivalent widths
;  Note that we're doing this using the original binning ... and prior to convolution
   if keyword_set(weq) then begin
      if nabs gt 0 then begin
         absweq = dblarr(nabs+1)
         for i=0,nabs-1 do $
            absweq[i+1] = total((1d - modabs[1:nwave-1,i])*dwave)
;        Total equivalent width has to be calculated from full absorption
;        spectrum, rather than taken as the total of separate components.
         absweq[0] = total((1d - modtmp[1:nwave-1])*dwave)
      endif else absweq=0d
      weq = {abs: absweq}
   endif

   ;  Convolve model with line spread function, represented by a Gaussian
   if keyword_set(specres) then begin
      fwhm_pix = specres*fwhmtosig/mean(dwaveuse)
      npix = fix(fwhm_pix*5d)
      if not npix then npix++
      kernel = psf_gaussian(npix=npix,ndim=1,/double,fwhm=fwhm_pix)
      modflux_con = convol(modflux,kernel,/normalize,/edge_mirror)
   endif else $
      modflux_con = modflux
   
   if keyword_set(upsample) then begin
   ;  Downsample back to original resolution
      modflux_con_ds = rebin(modflux_con[dslo:dshi],nwave-1)
      modflux_con_ds = [modflux_con_ds[0],modflux_con_ds]
   endif else $
      modflux_con_ds = modflux_con

   ;  Weighted average velocity and weighted RMS width
   ;  Trump et al. 2006, ApJS, 165, 1
   ;  Section 4.5
   if keyword_set(vwtabs) AND keyword_set(vels) AND nabs gt 0 then begin
      ;     Compute singlet version of profile
      ilo = 1
      modtmp = dblarr(nwave)+1d
      for i=0,nabs-1 do $
         modtmp*=ifsf_cmpsinglet(wave,param[ilo+i*4:ilo+(i+1)*4-1],$
         tratio,lratio)
      igdvels = where(vels gt -2.9d4,ctgdvels)
      gdvels = vels[igdvels]
      gdmod = modtmp[igdvels]
      dgdvels = gdvels[1:ctgdvels-1] - gdvels[0:ctgdvels-2]
      AIabs = total((1d - gdmod[1:ctgdvels-1])*dgdvels)
      if AIabs gt 0d then begin
         vwtavgabs = total((1d - gdmod[1:ctgdvels-1])*dgdvels*gdvels)/AIabs
         vwtrmsabs = sqrt(total((1d - gdmod[1:ctgdvels-1])*dgdvels*(gdvels-vwtavgabs)^2)/AIabs)
      endif else begin
         vwtavgabs = 0d
         vwtrmsabs = 0d
      endelse
      vwtabs = [vwtavgabs,vwtrmsabs]
   endif
      
   return,modflux_con_ds
   
end
