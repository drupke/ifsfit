; docformat = 'rst'
;
;+
; This function is input to MPFIT, which uses it to compute the
; doublet spectrum.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    A doublet absorption spectrum.
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
;    doubletname: in, optional, type=string
;      Name of doublet
;    doubletemflux: in, optional, type=dblarr(ndoubletem+1)
;      Returns flux of doublet emission. The first element is the total flux and
;      subsequent elements are the fluxes of separate components. Requires
;      continuum, as well.
;    moddoubletabs: in, optional, type=dblarr(4,ndoubletabs)
;      Returns doublet absorption spectrum.
;    moddoubletem: in, optional, type=dblarr(4,ndoubletem)
;      Returns doublet emission spectrum.
;    vels: in, optional, type=dblarr(N)
;      Velocity array w.r.t. the red, low-tau line of the doublet.
;    vwtabs: out, optional, type=dblarr(2)
;      The depth-weighted average velocity and velocity RMS.
;    weq: out, optional, type=structure
;      Returns equivalent widths. Structure tags are 'em' and 'abs,' and each is
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
;    
; :Copyright:
;    Copyright (C) 2014--2018 David S. N. Rupke, Anthony To
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
function ifsf_doubletfcn, wave, param,$
                          cont=cont,doubletname=doubletname,weq=weq,$
                          modabs=modabs, $
                          modem=modem, $
                          emflux=emflux, $
                          vels=vels,vwtabs=vwtabs, $
                          specres=specres,upsample=upsample,$
                          cfcorr=cfcorr

;  Default to NV unless otherwise specified
   if keyword_set(doubletname) then begin
      if doubletname eq 'NV' then begin
;        Data from Morton 2003
;        optical depth ratio
         LogLF = 2.286d -1.985d
         tratio = 10d^LogLF
;        corresponding wavelength ratio
         lratio = 1242.804d/1238.821d
      endif else if doubletname eq 'OVI' then begin
;        Data from Morton 2003
         LogLF = 2.136d -1.834d
         tratio = 10d^LogLF
         lratio = 1037.6167d/1031.9261d
      endif else if doubletname eq 'PV' then begin
;        Data from Morton 2003
         LogLF = 2.722d -2.420d
         tratio = 10d^LogLF
         lratio = 1128.0078d/1117.9774d
      endif else if doubletname eq 'MgII' then begin
;        Data from Morton 2003
         LogLF = 3.236d - 2.933d
         tratio = 10d^LogLF
         lratio = 2803.530d/2796.352d
         endif else if doubletname eq 'CaII' then begin
;        Data from Morton 1991
         LogLF = 3.397d - 3.096d
         tratio = 10d^LogLF
         lratio = 3969.591d/3934.777d
      endif else if doubletname eq 'FeIIUV1' then begin
;        Data from Morton 2003
         LogLF = 2.793d - 2.252d
         tratio = 10d^LogLF
         lratio = 2585.8762d/2599.3959d
      endif else if doubletname eq 'FeIIUV2' then begin
;        Data from Morton 2003
         LogLF = 2.882d - 1.871d
         tratio = 10d^LogLF
         lratio = 2373.7365d/2382.0386d
      endif
   endif else begin
      doubletname = 'NV'
;     Data from Morton 2003
      LogLF = 2.286d -1.985d
      tratio = 10d^LogLF
      lratio = 1242.804d/1238.821d
   endelse

   c = 299792.458d
   fwhmtosig = 2d*sqrt(2d*alog(2d))

;  Numbers of components
   nabs = param[0]
   nem = param[1]

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
      dslo=1 ; come back to this;  should be 0 or 1, not sure
      dshi=nwave-1 ; come back to this; should be nwave-2 or nwave-1, not sure
   endelse

   modflux = dblarr(nwaveuse)+1d
;  doublet emission
   ilo = 2+nabs*4
   if nem gt 0 then begin
      modemuse=dblarr(nwaveuse,nem)+1d
      modem=dblarr(nwave,nem)+1d
   endif else modem=0
   for i=0,nem-1 do begin
      arg1 = ((waveuse-param[ilo+i*4])/$
             (2d*param[ilo+i*4]*param[1+ilo+i*4]/c))^2d
      arg2 = ((lratio * waveuse-param[ilo+i*4])/$
             (2d*param[ilo+i*4]*param[1+ilo+i*4]/c))^2d
      mask1 = (arg1 lt 80)
      mask2 = (arg2 lt 80)
      modemuse[*,i]+=param[2+ilo+i*4]*(mask1*exp(-arg1*mask1) + $
         param[3+ilo+i*4]*mask2*exp(-arg2*mask2))
      modflux+=param[2+ilo+i*4]*(mask1*exp(-arg1*mask1) + $
         param[3+ilo+i*4]*mask2*exp(-arg2*mask2))      
      modtmpem = rebin(modemuse[dslo:dshi,i],nwave-1)
      modem[*,i] = [modtmpem[0],modtmpem]
   endfor
;  doublet absorption
   ilo = 2
   if nabs gt 0 then begin
      if keyword_set(cfcorr) then begin
         ccmodabsuse = dblarr(nwaveuse,nabs)+1d
         ccmodabs = dblarr(nwave,nabs)+1d
         ccmodtmpuse = dblarr(nwaveuse)+1d
         ccmodtmp = dblarr(nwave)+1d
      endif
      modabsuse = dblarr(nwaveuse,nabs)+1d
      modabs = dblarr(nwave,nabs)+1d
      modtmpuse = dblarr(nwaveuse)+1d
      modtmp = dblarr(nwave)+1d
      for i=0,nabs-1 do begin
         modabsuse[*,i]*=ifsf_cmpdoublet(waveuse,param[ilo+i*4:ilo+(i+1)*4-1],$
            tratio,lratio)
         modtmpuse*=ifsf_cmpdoublet(waveuse,param[ilo+i*4:ilo+(i+1)*4-1],$
            tratio,lratio)
         modtmpabs = rebin(modabsuse[dslo:dshi,i],nwave-1)
         modabs[*,i] = [modtmpabs[0],modtmpabs]
         if keyword_set(cfcorr) then begin
            ccmodabsuse[*,i]*=ifsf_cmpdoublet(waveuse,$
               [1d,param[1+ilo+i*4:ilo+(i+1)*4-1]],$
               tratio,lratio)
            ccmodtmpuse*=ifsf_cmpdoublet(waveuse,$
               [1d,param[1+ilo+i*4:ilo+(i+1)*4-1]],$
               tratio,lratio)
            ccmodtmpabs = rebin(ccmodabsuse[dslo:dshi,i],nwave-1)
            ccmodabs[*,i] = [ccmodtmpabs[0],ccmodtmpabs]
         endif
      endfor
      modflux *= modtmpuse
      modtmptmp = rebin(modtmpuse[dslo:dshi],nwave-1)
      modtmp = [modtmptmp[0],modtmptmp]
      if keyword_set(cfcorr) then begin
         ccmodtmptmp = rebin(ccmodtmpuse[dslo:dshi],nwave-1)
         ccmodtmp = [ccmodtmptmp[0],ccmodtmptmp]
      endif
   endif else modabs=0

;  Optionally, compute equivalent widths
;  Note that we're doing this using the original binning ... and prior to convolution
   if keyword_set(weq) then begin
      if nem gt 0 then begin
         emweq = dblarr(nem+1)
         for i=0,nem-1 do $
            emweq[i+1] = total((1d - modem[1:nwave-1,i])*dwave)
         emweq[0] = total(emweq[1:nem])
      endif else emweq=0d
      if keyword_set(cfcorr) then begin
         modabsuse = ccmodabs
         modtmpuse = ccmodtmp
      endif else begin
         modabsuse = modabs
         modtmpuse = modtmp
      endelse
      if nabs gt 0 then begin
         absweq = dblarr(nabs+1)
         for i=0,nabs-1 do $
            absweq[i+1] = total((1d - modabsuse[1:nwave-1,i])*dwave)
;        Total doublet equivalent width has to be calculated from full absorption
;        spectrum, rather than taken as the total of separate components.
         absweq[0] = total((1d - modtmpuse[1:nwave-1])*dwave)
      endif else absweq=0d
      weq = {em: emweq,abs: absweq}
   endif
   if keyword_set(emflux) AND keyword_set(cont) then begin
      if nem gt 0 then begin
         emflux = dblarr(nem+1)
         for i=0,nem-1 do $
            emflux[i+1] = $
               total((modem[1:nwave-1,i]*cont[1:nwave-1]-cont[1:nwave-1])*dwave)
         emflux[0] = total(emflux[1:nem])
      endif else emflux=0d
   endif else emflux=0d

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

   ;  Weighted absorption average velocity and weighted RMS width
   ;  Note that we're doing this using the original binning ... and w/o convolution
   ;  Trump et al. 2006, ApJS, 165, 1
   ;  Section 4.5
   if keyword_set(vwtabs) AND keyword_set(vels) AND nabs gt 0 then begin
      ;     Compute singlet version of profile
      ilo = 2
      modtmp = dblarr(nwave)+1d
      for i=0,nabs-1 do begin
         if keyword_set(cfcorr) then paruse = [1d,param[1+ilo+i*4:ilo+(i+1)*4-1]] $
         else paruse = param[ilo+i*4:ilo+(i+1)*4-1]
         modtmp*=ifsf_cmpsinglet(wave,paruse,tratio,lratio)
      endfor
      ;      igdvels = where(vels lt 0 AND vels gt -2.9d4,ctgdvels)
      igdvels = where(vels gt -2.9d4,ctgdvels)
      ;      gdvels = abs(vels[igdvels])
      gdvels = vels[igdvels]
      gdmod = modtmp[igdvels]
      ;      dgdvels = gdvels[0:ctgdvels-2] - gdvels[1:ctgdvels-1]
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
