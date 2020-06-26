; docformat = 'rst'
;
;+
; This function is input to MPFIT, which uses it to compute the
; NaD + HeI spectrum.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    An NaD + HeI spectrum.
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
;    modhei: in, optional, type=dblarr(3,nhei)
;      Returns HeI emission spectrum.
;    modnadabs: in, optional, type=dblarr(4,nnadabs)
;      Returns NaD absorption spectrum.
;    modnadem: in, optional, type=dblarr(4,nnadem)
;      Returns NaD emission spectrum.
;    nademflux: in, optional, type=dblarr(nnadem+1)
;      Returns flux of NaD emission. The first element is the total flux and
;      subsequent elements are the fluxes of separate components. Requires 
;      continuum, as well.
;    specres: in, optional, type=double, def=0.64d
;      Estimated spectral resolution in wavelength units (sigma).
;    wavhei: out, optional, type=dblarr(nhei)
;      Wavelengths of HeI emission lines in obs. frame
;    wavnadabs: out, optional, type=dblarr(nnadabs)
;      Wavelengths of NaD absorption lines in obs. frame (NaD_1, red line 5896)
;    wavnadem: out, optional, type=dblarr(nnadem)
;      Wavelengths of NaD emission lines in obs. frame (NaD_1, red line 5896)
;    weq: in, optional, type=structure
;      Returns equivalent widths. Structure tags are 'em' and 'abs,' and each is
;      an array with the first element being the total equivalent width and
;      subsequent elements being the equivalent widths of separate components.
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
;      2014may09, DSNR, created
;      2014may14, DSNR, removed floating underflow
;      2014jun10, DSNR, cleaned up treatment of optional individual component 
;                       outputs; added optional computation of Weq and emission
;                       line flux
;      2016nov03, DSNR, added convolution with spectral resolution
;      2017may18, DSNR, upsample spectra to avoid undersampling model in case of
;                       v. narrow components
;      2020jun25, DSNR, option to output lambda(NaD em)
;    
; :Copyright:
;    Copyright (C) 2014--2020 David S. N. Rupke
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
function ifsf_nadfcn, wave, param, modhei=modhei, modnadabs=modnadabs, $
                      modnadem=modnadem, weq=weq, nademflux=nademflux, $
                      cont=cont, specres=specres, wavnadabs=wavnadabs, $
                      wavnadem=wavnadem,wavhei=wavhei



   c = 299792.458d
   fwhmtosig = 2d*sqrt(2d*alog(2d))
;  NaD wavelength ratio (red to blue)
   lratio = 1.001014158d

;  Numbers of components
   nhei = param[0]
   nnadabs = param[1]
   nnadem = param[2]

;  Upsample wave
   nwave = n_elements(wave)
;  [Upsample must be odd.] This is the factor by which to "upsample" the 
;  spectra (opposite of bin!) for the model computation. Problem for cases 
;  where sigma << spectral resolution.
;  Factor here should be high enough that model contains say ~10 points within
;  mean +/- 2.5 sigma (99.8% of area). So if minimum possible sigma = 5 km/s,
;  upsample = 7, native dispersion = 0.44 (as for WiFeS), and a line at 6000 A,
;  then there is a point every 3 km/s, vs. 5-sigma = 25 km/s for the line.
   fac_upsample = 7
   waveuse = rebin(reform(wave,nwave),nwave*fac_upsample)
   nwaveuse = n_elements(waveuse)
;   waveuse = waveuse[0:nwaveuse-fac_upsample]
;   nwaveuse -= fac_upsample-1

   dwave = wave[1:nwave-1] - wave[0:nwave-2]
   dwaveuse = waveuse[1:nwaveuse-1] - waveuse[0:nwaveuse-2]

;  Indices for downsampling back to original resolution
   dslo = ceil(double(fac_upsample)/2d)
   dshi = nwaveuse-(1+floor(double(fac_upsample)/2d))

;  HeI emission
   modflux = dblarr(nwaveuse)+1d
   if nhei gt 0 then begin
      wavhei=dblarr(nhei)
      modheiuse=dblarr(nwaveuse,nhei)+1d
      modhei=dblarr(nwave,nhei)+1d
      for i=0,nhei-1 do begin
         wavhei[i]=param[3+i*3]
;         sigma = sqrt((param[3+i*3]*param[4+i*3]/c)^2d + specres^2d)
         sigma = sqrt((param[3+i*3]*param[4+i*3]/c)^2d)
         arg = ((waveuse-param[3+i*3])/(2d*sigma))^2d
         mask = (arg lt 80)
         modheiuse[*,i]+=param[5+i*3]*mask*exp(-arg*mask)
         modflux+=param[5+i*3]*mask*exp(-arg*mask)
         modtmphei = rebin(modheiuse[dslo:dshi,i],nwave-1)
         modhei[*,i] = [modtmphei[0],modtmphei]
      endfor
   endif else begin
      wavhei=0
      modhei=0
   endelse
;  NaD emission
   ilo = 3+nhei*3+nnadabs*4
   if nnadem gt 0 then begin
      modnademuse=dblarr(nwaveuse,nnadem)+1d
      modnadem=dblarr(nwave,nnadem)+1d
      wavnadem = dblarr(nnadem)
   endif else begin
      wavnadem=0
      modnadem=0d
   endelse
   for i=0,nnadem-1 do begin
      wavnadem[i] = param[ilo+i*4]
;      sigma = sqrt((param[ilo+i*4]*param[1+ilo+i*4]/c)^2d + specres^2d)
      sigma = sqrt((param[ilo+i*4]*param[1+ilo+i*4]/c)^2d)
      arg1 = ((waveuse-param[ilo+i*4])/(2d*sigma))^2d
      arg2 = ((lratio * waveuse-param[ilo+i*4])/(2d*sigma))^2d
      mask1 = (arg1 lt 80)
      mask2 = (arg2 lt 80)
      modnademuse[*,i]+=param[2+ilo+i*4]*(mask1*exp(-arg1*mask1) + $
                                          param[3+ilo+i*4]*mask2*exp(-arg2*mask2))
      modflux+=param[2+ilo+i*4]*(mask1*exp(-arg1*mask1) + $
                                 param[3+ilo+i*4]*mask2*exp(-arg2*mask2))
      modtmpnadem = rebin(modnademuse[dslo:dshi,i],nwave-1)
      modnadem[*,i] = [modtmpnadem[0],modtmpnadem]
   endfor
;  NaD absorption
   ilo = 3+nhei*3
   if nnadabs gt 0 then begin
      modnadabsuse = dblarr(nwaveuse,nnadabs)+1d
      modnadabs = dblarr(nwave,nnadabs)+1d
      wavnadabs = dblarr(nnadabs)
      modtmpuse = dblarr(nwaveuse)+1d
      modtmp = dblarr(nwave)+1d
      for i=0,nnadabs-1 do begin
         wavnadabs[i] = param[ilo+i*4+2]
         modnadabsuse[*,i]*=ifsf_cmpnad(waveuse,param[ilo+i*4:ilo+(i+1)*4-1])
         modtmpuse*=ifsf_cmpnad(waveuse,param[ilo+i*4:ilo+(i+1)*4-1])
         modtmpnadabs = rebin(modnadabsuse[dslo:dshi,i],nwave-1)
         modnadabs[*,i] = [modtmpnadabs[0],modtmpnadabs]
      endfor
      modflux *= modtmpuse
      modtmptmp = rebin(modtmpuse[dslo:dshi],nwave-1)
      modtmp = [modtmptmp[0],modtmptmp]
   endif else begin
      wavnadabs=0
      modnadabs=0
   endelse

;  Optionally, compute equivalent widths
   if keyword_set(weq) then begin
      if nnadem gt 0 then begin
         emweq = dblarr(nnadem+1)
         for i=0,nnadem-1 do $
            emweq[i+1] = total((1d - modnadem[1:nwave-1,i])*dwave)
         emweq[0] = total(emweq[1:nnadem])
      endif else emweq=0d
      if nnadabs gt 0 then begin
         absweq = dblarr(nnadabs+1)
         for i=0,nnadabs-1 do $
            absweq[i+1] = total((1d - modnadabs[1:nwave-1,i])*dwave)
;        Total NaD equivalent width has to be calculated from full absorption
;        spectrum, rather than taken as the total of separate components.
         absweq[0] = total((1d - modtmp[1:nwave-1])*dwave)
      endif else absweq=0d
      weq = {em: emweq,abs: absweq}
   endif
   if keyword_set(nademflux) AND keyword_set(cont) then begin
      if nnadem gt 0 then begin
         nademflux = dblarr(nnadem+1)
         for i=0,nnadem-1 do $
            nademflux[i+1] = $
               total((modnadem[1:nwave-1,i]*cont[1:nwave-1]-cont[1:nwave-1])*dwave)
         nademflux[0] = total(nademflux[1:nnadem])
      endif else nademflux=0d
   endif else nademflux=0d

;  Convolve model with line spread function, represented by a Gaussian
   fwhm_pix = specres*fwhmtosig/mean(dwaveuse)
   npix = fix(fwhm_pix*5d)
   if not npix then npix++
   kernel = psf_gaussian(npix=npix,ndim=1,/double,fwhm=fwhm_pix)
   modflux_con = convol(modflux,kernel,/normalize,/edge_mirror)
   
;  Downsample back to original resolution
   modflux_con_ds = rebin(modflux_con[dslo:dshi],nwave-1)
   modflux_con_ds = [modflux_con_ds[0],modflux_con_ds]
   
   return,modflux_con_ds

end
