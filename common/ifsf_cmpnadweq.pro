; docformat = 'rst'
;
;+
;
; Compute absorption and emission-line equivalent widths for NaD.
; 
; This algorithm will detect emission and absorption line regions regardless of
; their relative position in wavelength space. However, it will only give correct
; results for the equivalent width if the emission and absorption-line regions
; do not overlap; otherwise any excess emission will be included in the absorption
; line equivalent width calculation, and vice versa.
;
; Weq is set to 0 if AUTOWAVELIM is set and no lines are detected. Weq(EM) is set
; to 0 if neither WAVELIM nor AUTOWAVELiM is set.
;
; :Categories:
;    IFSFIT
;
; :Returns: 
;    Array of equivalent widths and their errors. If AUTOWAVELIM is
;    selected, also returns the indices into the input arrays that
;    give the lower and upper boundaries of the absorption and
;    emission lines.
;
; :Params:
;    wave: in, required, type=dblarr(N)
;      Wavelengths.
;    flux: in, required, type=dblarr(N)
;      Normalized fluxes.
;    err: in, required, type=dblarr(N)
;      Flux errors.
;
; :Keywords:
;    wavelim: in, optional, type=dblarr(4)
;      Limits of integration for equivalent widths. First two elements are 
;      wavelength ranges of absorption, second two are for emission. Only applies
;      if autowavelim not set
;    autowavelim: in, optional, type=dblarr(4)
;      Same as wavelim, except now the ranges are for finding the actual
;      absorption and emission limits using the automatic algorithm. If this 
;      keyword is selected, the output array has a second dimension
;      containing the indices into the input arrays that give the
;      lower and upper boundaries of the absorption and emission lines.
;    smoothkernel: in, optional, type=long, default=5
;      Kernel size for boxcar smoothing the spectrum to look for features. The
;      default is optimized for the GMOS B600 grating, which has a spectral 
;      resolution element of ~5 pixels at 6000 A (based on 0.46 A dispersion 
;      and interpolating blaze resolutions of 1688 and 3744 at 461 and 926 nm, 
;      respectively).
;    emflux: out, optional, type=dblarr(2)
;      Emission line flux and error. The former also requires the subtraction-
;      normalized spectrum (SNFLUX) and the latter the un-normalized flux
;      (UNERR).
;    emul: out, optional, type=dblarr(4)
;      Estimates for upper limits to emission line equivalent width, error, 
;      flux, and error. If absorption line is detected but no emission line, 
;      integrate over specified wavelength range to get upper limits.
;    emwid: in, optional, type=double, default=15d
;      Wavelength range over which to integrate to estimate upper limit for 
;      emission line equivalent width and flux.
;    iabsoff: in, optional, type=long, default=4l
;      Index offset from absorption line for calculating emission line upper limit.
;    snflux: in, optional, type=dblarr(N)
;      Subtraction-normalized flux array.
;    snrabsthresh: in, optional, type=double, default=1.6
;      Threshold S/N per pixel for detecting a pixel that can be attributable 
;      to possible absorption.
;    snremthresh: in, optional, type=double, default=1.6
;    unerr: in, optional, type=dblarr(N)
;      Un-normalized error array.
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
;      2014may14, DSNR, created
;      2014junXY, DSNR, added ability to output emission line flux
;      2014jun18, DSNR, added ability to output emission line flux upper limits
;      2014jul29, DSNR, updated input keywords for emission line limits
;      2016may03, DSNR, added trigger to deal with very noisy data or data
;                       where continuum goes below 0
;      2021jan07, DSNR, updated documentation
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
function ifsf_cmpnadweq,wave,flux,err,$
                        autowavelim=autowavelim,emflux=emflux,emul=emul,$
                        smoothkernel=smoothkernel,snflux=snflux,unerr=unerr,$
                        wavelim=wavelim,emwid=emwid,snrabsthresh=snrabsthresh,$
                        snremthresh=snremthresh

;  Thresholds for line detection   

   if ~ keyword_set(snrabsthresh) then snrabsthresh=1.6d
   if ~ keyword_set(snremthresh) then snremthresh=-1.6d
   
   if ~ keyword_set(smoothkernel) then smoothkernel=5l
   if ~ keyword_set(iabsoff) then iabsoff = 4l
   if ~ keyword_set(emwid) then emwid = 20d 

;  Case if IFSF_NORMNAD decided the data was bad   
   inz = where(flux gt 0,ctnz)
   if ctnz eq 0 then begin
      out = [0d,0d,0d,0d]
      if keyword_set(autowavelim) then begin
         autoindices=[0d,0d,0d,0d]
         out=[[out],[autoindices]]
      endif
      goto,baddata
   endif
   
;  If neither WAVELIM nor AUTOWAVELIM is set, then integration defaults to entire wavelength 
;  range and absorption only.
   if ~ keyword_set(wavelim) then begin
      iabslo = 1l
      iabsup = n_elements(wave)-1l
      iemlo = -1l
      iemup = -1l
;  If WAVELILM is set but AUTOWAVELIM is not set, then use these limits
   endif else if ~ keyword_set(autowavelim) then begin
      iabslo = value_locate(wave,wavelim[0])
      iabsup = value_locate(wave,wavelim[1])
      iemlo = value_locate(wave,wavelim[2])
      iemup = value_locate(wave,wavelim[3])
   endif

;  If AUTOWAVELIM is set:
;  Algorithm for automatically finding wavelength ranges for absorption and 
;  emission line integration. The algorithm boxcar smooths the spectrum and 
;  error using a kernel equal to the size of a spectral resolution element. This
;  is akin to averaging -smoothkernel- adjacent points. If this average yields a point with
;  SNR > snrabsthresh, then it's added to the list of absorption line indices.
;  Same for emission, but SNR < snremthresh.
   if keyword_set(autowavelim) then begin
;     subtract/add 2 pixels to avoid edge effects
      i1abs = value_locate(wave,autowavelim[0])-2l
      i2abs = value_locate(wave,autowavelim[1])+2l
      i1em = value_locate(wave,autowavelim[2])-2l
      i2em = value_locate(wave,autowavelim[3])+2l
      sflux = smooth(flux,smoothkernel,/edge_truncate)
      serr = sqrt(smooth(err^2d,smoothkernel,/edge_truncate))
;     1 minus flux, since we care about features above/below the continuum
      snr = (1d -sflux)/serr
;     truncate smoothed arrays for easier bookkeeping
      sflux_abs=sflux[i1abs:i2abs]
      serr_abs=serr[i1abs:i2abs]
      snr_abs=snr[i1abs:i2abs]
      sflux_em=sflux[i1em:i2em]
      serr_em=serr[i1em:i2em]
      snr_em=snr[i1em:i2em]
;     find significant absorption/emission features
      iabs = where(snr_abs ge snrabsthresh,ctabs)
      iem = where(snr_em le snremthresh,ctem)
      ; Make sure more than one point is found, and then assign lower and upper
      ; indices of range based on lowest wavelength and highest wavelength points
      ; found. The extra bit of logic makes sure that the indices don't stray
      ; outside the ranges of the truncated arrays.
      if ctabs gt 1 then begin
         iabslo = iabs[0]-2l lt 1 ? 1 : iabs[0]-2l
         iabsup = $
            iabs[ctabs-1]+2l gt i2abs-i1abs ? i2abs-i1abs : iabs[ctabs-1]+2l
         iabslo += i1abs
         iabsup += i1abs
      endif else begin
         iabslo=-1l
         iabsup=-1l
      endelse
;     Same for emission.
      if ctem gt 1 then begin
         iemlo=iem[0]-2l lt 1 ? 1 : iem[0]-2l
         iemup = iem[ctem-1]+2l gt i2em-i1em ? i2em-i1em : iem[ctem-1]+2l
         iemlo += i1em
         iemup += i1em
      endif else begin
         iemlo=-1l
         iemup=-1l
      endelse
   endif else begin
      ctabs=0
      ctem=0
   endelse
      
   ;  Compute equivalent widths.
   ;  Absorption equivalent width
   if iabslo ne -1l AND iabsup ne -1l then begin
      weq_abs = total((1d -flux[iabslo:iabsup])*$
                      (wave[iabslo:iabsup]-$
                       wave[iabslo-1:iabsup-1]))
      weq_abs_e = sqrt(total(err[iabslo:iabsup]^2*$
                             (wave[iabslo:iabsup]-$
                              wave[iabslo-1:iabsup-1])))
   endif else begin
      weq_abs = 0d
      weq_abs_e = 0d
   endelse
   ;  Emission equivalent width
   if iemlo ne -1l AND iemup ne -1l then begin
      weq_em = total((1d -flux[iemlo:iemup])*$
                     (wave[iemlo:iemup]-wave[iemlo-1:iemup-1]))
      weq_em_e = sqrt(total(err[iemlo:iemup]^2*$
                            (wave[iemlo:iemup]-wave[iemlo-1:iemup-1])))
      if keyword_set(snflux) AND keyword_set(emflux) then begin
         fl_em = total(snflux[iemlo:iemup]*$
                       (wave[iemlo:iemup]-wave[iemlo-1:iemup-1]))
         if keyword_set(unerr) then $
            fl_em_e = sqrt(total(unerr[iemlo:iemup]^2*$
                                 (wave[iemlo:iemup]-wave[iemlo-1:iemup-1]))) $
         else fl_em_e = 0d
         emflux = [fl_em,fl_em_e]
      endif
   ;  Emission equivalent width upper limit
   endif else if keyword_set(emul) AND ctabs gt 1 then begin
      weq_em = 0d
      weq_em_e = 0d
      emul=dblarr(4)
      ilo = iabsup+iabsoff
      iup = value_locate(wave,wave[ilo]+emwid)
      emul[0] = total((1d -flux[ilo:iup])*$
                      (wave[ilo:iup]-wave[ilo-1:iup-1]))
      emul[1] = sqrt(total(err[ilo:iup]^2*$
                           (wave[ilo:iup]-wave[ilo-1:iup-1])))
      if keyword_set(snflux) AND keyword_set(emflux) then begin
         emul[2] = total(snflux[ilo:iup]*$
                         (wave[ilo:iup]-wave[ilo-1:iup-1]))
         if keyword_set(unerr) then $
            emul[3] = sqrt(total(unerr[ilo:iup]^2*$
                                 (wave[ilo:iup]-wave[ilo-1:iup-1])))
      endif
   endif else begin
      weq_em = 0d
      weq_em_e = 0d
   endelse

   out = [weq_abs,weq_abs_e,weq_em,weq_em_e]
;  indices into original arrays of auto-detect regions
   if keyword_set(autowavelim) then begin
      autoindices=[iabslo,iabsup,iemlo,iemup]
      out=[[out],[autoindices]]
   endif

baddata:

   return,out
  
end
