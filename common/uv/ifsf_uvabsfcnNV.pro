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
;    modhei: in, optional, type=dblarr(3,nhei)
;      Returns HeI emission spectrum.
;    modnadabs: in, optional, type=dblarr(4,nnadabs)
;      Returns NaD absorption spectrum.
;    modnadem: in, optional, type=dblarr(4,nnadem)
;      Returns NaD emission spectrum.
;    weq: in, optional, type=structure
;      Returns equivalent widths. Structure tags are 'em' and 'abs,' and each is
;      an array with the first element being the total equivalent width and 
;      subsequent elements being the equivalent widths of separate components.
;    nademflux: in, optional, type=dblarr(nnadem+1)
;      Returns flux of NaD emission. The first element is the total flux and
;      subsequent elements are the fluxes of separate components. Requires 
;      continuum, as well.
;    cont: in, optional, type=dblarr(N)
;      Continuum used to originally normalize data for fits.
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
;    
; :Copyright:
;    Copyright (C) 2014 David S. N. Rupke
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
function ifsf_uvabsfcnNV, wave, param, lratio, tratio, modhei=modhei, modnadabs=modnadabs, $
                      modnadem=modnadem, weq=weq, nademflux=nademflux, $
                      cont=cont
   c = 299792.458d
;  NV wavelength ratio (red to blue)
   lratio = 1242.804/1238.821d
;  OVI wavelength ratio (red to blue)  
;   lratio = 1037.613/1031.912d


;  Numbers of components
;   nhei = param[0]
   nnadabs = param[0]
   nnadem = param[1]

   nwave = n_elements(wave)

;  HeI emission
   modflux = dblarr(nwave)+1d
;   if nhei gt 0 then modhei=dblarr(nwave,nhei)+1d else modhei=0
;   for i=0,nhei-1 do begin
;      arg = ((wave-param[3+i*3])/$
;            (2d*param[3+i*3]*param[4+i*3]/c))^2d
;      mask = (arg lt 80)
;      modhei[*,i]+=param[5+i*3]*mask*exp(-arg*mask)
;      modflux+=param[5+i*3]*mask*exp(-arg*mask)
;   endfor
;  NaD emission
;   ilo = 3+nhei*3+nnadabs*4
   ilo = 2+nnadabs*4
   if nnadem gt 0 then modnadem=dblarr(nwave,nnadem)+1d else modnadem=0
   for i=0,nnadem-1 do begin
      arg1 = ((wave-param[ilo+i*4])/$
             (2d*param[ilo+i*4]*param[1+ilo+i*4]/c))^2d
      arg2 = ((lratio * wave-param[ilo+i*4])/$
             (2d*param[ilo+i*4]*param[1+ilo+i*4]/c))^2d
      mask1 = (arg1 lt 80)
      mask2 = (arg2 lt 80)
      modnadem[*,i]+=param[2+ilo+i*4]*(mask1*exp(-arg1*mask1) + $
                                       param[3+ilo+i*4]*mask2*exp(-arg2*mask2))
      modflux+=param[2+ilo+i*4]*(mask1*exp(-arg1*mask1) + $
                                 param[3+ilo+i*4]*mask2*exp(-arg2*mask2))      
   endfor
;  NaD absorption
;   ilo = 3+nhei*3
   ilo = 2
   if nnadabs gt 0 then begin
      modnadabs = dblarr(nwave,nnadabs)+1d
      modtmp = dblarr(nwave)+1d
      for i=0,nnadabs-1 do begin
         modnadabs[*,i]*=ifsf_cmpuvabs(wave,param[ilo+i*4:ilo+(i+1)*4-1])
         modtmp*=ifsf_cmpuvabs(wave,param[ilo+i*4:ilo+(i+1)*4-1])
      endfor
      modflux *= modtmp
   endif else modnadabs=0

;  Optionally, compute equivalent widths
   if keyword_set(weq) then begin
      dwave = wave[1:nwave-1] - wave[0:nwave-2]
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
      dwave = wave[1:nwave-1] - wave[0:nwave-2]
      if nnadem gt 0 then begin
         nademflux = dblarr(nnadem+1)
         for i=0,nnadem-1 do $
            nademflux[i+1] = $
               total((modnadem[1:nwave-1,i]*cont[1:nwave-1]-cont[1:nwave-1])*dwave)
         nademflux[0] = total(nademflux[1:nnadem])
      endif else nademflux=0d
   endif else nademflux=0d
   
   return,modflux
   
end
