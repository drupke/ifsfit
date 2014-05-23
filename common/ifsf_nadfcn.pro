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
function ifsf_nadfcn, wave, param, modhei=modhei, modnadabs=modnadabs, $
                      modnadem=modnadem

   c = 299792.458d
;  NaD wavelength ratio (red to blue)
   lratio = 1.001014158d

;  Numbers of components
   nhei = param[0]
   nnadabs = param[1]
   nnadem = param[2]

   nwave = n_elements(wave)

;  Emission
   modflux = dblarr(nwave)+1d
   for i=0,nhei-1 do begin
      arg = ((wave-param[3+i*3])/$
            (2d*param[3+i*3]*param[4+i*3]/c))^2d
      mask = (arg lt 80)
      modflux += param[5+i*3]*mask*exp(-arg*mask)
   endfor
   if keyword_set(modhei) then begin
      if nhei gt 0 then begin
         modhei=dblarr(nwave,nhei)+1d
         for i=0,nhei-1 do begin
            arg = ((wave-param[3+i*3])/$
                  (2d*param[3+i*3]*param[4+i*3]/c))^2d
            mask = (arg lt 80)
            modhei[*,i] += param[5+i*3]*mask*exp(-arg*mask)
         endfor
      endif else modhei=0
   endif
   ilo = 3+nhei*3+nnadabs*4
   for i=0,nnadem-1 do begin
      arg1 = ((wave-param[ilo+i*4])/$
             (2d*param[ilo+i*4]*param[1+ilo+i*4]/c))^2d
      arg2 = ((wave-lratio*param[ilo+i*4])/$
             (2d*lratio*param[ilo+i*4]*param[1+ilo+i*4]/c))^2d
      mask1 = (arg1 lt 80)
      mask2 = (arg2 lt 80)
      modflux += param[2+ilo+i*4]*(mask1*exp(-arg1*mask1) + $
                                   param[3+ilo+i*4]*mask2*exp(-arg2*mask2))
   endfor
   if keyword_set(modnadem) then begin
      if nnadem gt 0 then begin
         modnadem=dblarr(nwave,nnadem)+1d
         for i=0,nnadem-1 do begin
            arg1 = ((wave-param[ilo+i*4])/$
                   (2d*param[ilo+i*4]*param[1+ilo+i*4]/c))^2d
            arg2 = ((wave-lratio*param[ilo+i*4])/$
                   (2d*lratio*param[ilo+i*4]*param[1+ilo+i*4]/c))^2d
            mask1 = (arg1 lt 80)
            mask2 = (arg2 lt 80)
            modnadem[*,i] += $
               param[2+ilo+i*4]*(mask1*exp(-arg1*mask1) + $
                                 param[3+ilo+i*4]*mask2*exp(-arg2*mask2))
         endfor
      endif else modnadem=0
   endif
;  Absorption
   modabs = dblarr(nwave)+1d
   ilo = 3+nhei*3
   for i=0,nnadabs-1 do modabs *= ifsf_cmpnad(wave,param[ilo+i*4:ilo+(i+1)*4-1])
   modflux *= modabs
   if keyword_set(modnadabs) then begin
      if nnadabs gt 0 then begin
         modnadabs=dblarr(nwave,nnadabs)+1d
         for i=0,nnadabs-1 do $
            modnadabs[*,i] *= ifsf_cmpnad(wave,param[ilo+i*4:ilo+(i+1)*4-1])
      endif else modnadabs=0
   endif
   
   return,modflux

end
