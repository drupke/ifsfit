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
;    A doublet absorption specturm.
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
                          moddoubletabs=moddoubletabs, $
                          moddoubletem=moddoubletem, $
                          doubletemflux=doubletemflux, $
                          vels=vels,vwtabs=vwtabs

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
         lratio = 2802.7056d/2795.5301d
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

;  Numbers of components
   ndoubletabs = param[0]
   ndoubletem = param[1]

   nwave = n_elements(wave)

   modflux = dblarr(nwave)+1d
;  doublet emission
   ilo = 2+ndoubletabs*4
   if ndoubletem gt 0 then moddoubletem=dblarr(nwave,ndoubletem)+1d $
   else moddoubletem=0
   for i=0,ndoubletem-1 do begin
      arg1 = ((wave-param[ilo+i*4])/$
             (2d*param[ilo+i*4]*param[1+ilo+i*4]/c))^2d
      arg2 = ((lratio * wave-param[ilo+i*4])/$
             (2d*param[ilo+i*4]*param[1+ilo+i*4]/c))^2d
      mask1 = (arg1 lt 80)
      mask2 = (arg2 lt 80)
      moddoubletem[*,i]+=param[2+ilo+i*4]*(mask1*exp(-arg1*mask1) + $
                                       param[3+ilo+i*4]*mask2*exp(-arg2*mask2))
      modflux+=param[2+ilo+i*4]*(mask1*exp(-arg1*mask1) + $
                                 param[3+ilo+i*4]*mask2*exp(-arg2*mask2))      
   endfor
;  doublet absorption
   ilo = 2
   if ndoubletabs gt 0 then begin
      moddoubletabs = dblarr(nwave,ndoubletabs)+1d
      modtmp = dblarr(nwave)+1d
      for i=0,ndoubletabs-1 do begin
         moddoubletabs[*,i]*=ifsf_cmpdoublet(wave,param[ilo+i*4:ilo+(i+1)*4-1],$
                                             tratio,lratio)
         modtmp*=ifsf_cmpdoublet(wave,param[ilo+i*4:ilo+(i+1)*4-1],$
                                 tratio,lratio)
      endfor
      modflux *= modtmp
   endif else moddoubletabs=0

;  Optionally, compute equivalent widths
   if keyword_set(weq) then begin
      dwave = wave[1:nwave-1] - wave[0:nwave-2]
      if ndoubletem gt 0 then begin
         emweq = dblarr(ndoubletem+1)
         for i=0,ndoubletem-1 do $
            emweq[i+1] = total((1d - moddoubletem[1:nwave-1,i])*dwave)
         emweq[0] = total(emweq[1:ndoubletem])
      endif else emweq=0d
      if ndoubletabs gt 0 then begin
         absweq = dblarr(ndoubletabs+1)
         for i=0,ndoubletabs-1 do $
            absweq[i+1] = total((1d - moddoubletabs[1:nwave-1,i])*dwave)
;        Total doublet equivalent width has to be calculated from full absorption
;        spectrum, rather than taken as the total of separate components.
         absweq[0] = total((1d - modtmp[1:nwave-1])*dwave)
      endif else absweq=0d
      weq = {em: emweq,abs: absweq}
   endif
   if keyword_set(doubletemflux) AND keyword_set(cont) then begin
      dwave = wave[1:nwave-1] - wave[0:nwave-2]
      if ndoubletem gt 0 then begin
         doubletemflux = dblarr(ndoubletem+1)
         for i=0,ndoubletem-1 do $
            doubletemflux[i+1] = $
               total((moddoubletem[1:nwave-1,i]*cont[1:nwave-1]-cont[1:nwave-1])*dwave)
         doubletemflux[0] = total(doubletemflux[1:ndoubletem])
      endif else doubletemflux=0d
   endif else doubletemflux=0d

;  Weighted average velocity and weighted RMS width
;  Trump et al. 2006, ApJS, 165, 1
;  Section 4.5
   if keyword_set(vwtabs) AND keyword_set(vels) AND ndoubletabs gt 0 then begin
;     Compute singlet version of profile
      ilo = 2
      modtmp = dblarr(nwave)+1d
      for i=0,ndoubletabs-1 do $
         modtmp*=ifsf_cmpsinglet(wave,param[ilo+i*4:ilo+(i+1)*4-1],$
                                 tratio,lratio)
;      igdvels = where(vels lt 0 AND vels gt -2.9d4,ctgdvels)
;      gdvels = abs(vels[igdvels])
      igdvels = where(vels gt -2.9d4,ctgdvels)
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

   return,modflux
   
end
