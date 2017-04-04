; docformat = 'rst'
;
;+
;
; Initialize parameters for fitting NaD + HeI.
;
; :Categories:
;    IFSFIT/INIT
;
; :Returns:
;    Parameter array (PARINFO) for input to MPFIT. Contents of PARINFO:
;      Elements 1-3: Number of HeI emission, NaD absorption, and NaD emission 
;        components.
;      If HeI present, first set of elements is triplets of HeI parameters:
;        central wavelength (in A), sigma (in km/s), and peak flux.
;      The next set of elements (or the first set if HeI is not present) is
;        quadruplets of NaD abs. parameters: covering factor, optical depth, 
;        central wavelength (A), and sigma (km/s), all applying to the D1 (red)
;        line.
;      The next set of elements (or the first set if HeI and NaD absorption are 
;        both absent) is quadruplets of NaD emission parameters: 
;        central wavelength (A), and sigma (km/s), flux (all applying to the D1
;        line), and the D2/D1 flux ratio
;        
;
; :Params:
;    inithei: in, required, type=dblarr(N,3)
;      Initial values for HeI 5876 line parameters (wavelength in A, sigma in 
;      km/s, and peak flux) for each of N components.
;    initnadabs: in, required, type=dblarr(M,4)
;      Initial values for NaD absorption line parameters (covering factor, 
;      optical depth of D2, wavelength in A of D2, sigma in km/s) for each of 
;      M components.
;    initnadem: in, required, type=dblarr(L,4)
;      Initial values for NaD emission line parameters (wavelength in A of D2, 
;      sigma in km/s, peak flux of D2, and flux ratio D2/D1) for each of L 
;      components.
;    siglimnadabs: in, required, type=dblarr(2)
;      Limits to sigma for NaD absorption.
;    siglimnadem: in, required, type=dblarr(2)
;      Limits to sigma for NaD emission.
;      
; :Keywords:
;    taumax: in, optional, type=double, default=5d
;      Upper limit to optical depth.
;    heifix: in, optional, type=bytarr(nhei,3)
;      Input this array with parameters to be fixed set to 1.
;    siglimhei: in, required, type=dblarr(2)
;      Limits to sigma for HeI emission.
;    nadabsfix: in, required, type=bytarr(nnadabs,4)
;      Input this array with parameters to be fixed set to 1.
;    nademfix: in, required, type=bytarr(nnadem,4)
;      Input this array with parameters to be fixed set to 1.
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
;      2014may28, DSNR, added NADEMFIX parameter
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
function ifsf_initnad,inithei,initnadabs,initnadem,siglimnadabs,siglimnadem,$
                      taumax=taumax,heifix=heifix,siglimhei=siglimhei,$
                      nadabsfix=nadabsfix,nademfix=nademfix
                      

   c = 299792.458d
;  NaD optical depth ratio (blue to red)
   tratio = 2.0093d

   if not keyword_set(taumax) then taumax=5d

;  Get numbers of components
   size_hei = size(inithei)
   size_nadabs = size(initnadabs)
   size_nadem = size(initnadem)
   if size_hei[0] eq 1 then nhei = 1 $
   else if size_hei[0] gt 1 then nhei = fix(size_hei[1]) else nhei=0l
   if size_nadabs[0] eq 1 then nnadabs = 1 $
   else if size_nadabs[0] gt 1 then nnadabs = fix(size_nadabs[1]) else nnadabs=0l
   if size_nadem[0] eq 1 then nnadem = 1 $
   else if size_nadem[0] gt 1 then nnadem = fix(size_nadem[1]) else nnadem=0l
  
;  Initialize PARINFO  
   parinfo = REPLICATE({value:0d, fixed:0b, limited:[0B,0B], tied:'', $
                       limits:[0d,0d], step:0d, mpprint:0b, mpside:2, $
                       parname:'', line:'', comp:0d}, $
                       3+nhei*3+nnadabs*4+nnadem*4)

;  Record numbers of components
   parinfo[0].value = nhei
   parinfo[0].fixed = 1B
   parinfo[0].parname = 'No. of HeI components'
   parinfo[1].value = nnadabs
   parinfo[1].fixed = 1B
   parinfo[1].parname = 'No. of NaD absorption components'
   parinfo[2].value = nnadem
   parinfo[2].fixed = 1B
   parinfo[2].parname = 'No. of NaD emission components'

;  HeI
   if nhei gt 0 then begin
;     Initial values
      ilo = 3
      parinfo[ilo:ilo+nhei*3-1].value = reform(transpose(inithei),nhei*3)
;     Index arrays
      ind_w = 3 + indgen(nhei)*3
      ind_s = ind_w + 1
      ind_f = ind_w + 2
;     Limits + fixed/free
      parinfo[ind_f].limited[0] = 1B
      parinfo[ind_f].limits[0]  = 0d
      parinfo[ind_w].limited[0] = 1B
      parinfo[ind_w].limited[1] = 1B
      parinfo[ind_w].limits[0] = inithei[*,0]-10d
      parinfo[ind_w].limits[1]  = inithei[*,0]+10d
      parinfo[ind_s].limited[0] = 1B
      parinfo[ind_s].limited[1] = 1B
      if ~ keyword_set(siglimhei) then siglimhei=[5d,1000d]
      parinfo[ind_s].limits[0] = siglimhei[0]
      parinfo[ind_s].limits[1]  = siglimhei[1]
      if keyword_set(heifix) then $
         parinfo[ilo:ilo+nhei*3-1].fixed = reform(transpose(heifix),nhei*3)
;     Labels
      parinfo[ind_f].parname = 'flux_peak'
      parinfo[ind_w].parname = 'wavelength'
      parinfo[ind_s].parname = 'sigma'
      parinfo[ilo:ilo+nhei*3-1].line = 'HeI5876'
      parinfo[ilo:ilo+nhei*3-1].comp = rebin(indgen(nhei)+1,nhei*3)
   endif

;  NaD absorption
   if nnadabs gt 0 then begin
;     Initial values
      ilo = 3+nhei*3
      parinfo[ilo:ilo+nnadabs*4-1].value = reform(transpose(initnadabs),nnadabs*4)
;     Index arrays
      ind_c = ilo + indgen(nnadabs)*4
      ind_t = ind_c + 1
      ind_w = ind_c + 2
      ind_s = ind_c + 3
;     Limits + fixed/free
      parinfo[ind_c].limited[0] = 1B
      parinfo[ind_c].limited[1] = 1B
      parinfo[ind_c].limits[0]  = 0d
      parinfo[ind_c].limits[1]  = 1d
      parinfo[ind_t].limited[0] = 1B
      parinfo[ind_t].limited[1] = 1B
      parinfo[ind_t].limits[0]  = 0d
      parinfo[ind_t].limits[1]  = taumax
      parinfo[ind_s].limited[0] = 1B
      parinfo[ind_s].limited[1] = 1B
      parinfo[ind_s].limits[0]  = siglimnadabs[0]
      parinfo[ind_s].limits[1]  = siglimnadabs[1]
      parinfo[ind_w].limited[0] = 1B
      parinfo[ind_w].limited[1] = 1B
      for i=0,nnadabs-1 do begin
         parinfo[ilo+2+i*4].limits[0] = initnadabs[i,2]-10d
         parinfo[ilo+2+i*4].limits[1] = initnadabs[i,2]+10d
      endfor
      if keyword_set(nadabsfix) then $
         parinfo[ilo:ilo+nnadabs*4-1].fixed = $
            reform(transpose(nadabsfix),nnadabs*4)
;     Labels
      parinfo[ind_c].parname = 'covering_factor'
      parinfo[ind_t].parname = 'optical_depth'
      parinfo[ind_w].parname = 'wavelength'
      parinfo[ind_s].parname = 'sigma'
      parinfo[ilo:ilo+nnadabs*4-1].line = 'NaD1'
      parinfo[ilo:ilo+nnadabs*4-1].comp = rebin(indgen(nnadabs)+1,nnadabs*4)
   endif

;  NaD emission
   if nnadem gt 0 then begin
;     Initial values
      ilo = 3+nhei*3+nnadabs*4
      parinfo[ilo:ilo+nnadem*4-1].value = reform(transpose(initnadem),nnadem*4)
;     Index arrays
      ind_w = ilo + indgen(nnadem)*4
      ind_s = ind_w + 1
      ind_f = ind_w + 2
      ind_r = ind_w + 3
;     Limits + fixed/free
      parinfo[ind_w].limited[0] = 1B
      parinfo[ind_w].limited[1] = 1B
      for i=0,nnadem-1 do begin
         parinfo[ilo+i*4].limits[0] = initnadem[i,0]-10d
         parinfo[ilo+i*4].limits[1] = initnadem[i,0]+10d
      endfor
      parinfo[ind_s].limited[0] = 1B
      parinfo[ind_s].limited[1] = 1B
      parinfo[ind_s].limits[0]  = siglimnadem[0]
      parinfo[ind_s].limits[1]  = siglimnadem[1]
      parinfo[ind_f].limited[0] = 1B
      parinfo[ind_f].limits[0]  = 0d
      parinfo[ind_r].limited[0] = 1B
      parinfo[ind_r].limited[1] = 1B
      parinfo[ind_r].limits[0]  = 1d
      parinfo[ind_r].limits[1]  = tratio
      if keyword_set(nademfix) then $
         parinfo[ilo:ilo+nnadem*4-1].fixed = reform(transpose(nademfix),nnadem*4)
;     Labels
      parinfo[ind_w].parname = 'wavelength'
      parinfo[ind_s].parname = 'sigma'
      parinfo[ind_f].parname = 'flux_peak'
      parinfo[ind_r].parname = 'flux_ratio_D2D1'
      parinfo[ilo:ilo+nnadem*4-1].line = 'NaD1'
      parinfo[ilo:ilo+nnadem*4-1].comp = rebin(indgen(nnadem)+1,nnadem*4)
   endif

;  Check parinit initial values vs. limits
   badpar = where((parinfo.limited[0] AND $
                   parinfo.value lt parinfo.limits[0]) OR $
                  (parinfo.limited[1] AND $
                   parinfo.value gt parinfo.limits[1]),ct)
   if ct gt 0 then begin
      print,'IFSF_INITNAD: Initial values are outside limits.'
      print,'Offending parameters:'
      print,'Quantity','Line','Comp','Value','Lower limit','Upper limit',$
            format='(2A15,A5,3A15)'
      for i=0,ct-1 do begin
         j = badpar[i]
         print,parinfo[j].parname,parinfo[j].line,parinfo[j].comp,$
               parinfo[j].value,parinfo[j].limits[0],$
               parinfo[j].limits[1],format='(2A15,I5,3E15.6)'
      endfor
      return,0
   endif else begin
      return,parinfo
   endelse

end
