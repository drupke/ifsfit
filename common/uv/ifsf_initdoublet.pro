; docformat = 'rst'
;
;+
;
; Initialize parameters for fitting doublets.
;
; :Categories:
;    IFSFIT/INIT
;
; :Returns:
;    Parameter array (PARINFO) for input to MPFIT. Contents of PARINFO:
;      Elements 1-2: Number doublet absorption and emission 
;        components.
;      The first set of elements is
;        quadruplets of doublet abs. parameters: covering factor, optical depth, 
;        central wavelength (A), and sigma (km/s), all applying to the red
;        line.
;      The next set of elements is quadruplets of doublet emission parameters: 
;        central wavelength (A), and sigma (km/s), flux (all applying to the red
;        line), and the blue/red flux ratio
;        
;
; :Params:
;    doublet: in, required, type=string
;      Doublet to fit.
;    initdoubletabs: in, required, type=dblarr(M,4)
;      Initial values for doublet absorption line parameters (covering factor, 
;      optical depth of red line, wavelength in A of red line, sigma in km/s) 
;      for each of M components.
;    initdoubletem: in, required, type=dblarr(L,4)
;      Initial values for doublet emission line parameters (wavelength in A of
;      red line, sigma in km/s, peak flux of red line, and flux ratio of red 
;      to blue line) for each of L components.
;    siglimdoubletabs: in, required, type=dblarr(2)
;      Limits to sigma for doublet absorption.
;    siglimdoubletem: in, required, type=dblarr(2)
;      Limits to sigma for doublet emission.
;      
; :Keywords:
;    taumax: in, optional, type=double, default=5d
;      Upper limit to optical depth.
;    doubletabsfix: in, required, type=bytarr(ndoubletabs,4)
;      Input this array with parameters to be fixed set to 1.
;    doubletemfix: in, required, type=bytarr(ndoubletem,4)
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
function ifsf_initdoublet,doublet,initdoubletabs,initdoubletem,siglimdoubletabs,$
                          siglimdoubletem,taumax=taumax,$
                          doubletabsfix=doubletabsfix,doubletemfix=doubletemfix
                      
   if doublet eq 'OVI' then begin
      LogLF = 2.136d -1.834d
      linename = 'OVI1037'
   endif else if doublet eq 'NV' then begin
      LogLF = 2.286d -1.985d
      linename = 'NV1242'
   endif else if doublet eq 'PV' then begin
      LogLF = 2.722d -2.420d
      linename = 'PV1128'
   endif else if doublet eq 'MgII' then begin
      LogLF = 3.236d - 2.933d
      linename = 'MgII2802'
   endif else if doublet eq 'FeIIUV1' then begin
         LogLF = 2.793d - 2.252d
         linename = 'FeII2585'
   endif else if doublet eq 'FeIIUV2' then begin
         LogLF = 2.882d - 1.871d
         linename = 'FeII2373'
   endif else begin
      print,'IFSF_INITDOUBLET: Error -- doublet not specified.'
      stop
   endelse
   
   tratio = 10d^LogLF
   c = 299792.458d
   
   if not keyword_set(taumax) then taumax=5d

;  Get numbers of components
   size_doubletabs = size(initdoubletabs)
   size_doubletem = size(initdoubletem)
   if size_doubletabs[0] eq 1 then ndoubletabs = 1 $
   else if size_doubletabs[0] gt 1 then ndoubletabs = fix(size_doubletabs[1]) $
   else ndoubletabs=0l
   if size_doubletem[0] eq 1 then ndoubletem = 1 $
   else if size_doubletem[0] gt 1 then ndoubletem = fix(size_doubletem[1]) $
   else ndoubletem=0l
  
;  Initialize PARINFO  
   parinfo = REPLICATE({value:0d, fixed:0b, limited:[0B,0B], tied:'', $
                       limits:[0d,0d], step:0d, mpprint:0b, mpside:2, $
                       parname:'', line:'', comp:0d}, $
                       2+ndoubletabs*4+ndoubletem*4)

;  Record numbers of components
   parinfo[0].value = ndoubletabs
   parinfo[0].fixed = 1B
   parinfo[0].parname = 'No. of Doublet absorption components'
   parinfo[1].value = ndoubletem
   parinfo[1].fixed = 1B
   parinfo[1].parname = 'No. of Doublet emission components'

;  doublet absorption
   if ndoubletabs gt 0 then begin
;     Initial values
      ilo = 2
      parinfo[ilo:ilo+ndoubletabs*4-1].value = $
         reform(transpose(initdoubletabs),ndoubletabs*4)
;     Index arrays
      ind_c = ilo + indgen(ndoubletabs)*4
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
      parinfo[ind_s].limits[0]  = siglimdoubletabs[0]
      parinfo[ind_s].limits[1]  = siglimdoubletabs[1]
      parinfo[ind_w].limited[0] = 1B
      parinfo[ind_w].limited[1] = 1B
      for i=0,ndoubletabs-1 do begin
         parinfo[ilo+2+i*4].limits[0] = initdoubletabs[i,2]-2d
         parinfo[ilo+2+i*4].limits[1] = initdoubletabs[i,2]+2d
      endfor
      if keyword_set(doubletabsfix) then $
         parinfo[ilo:ilo+ndoubletabs*4-1].fixed = $
            reform(transpose(doubletabsfix),ndoubletabs*4)
;     Labels
      parinfo[ind_c].parname = 'covering_factor'
      parinfo[ind_t].parname = 'optical_depth'
      parinfo[ind_w].parname = 'wavelength'
      parinfo[ind_s].parname = 'sigma'
      parinfo[ilo:ilo+ndoubletabs*4-1].line = linename
      parinfo[ilo:ilo+ndoubletabs*4-1].comp = rebin(indgen(ndoubletabs)+1,ndoubletabs*4)
   endif

;  doublet emission
   if ndoubletem gt 0 then begin
;     Initial values
;      ilo = 3+nhei*3+ndoubletabs*4
      ilo = 2+ndoubletabs*4
      parinfo[ilo:ilo+ndoubletem*4-1].value = reform(transpose(initdoubletem),ndoubletem*4)
;     Index arrays
      ind_w = ilo + indgen(ndoubletem)*4
      ind_s = ind_w + 1
      ind_f = ind_w + 2
      ind_r = ind_w + 3
;     Limits + fixed/free
      parinfo[ind_w].limited[0] = 1B
      parinfo[ind_w].limited[1] = 1B
      for i=0,ndoubletem-1 do begin
         parinfo[ilo+i*4].limits[0] = initdoubletem[i,0]-10d
         parinfo[ilo+i*4].limits[1] = initdoubletem[i,0]+10d
      endfor
      parinfo[ind_s].limited[0] = 1B
      parinfo[ind_s].limited[1] = 1B
      parinfo[ind_s].limits[0]  = siglimdoubletem[0]
      parinfo[ind_s].limits[1]  = siglimdoubletem[1]
      parinfo[ind_f].limited[0] = 1B
      parinfo[ind_f].limits[0]  = 0d
      parinfo[ind_r].limited[0] = 1B
      parinfo[ind_r].limited[1] = 1B
      parinfo[ind_r].limits[0]  = 1d
      parinfo[ind_r].limits[1]  = tratio
      if keyword_set(doubletemfix) then $
         parinfo[ilo:ilo+ndoubletem*4-1].fixed = $
            reform(transpose(doubletemfix),ndoubletem*4)
;     Labels
      parinfo[ind_w].parname = 'wavelength'
      parinfo[ind_s].parname = 'sigma'
      parinfo[ind_f].parname = 'flux_peak'
      parinfo[ind_r].parname = 'flux_ratio_redblue'
      parinfo[ilo:ilo+ndoubletem*4-1].line = linename
      parinfo[ilo:ilo+ndoubletem*4-1].comp = $
         rebin(indgen(ndoubletem)+1,ndoubletem*4)
   endif

;  Check parinit initial values vs. limits
   badpar = where((parinfo.limited[0] AND $
                   parinfo.value lt parinfo.limits[0]) OR $
                  (parinfo.limited[1] AND $
                   parinfo.value gt parinfo.limits[1]),ct)
   if ct gt 0 then begin
      print,'IFSF_INITDOUBLET: Initial values are outside limits.'
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
