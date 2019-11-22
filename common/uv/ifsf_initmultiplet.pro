; docformat = 'rst'
;
;+
;
; Initialize parameters for fitting multiplets.
;
; :Categories:
;    IFSFIT/INIT
;
; :Returns:
;    Parameter array (PARINFO) for input to MPFIT. Contents of PARINFO:
;      Elements 1-2: Number absorption components.
;      The first set of elements is
;        quadruplets of  abs. parameters: covering factor, optical depth, 
;        central wavelength (A), and sigma (km/s), all applying to the red
;        line.
;        
;
; :Params:
;    initabs: in, required, type=dblarr(M,4)
;      Initial values for  absorption line parameters (covering factor, 
;      optical depth of red line, wavelength in A of red line, sigma in km/s) 
;      for each of M components.
;    siglimabs: in, required, type=dblarr(2)
;      Limits to sigma for  absorption.
;      
; :Keywords:
;    taumax: in, optional, type=double, default=5d
;      Upper limit to optical depth.
;    absfix: in, required, type=bytarr(nabs,4)
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
;      2019nov20, DSNR, copied from IFSF_INITDOUBLET
;    
; :Copyright:
;    Copyright (C) 2014--2019 David S. N. Rupke
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
function ifsf_initmultiplet,initabs,siglimabs,reflinename,$
                            taumax=taumax,absfix=absfix
                          
   c = 299792.458d
   
   if not keyword_set(taumax) then taumax=5d

;  Get numbers of components
   size_abs = size(initabs)
   if size_abs[0] eq 1 then nabs = 1 $
   else if size_abs[0] gt 1 then nabs = fix(size_abs[1]) $
   else nabs=0l
  
;  Initialize PARINFO  
   parinfo = REPLICATE({value:0d, fixed:0b, limited:[0B,0B], tied:'', $
                       limits:[0d,0d], step:0d, mpprint:0b, mpside:2, $
                       parname:'', line:'', comp:0d}, $
                       1+nabs*4)

;  Record numbers of components
   parinfo[0].value = nabs
   parinfo[0].fixed = 1B
   parinfo[0].parname = 'No. of absorption components'

;  absorption
   if nabs gt 0 then begin
;     Initial values
      ilo = 1
      parinfo[ilo:ilo+nabs*4-1].value = $
         reform(transpose(initabs),nabs*4)
;     Index arrays
      ind_c = ilo + indgen(nabs)*4
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
      parinfo[ind_s].limits[0]  = siglimabs[0]
      parinfo[ind_s].limits[1]  = siglimabs[1]
      parinfo[ind_w].limited[0] = 1B
      parinfo[ind_w].limited[1] = 1B
      for i=0,nabs-1 do begin
         parinfo[ilo+2+i*4].limits[0] = initabs[i,2]-2d
         parinfo[ilo+2+i*4].limits[1] = initabs[i,2]+2d
      endfor
      if keyword_set(absfix) then $
         parinfo[ilo:ilo+nabs*4-1].fixed = $
            reform(transpose(absfix),nabs*4)
;     Labels
      parinfo[ind_c].parname = 'covering_factor'
      parinfo[ind_t].parname = 'optical_depth'
      parinfo[ind_w].parname = 'wavelength'
      parinfo[ind_s].parname = 'sigma'
      parinfo[ilo:ilo+nabs*4-1].line = reflinename
      parinfo[ilo:ilo+nabs*4-1].comp = rebin(indgen(nabs)+1,nabs*4)
   endif

;  Check parinit initial values vs. limits
   badpar = where((parinfo.limited[0] AND $
                   parinfo.value lt parinfo.limits[0]) OR $
                  (parinfo.limited[1] AND $
                   parinfo.value gt parinfo.limits[1]),ct)
   if ct gt 0 then begin
      print,'IFSF_INITMULTILPLET: Initial values are outside limits.'
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
