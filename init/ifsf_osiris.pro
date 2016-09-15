; docformat = 'rst'
;
;+
;
; Initialize parameters for fitting.
;
; :Categories:
;    IFSFIT/INIT
;
; :Returns:
;    PARINFO structure for input into MPFIT.
;
; :Params:
;    linelist: in, required, type=hash(lines)
;      Emission line rest frame wavelengths.
;    linelistz: in, required, type=hash(lines\,maxncomp)
;      Emission line observed frame wavelengths.
;    linetie: in, required, type=hash(lines)
;      Name of emission line to which each emission line is tied
;      (in redshift and linewidth).
;    initflux: in, required, type=hash(lines\,maxncomp)
;      Initial guess for peak flux in each component.
;    initsig: in, required, type=hash(lines\,maxncomp)
;      Initial guess for emission lines widths.
;    maxncomp: in, required, type=double
;      Maximum no. of emission line components.
;    ncomp: in, required, type=hash(lines)
;      Number of velocity components.
;      
; :Keywords:
;    lratfix: in, optional, type=hash(lineratios,ncomp)
;      For each line ratio that should be fixed, input an array with each 
;      element set to either the BAD value (do not fix that component) or the 
;      value to which the line ratio will be fixed for that component.
;    siglim: in, optional, type=dblarr(2)
;      Lower and upper sigma limits in km/s.
;    sigfix: in, optional, type=hash(lines\,maxncomp)
;      Fix sigma at this value, for particular lines/components.
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
;      2016sep13, DSNR, copied from IFSF_GMOS
;    
; :Copyright:
;    Copyright (C) 2016 David S. N. Rupke
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
function ifsf_osiris,linelist,linelistz,linetie,$
                     initflux,initsig,maxncomp,ncomp,$
                     lratfix=lratfix,siglim=siglim,sigfix=sigfix

  bad = 1d99
  c = 299792.458d

; Sigma limits
  res = 3000d
  if ~ keyword_set(siglim) then siglim=[299792d/res/2.35d,2000d]

; Number of emission lines to fit
  nline = linelist->count()
  lines_arr = (linelist->keys())->toarray()
  
; Number of initial parameters before Gaussian parameters begin
  ppoff0 = 2
  ppoff = ppoff0
  
  parinfo = REPLICATE({value:0d, fixed:0b, limited:[0B,0B], tied:'', $
                       limits:[0d,0d], step:0d, mpprint:0b, mpside:2, $
                       parname:'', line:'', comp:0d}, $
                      ppoff+maxncomp*(nline*3))

; Number of initial parameters before Gaussian parameters begin
  parinfo[0].value = ppoff
  parinfo[0].fixed = 1B
  parinfo[0].parname = 'No. of non-Gaussian parameters'

; Maximum number of velocity components
  parinfo[1].value = maxncomp
  parinfo[1].fixed = 1B
  parinfo[1].parname = 'Maximum no. of velocity components'

; cycle through velocity components
  for i=0,maxncomp-1 do begin

;    index offsets for this component
     foff = ppoff+i*nline*3
     woff = foff+1
     soff = foff+2

;    cycle through lines
     iline=0
     foreach line,lines_arr do begin

;       indices
        ifoff = foff + iline*3
        iwoff = woff + iline*3
        isoff = soff + iline*3

        parinfo[ifoff].parname = 'flux_peak'
        parinfo[iwoff].parname = 'wavelength'
        parinfo[isoff].parname = 'sigma'
        parinfo[ifoff].line = line
        parinfo[iwoff].line = line
        parinfo[isoff].line = line
        parinfo[ifoff].comp = i+1
        parinfo[iwoff].comp = i+1
        parinfo[isoff].comp = i+1
        
;       if the number of components to be fit is exceeded, fix line fluxes to 0
        if i+1 gt ncomp[line] then begin 

           parinfo[ifoff].value = 0d
           parinfo[iwoff].value = 0d
           parinfo[isoff].value = 0d
           parinfo[ifoff].fixed = 1b
           parinfo[iwoff].fixed = 1b
           parinfo[isoff].fixed = 1b

        endif else begin

;             initial values
              parinfo[ifoff].value = initflux[line,i]
              parinfo[iwoff].value = linelistz[line,i]
              parinfo[isoff].value = initsig[line,i]
;             limits
              parinfo[ifoff].limited[0] = 1B
              parinfo[ifoff].limits[0]  = 0d
              parinfo[iwoff].limited = [1B,1B]
              parinfo[iwoff].limits[0] = linelistz[line,i]*0.997d
              parinfo[iwoff].limits[1] = linelistz[line,i]*1.003d
              parinfo[isoff].limited = [1B,1B]
              parinfo[isoff].limits = siglim
;             ties
              if (line eq linetie[line]) then begin
                 parinfo[iwoff].tied = ''
                 parinfo[isoff].tied = ''
              endif else begin
                 indtie = where(lines_arr eq linetie[line])
                 parinfo[iwoff].tied = $
                    string(linelist[line],'/',linelist[linetie[line]],$
                           '* P[',woff+indtie*3,']',$
                           format='(D0.2,A0,D0.2,A0,I0,A0)') 
                 parinfo[isoff].tied = $
                    string('P[',soff+indtie*3,']',format='(A0,I0,A0)')   
              endelse
;             fixed/free
              if keyword_set(sigfix) then $
                 if sigfix.haskey(line) then $
                    if sigfix[line,i] ne 0 then begin
                       parinfo[isoff].fixed=1B
                       parinfo[isoff].value=sigfix[line,i]
                    endif

                 
        endelse

        iline++

     endforeach
   
  endfor

; Check parinit initial values vs. limits
  badpar = where((parinfo.limited[0] AND $
                  parinfo.value lt parinfo.limits[0]) OR $
                 (parinfo.limited[1] AND $
                  parinfo.value gt parinfo.limits[1]),ct)
  if ct gt 0 then begin
     print,'IFSF_OSIRIS: Initial values are outside limits.'
     print,'Offending parameters:'
     print,'Quantity','Line','Comp','Value','Lower limit','Upper limit',$
           format='(2A10,A5,3A15)'
     for i=0,ct-1 do begin
        j = badpar[i]
        print,parinfo[j].parname,parinfo[j].line,parinfo[j].comp,$
              parinfo[j].value,parinfo[j].limits[0],$
              parinfo[j].limits[1],format='(2A10,I5,3E15.6)'
     endfor
     return,0
  endif else begin
     return,parinfo
  endelse

end
