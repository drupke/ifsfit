; docformat = 'rst'
;
;+
;
; Initialize parameters for fitting. Specific to GMOS instrument.
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
;      2009jun01, DSNR, created
;      2009jun08, DSNR, added multiple components
;      2010may27, DSNR, re-written to fit in observed frame
;      2013sep13, DSNR, re-written to allow more than one common redshift
;      2013dec12, DSNR, documented, renamed, added license and copyright 
;      2014jan13, DSNR, updated to use hashes, and to add parname, line, and 
;                       comp tags into output parinfo structure
;      2014apr10, DSNR, added if statements to remove IEEE exceptions
;      2014apr17, DSNR, adjusted upper limits for Ha/Hb and [NII]/Ha
;      2014apr23, DSNR, added SIGFIX keyword
;      2014jun05, DSNR, added LRATFIX keyword
;    
; :Copyright:
;    Copyright (C) 2013-2014 David S. N. Rupke
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
function ifsf_gmos,linelist,linelistz,linetie,$
                   initflux,initsig,maxncomp,ncomp,$
                   lratfix=lratfix,siglim=siglim,sigfix=sigfix

  bad = 1d99
  c = 299792.458d

; Sigma limits
  gmosres = 3000d
  if ~ keyword_set(siglim) then siglim=[299792d/gmosres/2.35d,2000d]

; Number of emission lines to fit
  nline = linelist->count()
  lines_arr = (linelist->keys())->toarray()
  
; Number of initial parameters before Gaussian parameters begin
  lratlim = 4 ; maximum number of line ratios to constrain
  ppoff0 = 2
  ppoff = ppoff0 + maxncomp*lratlim
  
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

  if ~ keyword_set(lratfix) then lratfix = hash()
     
; [SII] ratio
  ilratlim = 0
  lratlab = '[SII]6716/6731'
  if ncomp->haskey('[SII]6716') then tmp_ncomp = ncomp['[SII]6716'] $
  else tmp_ncomp=0
  if tmp_ncomp gt 0 then begin
     ip1 = ppoff0 + ilratlim*maxncomp
     ip2 = ip1+tmp_ncomp-1
     fa = initflux['[SII]6716',0:tmp_ncomp-1]
     fb = initflux['[SII]6731',0:tmp_ncomp-1]
     frat = dblarr(tmp_ncomp)+1d ; default if initial s2b flux = 0
     inz = where(fb gt 0,ctnz)
     if ctnz gt 0 then frat[inz] = fa[inz]/fb[inz]
     parinfo[ip1:ip2].value = frat      
     parinfo[ip1:ip2].limited = rebin([1b,1b],2,tmp_ncomp)
     parinfo[ip1:ip2].limits  = rebin([0.44d,1.43d],2,tmp_ncomp)
     parinfo[ip1:ip2].parname = '[SII]6716/6731 line ratio'
     parinfo[ip1:ip2].comp = indgen(tmp_ncomp)+1
;    Check to see if line ratio is fixed
     ilratfix = where(lratfix.keys() eq lratlab,ctlratfix)
     for i=0,tmp_ncomp-1 do begin
;       If line ratio is fixed, then fix it
        lratfixed = 0b
        if ctlratfix gt 0 then begin
           if lratfix[lratlab,i] ne bad then begin
              parinfo[ip1+i].value = lratfix[lratlab,i]
              parinfo[ip1+i].fixed = 1b
              parinfo[ip1+i].limited = [0b,0b]
              lratfixed = 1b
           endif
        endif
        if ~ lratfixed then begin
;          case of pegging at or exceeding upper limit
           if parinfo[ip1+i].value ge parinfo[ip1+i].limits[1] then $
              parinfo[ip1+i].value = parinfo[ip1+i].limits[1] - $
                                     (parinfo[ip1+i].limits[1] - $
                                      parinfo[ip1+i].limits[0])*0.1
;          case of pegging at or dipping below lower limit
           if parinfo[ip1+i].value le parinfo[ip1+i].limits[0] then $
              parinfo[ip1+i].value = parinfo[ip1+i].limits[0] + $
                                     (parinfo[ip1+i].limits[1] - $
                                      parinfo[ip1+i].limits[0])*0.1
        endif
     endfor
  endif
  
; [NI] ratio
; See Ferland+12 for collisional case, Bautista99 for other cases. Upper limit 
; was originally 3, but found that it would peg at that and then the error for 
; [NI]5200 would be artificially large, and it would be removed from the fit. 
; Can fix to 1.5 (low density collisional limit, applicable to n <~ 10^3 
; cm^-3; Ferland+12 Appendix A.3) to solve artificially large errors in [NI]5200. 
  ilratlim = 1
  lratlab = '[NI]5200/5198'
  if ncomp->haskey('[NI]5198') then tmp_ncomp = ncomp['[NI]5198'] $
  else tmp_ncomp = 0
  if tmp_ncomp gt 0 then begin
     ip1 = ppoff0 + ilratlim*maxncomp
     ip2 = ip1+tmp_ncomp-1
     fa = initflux['[NI]5200',0:tmp_ncomp-1]
     fb = initflux['[NI]5198',0:tmp_ncomp-1]
     frat = dblarr(tmp_ncomp)+2d ; default if initial n1a flux = 0
     inz = where(fb gt 0,ctnz)
     if ctnz gt 0 then frat[inz] = fa[inz]/fb[inz]
     parinfo[ip1:ip2].value = frat
     parinfo[ip1:ip2].limited = rebin([1b,1b],2,tmp_ncomp)
     parinfo[ip1:ip2].limits  = rebin([0.6d,4d],2,tmp_ncomp)
     parinfo[ip1:ip2].parname = '[NI]5200/5198 line ratio'
     parinfo[ip1:ip2].comp = indgen(tmp_ncomp)+1
;    Check to see if line ratio is fixed
     ilratfix = where(lratfix.keys() eq lratlab,ctlratfix)
     for i=0,tmp_ncomp-1 do begin
;       If line ratio is fixed, then fix it
        lratfixed = 0b
        if ctlratfix gt 0 then begin
           if lratfix[lratlab,i] ne bad then begin
              parinfo[ip1+i].value = lratfix[lratlab,i]
              parinfo[ip1+i].fixed = 1b
              parinfo[ip1+i].limited = [0b,0b]
              lratfixed = 1b
           endif
        endif
        if ~ lratfixed then begin
;          case of pegging at or exceeding upper limit
           if parinfo[ip1+i].value ge parinfo[ip1+i].limits[1] then $
              parinfo[ip1+i].value = parinfo[ip1+i].limits[1] - $
                                     (parinfo[ip1+i].limits[1] - $
                                      parinfo[ip1+i].limits[0])*0.1
;          case of pegging at or dipping below lower limit
           if parinfo[ip1+i].value le parinfo[ip1+i].limits[0] then $
              parinfo[ip1+i].value = parinfo[ip1+i].limits[0] + $
                                     (parinfo[ip1+i].limits[1] - $
                                      parinfo[ip1+i].limits[0])*0.1
        endif
     endfor
  endif

; [NII]/Ha ratio
  ilratlim = 2
  lratlab = '[NII]6583/Ha'
  if ncomp->haskey('Halpha') then tmp_ncomp = ncomp['Halpha'] $
  else tmp_ncomp = 0
  if tmp_ncomp gt 0 then begin
     ip1 = ppoff0 + ilratlim*maxncomp
     ip2 = ip1 + tmp_ncomp - 1
     fa = initflux['[NII]6583',0:tmp_ncomp-1]
     fb = initflux['Halpha',0:tmp_ncomp-1]
     frat = dblarr(tmp_ncomp)+1d ; default if initial ha flux = 0
     inz = where(fb gt 0,ctnz)
     if ctnz gt 0 then frat[inz] = fa[inz]/fb[inz]
     parinfo[ip1:ip2].value = frat
     parinfo[ip1:ip2].limited = rebin([1B,1B],2,tmp_ncomp)
;    This upper limit appears to be the maximum seen in Kewley+06 or 
;    Rich+14 ("Composite Spectra in ..."). The lower limit is appropriate 
;    for ULIRGs.
     parinfo[ip1:ip2].limits  = rebin([0.1d,4d],2,tmp_ncomp)
     parinfo[ip1:ip2].parname = '[NII]/Halpha line ratio'
     parinfo[ip1:ip2].comp = indgen(tmp_ncomp)+1
;    Check to see if line ratio is fixed
     ilratfix = where(lratfix.keys() eq lratlab,ctlratfix)
     for i=0,tmp_ncomp-1 do begin
;       If line ratio is fixed, then fix it
        lratfixed = 0b
        if ctlratfix gt 0 then begin
           if lratfix[lratlab,i] ne bad then begin
              parinfo[ip1+i].value = lratfix[lratlab,i]
              parinfo[ip1+i].fixed = 1b
              parinfo[ip1+i].limited = [0b,0b]
              lratfixed = 1b
           endif
        endif
        if ~ lratfixed then begin
;          case of pegging at or exceeding upper limit
           if parinfo[ip1+i].value ge parinfo[ip1+i].limits[1] then $
              parinfo[ip1+i].value = parinfo[ip1+i].limits[1] - $
                                     (parinfo[ip1+i].limits[1] - $
                                      parinfo[ip1+i].limits[0])*0.1
;          case of pegging at or dipping below lower limit
           if parinfo[ip1+i].value le parinfo[ip1+i].limits[0] then $
              parinfo[ip1+i].value = parinfo[ip1+i].limits[0] + $
                                     (parinfo[ip1+i].limits[1] - $
                                      parinfo[ip1+i].limits[0])*0.1
        endif
     endfor
  endif

; Ha/Hb ratio
  ilratlim = 3
  lratlab = 'Ha/Hb'
  if ncomp->haskey('Halpha') then tmp_ncomp = ncomp['Halpha'] $
  else tmp_ncomp = 0
  if tmp_ncomp gt 0 then begin
     ip1 = ppoff0 + ilratlim*maxncomp
     ip2 = ip1 + tmp_ncomp - 1
     fa = initflux['Halpha',0:tmp_ncomp-1]
     fb = initflux['Hbeta',0:tmp_ncomp-1]
     frat = dblarr(tmp_ncomp)+3d ; default if initial hb flux = 0
     inz = where(fb gt 0,ctnz)
     if ctnz gt 0 then frat[inz] = fa[inz]/fb[inz]
     parinfo[ip1:ip2].value = frat
     parinfo[ip1:ip2].limited = rebin([1B,1B],2,tmp_ncomp)
;    Upper limit of 50 corresponds to E(B-V) = 2.89 using CCM
     parinfo[ip1:ip2].limits  = rebin([2.86d,50d],2,tmp_ncomp)
     parinfo[ip1:ip2].parname = 'Halpha/Hbeta line ratio'
     parinfo[ip1:ip2].comp = indgen(tmp_ncomp)+1
;    Check to see if line ratio is fixed
     ilratfix = where(lratfix.keys() eq lratlab,ctlratfix)
     for i=0,tmp_ncomp-1 do begin
;       If line ratio is fixed, then fix it
        lratfixed = 0b
        if ctlratfix gt 0 then begin
           if lratfix[lratlab,i] ne bad then begin
              parinfo[ip1+i].value = lratfix[lratlab,i]
              parinfo[ip1+i].fixed = 1b
              parinfo[ip1+i].limited = [0b,0b]
              lratfixed = 1b
           endif
        endif
        if ~ lratfixed then begin
;          case of pegging at or exceeding upper limit
           if parinfo[ip1+i].value ge parinfo[ip1+i].limits[1] then $
              parinfo[ip1+i].value = parinfo[ip1+i].limits[1] - $
                                     (parinfo[ip1+i].limits[1] - $
                                      parinfo[ip1+i].limits[0])*0.1
;          case of pegging at or dipping below lower limit
           if parinfo[ip1+i].value le parinfo[ip1+i].limits[0] then $
              parinfo[ip1+i].value = parinfo[ip1+i].limits[0] + $
                                     (parinfo[ip1+i].limits[1] - $
                                      parinfo[ip1+i].limits[0])*0.1
        endif
     endfor
  endif

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

;          initial values
           parinfo[ifoff].value = initflux[line,i]
           parinfo[iwoff].value = linelistz[line,i]
           parinfo[isoff].value = initsig[line,i]
;          limits
           parinfo[ifoff].limited[0] = 1B
           parinfo[ifoff].limits[0]  = 0d
           parinfo[iwoff].limited = [1B,1B]
           parinfo[iwoff].limits[0] = linelistz[line,i]*0.997d
           parinfo[iwoff].limits[1] = linelistz[line,i]*1.003d
           parinfo[isoff].limited = [1B,1B]
           parinfo[isoff].limits = siglim
;          ties
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
;          fixed/free
           if keyword_set(sigfix) then $
              if sigfix.haskey(line) then $
                 if sigfix[line,i] ne 0 then begin
                    parinfo[isoff].fixed=1B
                    parinfo[isoff].value=sigfix[line,i]
                 endif
                 
        endelse

        iline++

     endforeach

;    the if statement here prevents MPFIT_TIE from issuing an IEEE exception,
;    since if we're not fitting [SII] then the ratio is set to 0 at the 
;    beginning of this routine
     if ncomp->haskey('[SII]6716') then begin
        if ncomp['[SII]6716'] gt 0 then begin
           ilratlim = 0
           linea = where(lines_arr eq '[SII]6716')
           lineb  = where(lines_arr eq '[SII]6731')
           parinfo[foff+lineb*3].tied = 'P['+$
                                        string(foff+linea*3,$
                                               format='(I0)')+']/P['+$
                                        string(ppoff0+maxncomp*ilratlim+i,$
                                               format='(I0)')+']'
        endif
     endif                                              

     ilratlim = 1
     linea = where(lines_arr eq '[NI]5198')
     lineb = where(lines_arr eq '[NI]5200')
     parinfo[foff+linea*3].tied = 'P['+$
                                  string(foff+lineb*3,$
                                         format='(I0)')+']/P['+$
                                  string(ppoff0+maxncomp*ilratlim+i,$
                                         format='(I0)')+']'
                                         
     ilratlim = 2
     linea = where(lines_arr eq 'Halpha')
     lineb = where(lines_arr eq '[NII]6583')
     parinfo[foff+lineb*3].tied = 'P['+$
                                  string(foff+linea*3,$
                                         format='(I0)')+']*P['+$
                                  string(ppoff0+maxncomp*ilratlim+i,$
                                         format='(I0)')+']'

;    the if statement here prevents MPFIT_TIE from issuing an IEEE exception,
;    since if we're not fitting Halpha then the ratio is set to 0 at the
;    beginning of this routine
     if ncomp->haskey('Halpha') then begin
        if ncomp['Halpha'] gt 0 then begin
           ilratlim = 3
           linea = where(lines_arr eq 'Halpha')
           lineb = where(lines_arr eq 'Hbeta')
           parinfo[foff+lineb*3].tied = 'P['+$
                                        string(foff+linea*3,$
                                               format='(I0)')+']/P['+$
                                        string(ppoff0+maxncomp*ilratlim+i,$
                                               format='(I0)')+']'
        endif
     endif

     linea = where(lines_arr eq '[OIII]4959')
     lineb = where(lines_arr eq '[OIII]5007')
     parinfo[foff+linea*3].tied = 'P['+$
                                  string(foff+lineb*3,$
                                         format='(I0)')+']/3.0d'
                                         
     linea = where(lines_arr eq '[OI]6300')
     lineb = where(lines_arr eq '[OI]6364')
     parinfo[foff+lineb*3].tied = 'P['+$
                                  string(foff+linea*3,$
                                         format='(I0)')+']/3.0d'
     linea = where(lines_arr eq '[NII]6548')
     lineb = where(lines_arr eq '[NII]6583')
     parinfo[foff+linea*3].tied = 'P['+ $
                                  string(foff+lineb*3,$
                                         format='(I0)')+']/3.0d'

  endfor

; Check parinit initial values vs. limits
  badpar = where((parinfo.limited[0] AND $
                  parinfo.value lt parinfo.limits[0]) OR $
                 (parinfo.limited[1] AND $
                  parinfo.value gt parinfo.limits[1]),ct)
  if ct gt 0 then begin
     print,'IFSF_GMOS: Initial values are outside limits.'
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
