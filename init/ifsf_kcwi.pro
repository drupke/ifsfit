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
;    blrcomp: in, optional, type=dblarr(N_BLR)
;      For each velocity component to model as a broad line region
;      (BLR), put the index (unity-offset) of that component into this
;      scalar (or array if more than one component) and all fluxes but
;      Balmer line fluxes will be zeroed.
;    blrlines: in, optional, type=strarr(N_lines)
;      List of lines to fit with BLR component.
;    lratfix: in, optional, type=hash(lineratios,ncomp)
;      For each line ratio that should be fixed, input an array with each 
;      element set to either the BAD value (do not fix that component) or the 
;      value to which the line ratio will be fixed for that component.
;    siglim: in, optional, type=dblarr(2)
;      Lower and upper sigma limits in km/s.
;    sigfix: in, optional, type=hash(lines\,maxncomp)
;      Fix sigma at this value, for particular lines/components.
;    specres: in, optional, type=double, def=0.64d
;      Estimated spectral resolution in wavelength units (sigma).
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
;      2018aug12, DSNR, created
;    
; :Copyright:
;    Copyright (C) 2018 David S. N. Rupke
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
function ifsf_kcwi,linelist,linelistz,linetie,$
                   initflux,initsig,maxncomp,ncomp,$
                   lratfix=lratfix,siglim=siglim,sigfix=sigfix,$
                   blrcomp=blrcomp,blrlines=blrlines,specres=specres

  bad = 1d99
  c = 299792.458d
  if ~ keyword_set(blrlines) then $
     blrlines = ['Halpha','Hbeta','Hgamma','Hdelta','Hepsilon',$
                 'H8','H9','H10','H11']

; Estimated spectral resolution for BL grating based on measurements.
; Website says R = 1800 at 4550 A for 0.7" slit. This gives 2.53 A FWHM. 
; Sigma is then 1.08. 
  if ~ keyword_set(specres) then specres = 1.08d
; A reasonable lower limit of 5d for physicality ... Assume line is resolved.
  if ~ keyword_set(siglim) then siglim=[5d,2000d]
  if ~ keyword_set(blrcomp) then blrcomp = -1

; Number of emission lines to fit
  nline = linelist->count()
  lines_arr = (linelist->keys())->toarray()
  
; Number of initial parameters before Gaussian parameters begin
  lratlim = 5 ; maximum number of line ratios to constrain
  ppoff0 = 3
  ppoff = ppoff0 + maxncomp*lratlim
  
  parinfo = REPLICATE({value:0d, fixed:0b, limited:[0B,0B], tied:'', $
                       limits:[0d,0d], step:0d, mpprint:0b, mpside:2, $
                       parname:'', line:'', comp:0d, sigmawave_tie:'', $
                       flux_tie:''}, $
                      ppoff+maxncomp*(nline*3))

; Number of initial parameters before Gaussian parameters begin
  parinfo[0].value = ppoff
  parinfo[0].fixed = 1B
  parinfo[0].parname = 'No. of non-Gaussian parameters'

; Maximum number of velocity components
  parinfo[1].value = maxncomp
  parinfo[1].fixed = 1B
  parinfo[1].parname = 'Maximum no. of velocity components'

; Spectral resolution
  parinfo[2].value = specres
  parinfo[2].fixed = 1B
  parinfo[2].parname = 'Spectral resolution in wavelength space [sigma]'

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
  if ncomp->haskey('Halpha') AND ncomp->haskey('[NII]6583') $
    then tmp_ncomp = ncomp['Halpha'] $
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
;    parinfo[ip1:ip2].limits  = rebin([0d,4d],2,tmp_ncomp)
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

;; Ha/Hb ratio
;  ilratlim = 3
;  lratlab = 'Ha/Hb'
;  if ncomp->haskey('Halpha') AND ncomp->haskey('Hbeta') $
;    then tmp_ncomp = ncomp['Halpha'] $
;  else tmp_ncomp = 0
;  if tmp_ncomp gt 0 then begin
;     ip1 = ppoff0 + ilratlim*maxncomp
;     ip2 = ip1 + tmp_ncomp - 1
;     fa = initflux['Halpha',0:tmp_ncomp-1]
;     fb = initflux['Hbeta',0:tmp_ncomp-1]
;     frat = dblarr(tmp_ncomp)+3d ; default if initial hb flux = 0
;     inz = where(fb gt 0,ctnz)
;     if ctnz gt 0 then frat[inz] = fa[inz]/fb[inz]
;     parinfo[ip1:ip2].value = frat
;     parinfo[ip1:ip2].limited = rebin([1B,0B],2,tmp_ncomp)
;;    Upper limit of 50 corresponds to E(B-V) = 2.89 using CCM
;     parinfo[ip1:ip2].limits  = rebin([2.86d,0d],2,tmp_ncomp)
;     parinfo[ip1:ip2].parname = 'Halpha/Hbeta line ratio'
;     parinfo[ip1:ip2].comp = indgen(tmp_ncomp)+1
;;    Check to see if line ratio is fixed
;     ilratfix = where(lratfix.keys() eq lratlab,ctlratfix)
;     for i=0,tmp_ncomp-1 do begin
;;       If line ratio is fixed, then fix it
;        lratfixed = 0b
;        if ctlratfix gt 0 then begin
;           if lratfix[lratlab,i] ne bad then begin
;              parinfo[ip1+i].value = lratfix[lratlab,i]
;              parinfo[ip1+i].fixed = 1b
;              parinfo[ip1+i].limited = [0b,0b]
;              lratfixed = 1b
;           endif
;        endif
;        if ~ lratfixed then begin
;;;          case of pegging at or exceeding upper limit
;;           if parinfo[ip1+i].value ge parinfo[ip1+i].limits[1] then $
;;              parinfo[ip1+i].value = parinfo[ip1+i].limits[1] - $
;;                                     (parinfo[ip1+i].limits[1] - $
;;                                      parinfo[ip1+i].limits[0])*0.1
;;          case of pegging at or dipping below lower limit
;           if parinfo[ip1+i].value le parinfo[ip1+i].limits[0] then $
;              parinfo[ip1+i].value = parinfo[ip1+i].limits[0] + $
;                                     parinfo[ip1+i].limits[0]*0.1d
;;                                     (parinfo[ip1+i].limits[1] - $
;;                                      parinfo[ip1+i].limits[0])*0.1
;        endif
;     endfor
;  endif

; [OII] ratio
; Limits from Pradhan et al. 2006, MNRAS, 366, L6
; 28aug2016, DSNR, changed limits to be more physically reasonable for AGN physics
; 28mar2019, DSNR, changed back to defaults
  ilratlim = 3
  lratlab = '[OII]3729/3726'
  if ncomp->haskey('[OII]3726') then tmp_ncomp = ncomp['[OII]3726'] $
  else tmp_ncomp=0
  if tmp_ncomp gt 0 then begin
    ip1 = ppoff0 + ilratlim*maxncomp
    ip2 = ip1+tmp_ncomp-1
    fa = initflux['[OII]3726',0:tmp_ncomp-1]
    fb = initflux['[OII]3729',0:tmp_ncomp-1]
    frat = dblarr(tmp_ncomp)+1d ; default if initial s2b flux = 0
    inz = where(fb gt 0,ctnz)
    if ctnz gt 0 then frat[inz] = fb[inz]/fa[inz]
    parinfo[ip1:ip2].value = frat
    parinfo[ip1:ip2].limited = rebin([1b,1b],2,tmp_ncomp)
    parinfo[ip1:ip2].limits  = rebin([0.35d,1.5d],2,tmp_ncomp)
;    parinfo[ip1:ip2].limits  = rebin([0.75d,1.4d],2,tmp_ncomp)
    parinfo[ip1:ip2].parname = '[OII]3729/3726 line ratio'
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


; MgII emission ratio
  ilratlim = 4
  lratlab = 'MgII2796/2803'
  if ncomp->haskey('MgII2803') then tmp_ncomp = ncomp['MgII2803'] $
  else tmp_ncomp=0
  if tmp_ncomp gt 0 then begin
    ip1 = ppoff0 + ilratlim*maxncomp
    ip2 = ip1+tmp_ncomp-1
    fa = initflux['MgII2796',0:tmp_ncomp-1]
    fb = initflux['MgII2803',0:tmp_ncomp-1]
    frat = dblarr(tmp_ncomp)+1d ; default if initial s2b flux = 0
    inz = where(fb gt 0,ctnz)
    if ctnz gt 0 then frat[inz] = fa[inz]/fb[inz]
    parinfo[ip1:ip2].value = frat
    parinfo[ip1:ip2].limited = rebin([1b,1b],2,tmp_ncomp)
    parinfo[ip1:ip2].limits  = rebin([1d,2d],2,tmp_ncomp)
;    parinfo[ip1:ip2].limits  = rebin([0.75d,1.4d],2,tmp_ncomp)
    parinfo[ip1:ip2].parname = 'MgII2796/2803 line ratio'
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
        if ((i+1 gt ncomp[line]) OR $
           (where(i+1 eq blrcomp) ge 0 AND $
            where(line eq blrlines) eq -1)) then begin 

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
                           format='(E0.6,A0,E0.6,A0,I0,A0)') 
                 parinfo[isoff].tied = $
                    string('P[',soff+indtie*3,']',format='(A0,I0,A0)')   
                 parinfo[iwoff].sigmawave_tie = linetie[line]
                 parinfo[isoff].sigmawave_tie = linetie[line]
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

;    the if statement here prevents MPFIT_TIE from issuing an IEEE exception,
;    since if we're not fitting [SII] then the ratio is set to 0 at the 
;    beginning of this routine
     if ncomp->haskey('[SII]6716') then begin
        if ncomp['[SII]6716'] gt 0 then begin
           ilratlim = 0
           linea = where(lines_arr eq '[SII]6716',cta)
           lineb  = where(lines_arr eq '[SII]6731',ctb)
           if cta gt 0 AND ctb gt 0 then begin
             parinfo[foff+lineb*3].tied = 'P['+$
                                          string(foff+linea*3,$
                                                 format='(I0)')+']/P['+$
                                          string(ppoff0+maxncomp*ilratlim+i,$
                                                 format='(I0)')+']'
             parinfo[foff+lineb*3].flux_tie = '[SII]6716'                    
           endif
        endif
     endif                                              

     ilratlim = 1
     linea = where(lines_arr eq '[NI]5198',cta)
     lineb = where(lines_arr eq '[NI]5200',ctb)
     if cta gt 0 AND ctb gt 0 then begin
       parinfo[foff+linea*3].tied = 'P['+$
                                    string(foff+lineb*3,$
                                           format='(I0)')+']/P['+$
                                    string(ppoff0+maxncomp*ilratlim+i,$
                                           format='(I0)')+']'
       parinfo[foff+linea*3].flux_tie = '[NI]5200'
     endif
                                         
     ilratlim = 2
     linea = where(lines_arr eq 'Halpha',cta)
     lineb = where(lines_arr eq '[NII]6583',ctb)
     if cta gt 0 AND ctb gt 0 then begin
       parinfo[foff+lineb*3].tied = 'P['+$
                                    string(foff+linea*3,$
                                           format='(I0)')+']*P['+$
                                    string(ppoff0+maxncomp*ilratlim+i,$
                                           format='(I0)')+']'
       parinfo[foff+lineb*3].flux_tie = 'Halpha'
     endif
     
;;    the if statement here prevents MPFIT_TIE from issuing an IEEE exception,
;;    since if we're not fitting Halpha then the ratio is set to 0 at the
;;    beginning of this routine
;     if ncomp->haskey('Halpha') then begin
;       if ncomp['Halpha'] gt 0 then begin
;         ilratlim = 3
;         linea = where(lines_arr eq 'Halpha',cta)
;         lineb = where(lines_arr eq 'Hbeta',ctb)
;         if cta gt 0 AND ctb gt 0 then begin
;           parinfo[foff+lineb*3].tied = 'P['+$
;                                        string(foff+linea*3,$
;                                               format='(I0)')+']/P['+$
;                                        string(ppoff0+maxncomp*ilratlim+i,$
;                                               format='(I0)')+']'
;           parinfo[foff+lineb*3].flux_tie = 'Halpha'
;         endif
;       endif
;     endif

     ilratlim = 3
     linea = where(lines_arr eq '[OII]3726',cta)
     lineb = where(lines_arr eq '[OII]3729',ctb)
     if cta gt 0 AND ctb gt 0 then begin
       parinfo[foff+lineb*3].tied = 'P['+$
         string(foff+linea*3,$
         format='(I0)')+']*P['+$
         string(ppoff0+maxncomp*ilratlim+i,$
         format='(I0)')+']'
       parinfo[foff+linea*3].flux_tie = '[OII]3729'
     endif


     ilratlim = 4
     linea = where(lines_arr eq 'MgII2796',cta)
     lineb = where(lines_arr eq 'MgII2803',ctb)
     if cta gt 0 AND ctb gt 0 then begin
       parinfo[foff+lineb*3].tied = 'P['+$
         string(foff+linea*3,$
         format='(I0)')+']/P['+$
         string(ppoff0+maxncomp*ilratlim+i,$
         format='(I0)')+']'
       parinfo[foff+linea*3].flux_tie = 'MgII2796'
     endif

     linea = where(lines_arr eq '[OIII]4959',cta)
     lineb = where(lines_arr eq '[OIII]5007',ctb)
     if cta gt 0 AND ctb gt 0 then begin
       parinfo[foff+linea*3].tied = 'P['+$
                                    string(foff+lineb*3,$
                                           format='(I0)')+']/3.0d'
       parinfo[foff+linea*3].flux_tie = '[OIII]5007'
;      Make sure initial value is correct
       parinfo[foff+linea*3].value = parinfo[foff+lineb*3].value/3.0d
     endif
                                         
     linea = where(lines_arr eq '[OI]6300',cta)
     lineb = where(lines_arr eq '[OI]6364',ctb)
     if cta gt 0 AND ctb gt 0 then begin
       parinfo[foff+lineb*3].tied = 'P['+$
                                    string(foff+linea*3,$
                                           format='(I0)')+']/3.0d'
       parinfo[foff+lineb*3].flux_tie = '[OI]6300'
;      Make sure initial value is correct
       parinfo[foff+lineb*3].value = parinfo[foff+linea*3].value/3.0d
     endif
     
     linea = where(lines_arr eq '[NII]6548',cta)
     lineb = where(lines_arr eq '[NII]6583',ctb)
     if cta gt 0 AND ctb gt 0 then begin
       parinfo[foff+linea*3].tied = 'P['+ $
                                    string(foff+lineb*3,$
                                           format='(I0)')+']/3.0d'
       parinfo[foff+linea*3].flux_tie = '[NII]6583'
;      Make sure initial value is correct
       parinfo[foff+linea*3].value = parinfo[foff+lineb*3].value/3.0d
    endif

    linea = where(lines_arr eq '[NeIII]3967',cta)
    lineb = where(lines_arr eq '[NeIII]3869',ctb)
    if cta gt 0 AND ctb gt 0 then begin
       parinfo[foff+linea*3].tied = 'P['+ $
          string(foff+lineb*3,$
          format='(I0)')+']/3.0d'
       parinfo[foff+linea*3].flux_tie = '[NeIII]3869'
       ;      Make sure initial value is correct
       parinfo[foff+linea*3].value = parinfo[foff+lineb*3].value/3.0d
    endif

    linea = where(lines_arr eq '[NeV]3345',cta)
    lineb = where(lines_arr eq '[NeV]3426',ctb)
    if cta gt 0 AND ctb gt 0 then begin
       parinfo[foff+linea*3].tied = 'P['+ $
          string(foff+lineb*3,$
          format='(I0)')+']/2.8d'
       parinfo[foff+linea*3].flux_tie = '[NeV]3426'
       ;      Make sure initial value is correct
       parinfo[foff+linea*3].value = parinfo[foff+lineb*3].value/2.8d
    endif
  endfor

; Check parinit initial values vs. limits
  badpar = where((parinfo.limited[0] AND $
                  parinfo.value lt parinfo.limits[0]) OR $
                 (parinfo.limited[1] AND $
                  parinfo.value gt parinfo.limits[1]),ct)
  if ct gt 0 then begin
     print,'Quantity','Line','Comp','Value','Lower limit','Upper limit',$
           format='(2A20,A5,3A15)'
     for i=0,ct-1 do begin
        j = badpar[i]
        print,parinfo[j].parname,parinfo[j].line,parinfo[j].comp,$
              parinfo[j].value,parinfo[j].limits[0],$
              parinfo[j].limits[1],format='(2A20,I5,3E15.6)'
     endfor
     message,'Initial values are outside limits.'
  endif else begin
     return,parinfo
  endelse

end
