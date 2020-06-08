; docformat = 'rst'
;
;+
;
; Convert output of MPFIT, with best-fit line parameters in a single
; array, into a structure with separate arrays for different line
; parameters. Compute total line fluxes from the best-fit line
; parameters.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    A structure with separate hashes for different line parameters. The hashes 
;    are indexed by line, and each value is an array over components. 
;    Tags: flux, fluxerr, fluxpk, fluxpkerr, nolines, wave, and sigma.
;
; :Params:
;    linelist: in, required, type=hash(lines)
;      List of emission line rest-frame wavelengths.
;    param: in, required, type=dblarr(N)
;      Best-fit parameter array output by MPFIT.
;    perror: in, optional, type=dblarr(N)
;      Errors in best fit parameters, output by MPFIT.
;    parinfo: in, required, type=structure
;      Structure input into MPFIT. Each tag has N values, one per parameter. 
;      Used to sort param and perror arrays.
;
; :Keywords:
;    waveran: in, optional, type=dblarr(2)
;      Set to upper and lower limits to return line parameters only
;      for lines within the given wavelength range. Lines outside this
;      range have fluxes set to 0.
;    tflux: out, optional, type=structure
;      A structure with separate hashes for total flux and error. The hashes 
;      are indexed by line. Tags: tflux, tfluxerr
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
;      2009may26, DSNR, created
;      2009jun07, DSNR, added error propagation and rewrote
;      2013nov01, DSNR, added documentation
;      2013nov25, DSNR, renamed, added copyright and license
;      2014jan13, DSNR, re-written to use hashes rather than arrays
;      2014feb26, DSNR, replaced ordered hashes with hashes
;      2015sep20, DSNR, compute total line flux and error
;      2016sep26, DSNR, account for new treatment of spectral resolution;
;                       fix flux errors for tied lines
;      2016oct05, DSNR, include wavelength and sigma errors
;      2016oct08, DSNR, turned off Ha/Hb limits b/c of issues with estimating
;                       Hbeta error in a noisy spectrum when pegged at lower
;                       limit
;      2016oct10, DSNR, added option to combine doublets; changed calculation
;                       of error when line ratio pegged
;      2020may11, DSNR, bug fix? set outstr = {nolines:1}
;      2020jun05, DSNR, added MgII case to doublets
;    
; :Copyright:
;    Copyright (C) 2013--2020 David S. N. Rupke
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
function ifsf_sepfitpars,linelist,param,perror,parinfo,waveran=waveran,$
                         tflux=tflux,doublets=doublets
                         
;  Return 0 if no lines were fit
   if n_elements(param) eq 1 then begin
      outstr = {nolines:1} ; bug fix? Changed 0 to 1; DSNR 2020may11
      goto,nolines
   endif

   c = 299792.458d
   maxncomp = param[1]

;  If this turns out to be slow, can try creating the hashes by concatenation 
;  instead.
   flux = hash(linelist->keys())
   fluxerr = hash(linelist->keys())
   fluxpk = hash(linelist->keys())
   fluxpkerr = hash(linelist->keys())
   fluxpk_obs = hash(linelist->keys())
   fluxpkerr_obs = hash(linelist->keys())
   sigma = hash(linelist->keys())
   sigmaerr = hash(linelist->keys())
   sigma_obs = hash(linelist->keys())
   sigmaerr_obs = hash(linelist->keys())
   wave = hash(linelist->keys())
   waveerr = hash(linelist->keys())

   if keyword_set(tflux) then begin
      tf = hash(linelist->keys())
      tfe = hash(linelist->keys())
   endif

   in2ha = where(parinfo.parname eq '[NII]/Halpha line ratio',ctn2ha)
   in1rat = where(parinfo.parname eq '[NI]5200/5198 line ratio',ctn1rat)
   is2rat = where(parinfo.parname eq '[SII]6716/6731 line ratio',cts2rat)
;   ihahb = where(parinfo.parname eq 'Halpha/Hbeta line ratio',cthahb)
   io2rat = where(parinfo.parname eq '[OII]3729/3726 line ratio',cto2rat)
   img2rat = where(parinfo.parname eq 'MgII2796/2803 line ratio',ctmg2rat)

;  Populate hashes
   foreach line,linelist->keys() do begin
;     indices
      iline = where(parinfo.line eq line)
      ifluxpk = cgsetintersection(iline,where(parinfo.parname eq 'flux_peak'))
      isigma = cgsetintersection(iline,where(parinfo.parname eq 'sigma'))
      iwave = cgsetintersection(iline,where(parinfo.parname eq 'wavelength'))
;     indices for errors, which is more complicated if error comes from line
;     to which this one is tied
      sigmawave_tie = parinfo[isigma[0]].sigmawave_tie ; line tied to
      if sigmawave_tie eq '' then begin
         isigmaerr = isigma
         iwaveerr = iwave
      endif else begin
         ilinetie = where(parinfo.line eq sigmawave_tie)
         isigmaerr = $
            cgsetintersection(ilinetie,where(parinfo.parname eq 'sigma'))
         iwaveerr = $
            cgsetintersection(ilinetie,where(parinfo.parname eq 'wavelength'))
      endelse

      wave[line] = param[iwave]
      sigma[line] = param[isigma]
      fluxpk[line] = param[ifluxpk]
      sigma_obs[line] = param[isigma]
      fluxpk_obs[line] = param[ifluxpk]
;     This bit of jujitsu ensures waveerr doesn't go NaN ...
      p1tmp = perror[iwaveerr]
      p2tmp = param[iwaveerr]
      p3tmp = param[iwave]
      inz = where(p1tmp ne 0d AND p2tmp ne 0d AND p3tmp ne 0d,ctnz)
      waveerrtmp = dblarr(n_elements(p1tmp))
      if ctnz gt 0 then waveerrtmp[inz] = p1tmp[inz]/p2tmp[inz]*p3tmp[inz]
      waveerr[line] = waveerrtmp
      sigmaerr[line] = perror[isigmaerr]
      sigmaerr_obs[line] = perror[isigmaerr]
      fluxpkerr[line] = perror[ifluxpk]
      fluxpkerr_obs[line] = perror[ifluxpk]

;     Because of the way these lines are tied to others (with a division!) they
;     can yield NaNs in components that aren't fit. Correct this.
;      if line eq '[SII]6731' OR line eq 'Hbeta' OR line eq '[NI]5189' then begin
      if line eq '[SII]6731' OR line eq '[NI]5189' then begin
         inan = where(~finite(fluxpk[line]),ctnan)
         if ctnan gt 0 then begin
            fluxpk[line,inan] = 0d
            fluxpkerr[line,inan] = 0d
            fluxpk_obs[line,inan] = 0d
            fluxpkerr_obs[line,inan] = 0d
         endif
      endif
   endforeach

;  Loop twice over lines to be able to do flux error corrections below.
   foreach line,linelist->keys() do begin

;     Fix flux errors associated with line ratios. E.g., [NII]/Halpha is a fitted
;     parameter and [NII]6583 is tied to it, so the formal error in [NII]6583
;     flux is 0. Add errors in Halpha and [NII]/Halpha in quadrature to get
;     error in [NII]6583.
      if line eq '[NII]6583' AND ctn2ha gt 0 then begin
         fluxpkerr_obs[line,0:ctn2ha-1] = $
            fluxpk_obs[line,0:ctn2ha-1]*sqrt($
            (perror[in2ha]/param[in2ha])^2d + $
            (fluxpkerr_obs['Halpha',0:ctn2ha-1]/$
            fluxpk_obs['Halpha',0:ctn2ha-1])^2d)
;        In pegged case, set errors equal to each other
         ipegged = where(perror[in2ha] eq 0d AND param[in2ha] ne 0d,ctpegged)
         if ctpegged gt 0 then $
            fluxpkerr_obs['[NII]6583',ipegged] = $
               fluxpkerr_obs['Halpha',ipegged]
         fluxpkerr[line] = fluxpkerr_obs[line]
      endif
      if line eq '[SII]6731' AND cts2rat gt 0 then begin
         fluxpkerr_obs[line,0:cts2rat-1] = $
            fluxpk_obs[line,0:cts2rat-1]*sqrt($
            (perror[is2rat]/param[is2rat])^2d + $
            (fluxpkerr_obs['[SII]6716',0:cts2rat-1]/$
            fluxpk_obs['[SII]6716',0:cts2rat-1])^2d)
;        In pegged case, set errors equal to each other
         ipegged = where(perror[is2rat] eq 0d AND param[is2rat] ne 0d,ctpegged)
         if ctpegged gt 0 then $
            fluxpkerr_obs['[SII]6731',ipegged] = $
               fluxpkerr_obs['[SII]6716',ipegged]
         fluxpkerr[line] = fluxpkerr_obs[line]
      endif
      if line eq '[NI]5198' AND ctn1rat gt 0 then begin
         fluxpkerr_obs[line,0:ctn1rat-1] = $
            fluxpk_obs[line,0:ctn1rat-1]*sqrt($
            (perror[in1rat]/param[in1rat])^2d + $
            (fluxpkerr_obs['[NI]5200',0:ctn1rat-1]/$
            fluxpk_obs['[NI]5200',0:ctn1rat-1])^2d)
         fluxpkerr[line] = fluxpkerr_obs[line]
;        In pegged case, set errors equal to each other
         ipegged = where(perror[in1rat] eq 0d AND param[in1rat] ne 0d,ctpegged)
         if ctpegged gt 0 then $
            fluxpkerr_obs['[NI]5198',ipegged] = $
               fluxpkerr_obs['[NI]5200',ipegged]
         fluxpkerr[line] = fluxpkerr_obs[line]
      endif
;      if line eq 'Hbeta' AND cthahb gt 0 then begin
;         fluxpkerr_obs[line,0:cthahb-1] = $
;            fluxpk_obs[line,0:cthahb-1]*sqrt($
;            (perror[ihahb]/param[ihahb])^2d + $
;            (fluxpkerr_obs['Halpha',0:cthahb-1]/$
;            fluxpk_obs['Halpha',0:cthahb-1])^2d)
;         fluxpkerr[line] = fluxpkerr_obs[line]
;;        If Halpha/Hbeta gets too high, MPFIT sees it as an "upper limit" and
;;        sets perror = 0d. Then the errors seem too low, and it registers as a
;;        detection. This bit of code corrects that.
;         ipeggedupper = $
;            where(param[ihahb] gt 1d1 AND perror[ihahb] eq 0d,ctpegged)
;         if ctpegged gt 0 then begin
;            fluxpk[line,ipeggedupper] = 0d
;            fluxpk_obs[line,ipeggedupper] = 0d
;         endif
      if line eq 'Hbeta' AND linelist.haskey('Halpha') then begin
;        If Halpha/Hbeta goes belowlower limit, then we re-calculate the errors
;        add discrepancy in quadrature to currently calculated error. Assume
;        error in fitting is in Hbeta and adjust accordingly.
         fha = fluxpk['Halpha']
         ihahb = where(fluxpk['Halpha'] gt 0d AND fluxpk['Hbeta'] gt 0d,cthahb)
         if cthahb gt 0 then begin
            itoolow = $
               where(fluxpk['Halpha',ihahb]/fluxpk['Hbeta'] lt 2.86d,cttoolow)
            if cttoolow gt 0 then begin
               fluxpkdiff = fluxpk[line,itoolow] - fluxpk['Halpha',itoolow]/2.86d
               fluxpk[line,itoolow] -= fluxpkdiff
               fluxpk_obs[line,itoolow] -= fluxpkdiff
               fluxpkerr[line,itoolow] = $
                  sqrt(fluxpkerr[line,itoolow]^2d + fluxpkdiff^2d)
               fluxpkerr_obs[line,itoolow] = $
                  sqrt(fluxpkerr_obs[line,itoolow]^2d + fluxpkdiff^2d)
            endif
         endif
      endif
      if line eq '[OII]3729' AND cto2rat gt 0 then begin
         fluxpkerr_obs[line,0:cto2rat-1] = $
            fluxpk_obs[line,0:cto2rat-1]*sqrt($
            (perror[io2rat]/param[io2rat])^2d + $
            (fluxpkerr_obs['[OII]3726',0:cto2rat-1]/$
            fluxpk_obs['[OII]3726',0:cto2rat-1])^2d)
;        In pegged case, set errors equal to each other
         ipegged = where(perror[io2rat] eq 0d AND param[io2rat] ne 0d,ctpegged)
         if ctpegged gt 0 then $
            fluxpkerr_obs['[OII]3729',ipegged] = $
               fluxpkerr_obs['[OII]3726',ipegged]
         fluxpkerr[line] = fluxpkerr_obs[line]
      endif

      if line eq 'MgII2803' AND ctmg2rat gt 0 then begin
        fluxpkerr_obs[line,0:ctmg2rat-1] = $
          fluxpk_obs[line,0:ctmg2rat-1]*sqrt($
          (perror[img2rat]/param[img2rat])^2d + $
          (fluxpkerr_obs['MgII2796',0:ctmg2rat-1]/$
          fluxpk_obs['MgII2796',0:ctmg2rat-1])^2d)
        ;        In pegged case, set errors equal to each other
        ipegged = where(perror[img2rat] eq 0d AND param[img2rat] ne 0d,ctpegged)
        if ctpegged gt 0 then $
          fluxpkerr_obs['MgII2803',ipegged] = $
          fluxpkerr_obs['MgII2796',ipegged]
        fluxpkerr[line] = fluxpkerr_obs[line]
      endif

;     Add back in spectral resolution
;     Can't use sigma = 0 as criterion since the line could be fitted but unresolved.
      sigmatmp = sigma[line]/c*wave[line]
      inz = where(fluxpk[line] gt 0,ctnz)
      if ctnz gt 0 then begin
;        Make sure we're not adding something to 0 -- i.e. the component wasn't fit.
         sigmatmp[inz] = sqrt(sigmatmp[inz]^2d + param[2]^2d)
         sigma_obs[line,inz] = sigmatmp[inz]/wave[line,inz]*c ; in km/s
;        error propagation for adding in quadrature
         sigmaerr_obs[line,inz] *= $
            sigma[line,inz]/c*wave[line,inz] / sigmatmp[inz]
;        Correct peak flux and error for deconvolution
         fluxpk[line,inz] *= sigma_obs[line,inz]/sigma[line,inz]
         fluxpkerr[line,inz] *= sigma_obs[line,inz]/sigma[line,inz]
      endif

;     Compute total Gaussian flux
;     sigma and error need to be in wavelength space
      gflux = IFSF_GAUSSFLUX(fluxpk_obs[line],sigmatmp,$
                             normerr=fluxpkerr_obs[line],$
                             sigerr=sigmaerr_obs[line]/c*wave[line])
      flux[line] = gflux.flux
      fluxerr[line] = gflux.flux_err

;     Set fluxes to 0 outside of wavelength range, or if NaNs or infinite errors
      if keyword_set(waveran) then begin
         inoflux = where(waveran[0] gt wave[line]*(1 - 3d*sigma[line]/c) OR $
                         waveran[1] lt wave[line]*(1 + 3d*sigma[line]/c) OR $
                         finite(fluxerr[line]) eq 0 OR $
                         finite(fluxpkerr[line]) eq 0,ct)
         if ct gt 0 then begin
            flux[line,inoflux] = 0d
            fluxerr[line,inoflux] = 0d
            fluxpk[line,inoflux] = 0d
            fluxpkerr[line,inoflux] = 0d
            fluxpk_obs[line,inoflux] = 0d
            fluxpkerr_obs[line,inoflux] = 0d
         endif
      endif

;     Compute total fluxes summed over components
      igd = where(flux[line] gt 0d,ctgd)
      if keyword_set(tflux) then begin
         if ctgd gt 0 then begin
            tf[line] = total(flux[line,igd])
            tfe[line] = sqrt(total(fluxerr[line,igd]^2d))
         endif else begin
            tf[line] = 0d
            tfe[line] = 0d
         endelse
      endif
   endforeach
   
;  Special doublet cases: combine fluxes from each line
   if keyword_set(doublets) then begin
      sdoub = size(doublets)
      if sdoub[0] eq 1 then ndoublets = 1 else ndoublets = sdoub[2]
      for i=0,ndoublets-1 do begin
         if linelist.haskey(doublets[0,i]) AND linelist.haskey(doublets[1,i]) then begin
;           new line label
            dkey = doublets[0,i]+'+'+doublets[1,i]
;           add fluxes
            tf[dkey] = tf[doublets[0,i]]+tf[doublets[1,i]]
            flux[dkey] = flux[doublets[0,i]]+flux[doublets[1,i]]
            fluxpk[dkey] = fluxpk[doublets[0,i]]+fluxpk[doublets[1,i]]
            fluxpk_obs[dkey] = fluxpk_obs[doublets[0,i]]+fluxpk_obs[doublets[1,i]]
;           add flux errors in quadrature
            tfe[dkey] = sqrt(tfe[doublets[0,i]]^2d +$
                             tfe[doublets[1,i]]^2d)
            fluxerr[dkey] = sqrt(fluxerr[doublets[0,i]]^2d +$
                                 fluxerr[doublets[1,i]]^2d)
            fluxpkerr[dkey] = sqrt(fluxpkerr[doublets[0,i]]^2d +$
                                   fluxpkerr[doublets[1,i]]^2d)
            fluxpkerr_obs[dkey] = sqrt(fluxpkerr_obs[doublets[0,i]]^2d +$
                                       fluxpkerr_obs[doublets[1,i]]^2d)
;           average waves and sigmas and errors
            wave[dkey] = (wave[doublets[0,i]]+wave[doublets[1,i]])/2d
            waveerr[dkey] = (waveerr[doublets[0,i]]+waveerr[doublets[1,i]])/2d
            sigma[dkey] = (sigma[doublets[0,i]]+sigma[doublets[1,i]])/2d
            sigmaerr[dkey] = (sigmaerr[doublets[0,i]]+sigmaerr[doublets[1,i]])/2d
            sigma_obs[dkey] = (sigma_obs[doublets[0,i]]+sigma_obs[doublets[1,i]])/2d
            sigmaerr_obs[dkey] = (sigmaerr_obs[doublets[0,i]]+$
                                  sigmaerr_obs[doublets[1,i]])/2d
         endif
      endfor
   endif
   
   outstr = {nolines:0,$
             flux:flux,fluxerr:fluxerr,$
             fluxpk:fluxpk,fluxpkerr:fluxpkerr,$
             wave:wave,waveerr:waveerr,$
             sigma:sigma,sigmaerr:sigmaerr,$
             sigma_obs:sigma_obs,sigmaerr_obs:sigmaerr_obs,$
             fluxpk_obs:fluxpk_obs,fluxpkerr_obs:fluxpkerr_obs}
   if keyword_set(tflux) then tflux = {tflux:tf,tfluxerr:tfe}

nolines:

   return,outstr

end
