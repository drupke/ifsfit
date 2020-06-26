;
;+
;
; This procedure makes maps of various quantities. Contains one
; helper routine: IFSF_PA.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Postscript plots.
;
; :Params:
;    initproc: in, required, type=string
;      Name of procedure to initialize the fit.
;
; :Keywords:
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
;      2014jan24, DSNR, created
;      2014apr15, DSNR, moved colorbar division code to subroutine IFSF_CBDIV;
;                       added ability to automatically fix ranges separately
;                       for different components; cleared up some floating point
;                       errors
;      2014apr21, DSNR, added line ratio maps and VO plots
;      2014may23, DSNR, added NaD maps and continuum images
;      2014jun02, DSNR, allow use without previous emission-line fit from IFSF
;      2014jun04, DSNR, plots velocities using cumulative velocity distributions
;      2014julXY, DSNR  calls local plotting routine if requested;
;                       plots no. of components in NaD fit
;                       new helper subroutines for plotting axes and finding ranges
;      2014jul24, DSNR, plots NaD fitted velocities
;      2014dec07, DSNR, added compass to plots
;      2015may15, DSNR, fixed compass rose bug (apparent if aspect ratio of images
;                       differs much from 1)
;      2015jun03, DSNR, added IF statements to deal with case of fitted NaD
;                       absorption but no fitted NaD emission; added SNR threshold
;                       based on Weq to mapping of fitted NaD properties, not
;                       just empirical properties as was previously coded
;      2015jun08, DSNR, added peak fitting to continuum data in *cont.eps and
;                       *cont_rad.eps plots
;      2015sep21, DSNR, big changes to dereddening procedures, and other changes
;                       to line plotting procedures
;      2016jan24, DSNR, added 'diskline' and 'ofparline' tags to INITMAPS
;      2016feb04, DSNR, fixed treatment of HST image sizes
;      2016feb15, DSNR, fixed factor-of-2 error in PSF FWHM estimates
;      2018oct09, DSNR, added BUFFAC option to IFSF_HSTSUBIM calls
;    
; :Copyright:
;    Copyright (C) 2014--2018 David S. N. Rupke
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
function ifsf_pa,xaxis,yaxis

;  Computes PA east of north, assuming +y axis is north

;  assume input variables are 0-,1-,or 2-dimensional
   sizetmp = size(xaxis)
   if sizetmp[0] eq 0 then pa = dblarr(1) $
   else if sizetmp[0] eq 1 then pa = dblarr(sizetmp[1]) $
   else pa = dblarr(sizetmp[1],sizetmp[2])

   iul = where(xaxis le 0 AND yaxis gt 0,ctul)
   ill = where(xaxis lt 0 AND yaxis le 0,ctll)
   ilr = where(xaxis ge 0 AND yaxis lt 0,ctlr)
   iur = where(xaxis gt 0 AND yaxis ge 0,ctur)
   if ctul gt 0 then pa[iul] = atan(abs(xaxis[iul])/yaxis[iul])/!DPi*180d
   if ctll gt 0 then pa[ill] = 90d + atan(abs(yaxis[ill])/abs(xaxis[ill]))/!DPi*180d
   if ctlr gt 0 then pa[ilr] = 180d + atan(xaxis[ilr]/abs(yaxis[ilr]))/!DPi*180d
   if ctur gt 0 then pa[iur] = 270d + atan(yaxis[iur]/xaxis[iur])/!DPi*180d

   if sizetmp[0] eq 0 then pa = pa[0]
   return,pa
end
;-------------------------------------------------------------------------------
pro ifsf_makemaps,initproc

   fwhm2sig = 2d*sqrt(2d*alog(2d))
   plotquantum = 2.5d ; in inches
   bad = 1d99
   c_kms = 299792.458d
   ncbdivmax = 7
   maxnadabscomp = 3
   maxnademcomp = 3
   taumax = 5d

;  physical constants
   mump = 1.4d*1.672649d-24                      
                                ; (log) mass per particle, in grams
   mumpsm = 1.4d*1.672649d-24 / 1.989d33
                                ; (log) mass per particle, in solar masses
   speryr = 24d*3600d*365.25d    ; seconds in a year
   mperpc = 3.0856d16            ; meters in a parsec
   lsun = 3.826d33               ; solar luminosities, erg/s
   msun = 1.989e33               ; solar mass, g
   c_cms = 2.99792d10
   ionfrac = 0.9d                      ; Na ionization fraction
   oneminusionfrac_relerr = [0d,0d]
   naabund = -5.69d                    ; Na abundance
   nadep = 0.95d                    ; Na depletion
   volemis = 2.63d-25            ; volume emissivity of Ha = product of recomb. coeff. 
                                 ; and photon energy; units erg cm^3 s^-1
   elecden_default = 100d                 ; electron density, cm^-3
   elecden_err_default = [50d,100d]

   lineofdashes = strjoin(replicate('-',62))
   
;  factor by which to resample images for PS-to-PDF conversion
   samplefac = 10
   resampthresh = 500

;  Values for computing electron density from [SII] ratio
;  from Sanders, Shapley, et al. 2015
   s2_minrat = 0.4375d
   s2_maxrat = 1.4484d
   s2_a = 0.4315d
   s2_b = 2107d
   s2_c = 627.1d
;   s2_maxden = (s2_c * s2_minrat - s2_a*s2_b)/(s2_a - s2_minrat)
;   s2_minden = (s2_c * s2_maxrat - s2_a*s2_b)/(s2_a - s2_maxrat)
   s2_maxden = 1d4
   s2_minden = 1d1
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Load initialization parameters and line data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  Get initialization structures
   initmaps={dumy: 1}
   initdat=call_function(initproc,initmaps=initmaps)
   if tag_exist(initdat,'donad') then begin
      initnad={dumy: 1}
      initdat=call_function(initproc,initnad=initnad)
   endif

;  Get galaxy-specific parameters from initialization file
   center_axes = -1
   center_nuclei = -1
   if tag_exist(initmaps,'center_axes') then $
      center_axes = initmaps.center_axes
   if tag_exist(initmaps,'center_nuclei') then $
      center_nuclei = initmaps.center_nuclei

;  Get linelist
   if ~ tag_exist(initdat,'noemlinfit') then begin
      linelabels=1b
      if tag_exist(initdat,'argslinelist') then $
         linelist = ifsf_linelist(initdat.lines,linelab=linelabels,$
                                  _extra=initdat.argslinelist) $
      else $
         linelist = ifsf_linelist(initdat.lines,linelab=linelabels)
;     Linelist with doublets to combine
     emldoublets = [['[SII]6716','[SII]6731'],$
                    ['[OII]3726','[OII]3729'],$
                    ['[NI]5198','[NI]5200'],$
                    ['[NeIII]3869','[NeIII]3967'],$
                    ['[NeV]3345','[NeV]3426'],$
                    ['MgII2796','MgII2803']]
      sdoub = size(emldoublets)
      if sdoub[0] eq 1 then ndoublets = 1 else ndoublets = sdoub[2]
      lines_with_doublets = initdat.lines
      for i=0,ndoublets-1 do begin
         if linelist.haskey(emldoublets[0,i]) AND $
            linelist.haskey(emldoublets[1,i]) then begin
            dkey = emldoublets[0,i]+'+'+emldoublets[1,i]
            lines_with_doublets = [lines_with_doublets,dkey]
         endif
      endfor
      if tag_exist(initdat,'argslinelist') then $
         linelist_with_doublets = $
            ifsf_linelist(lines_with_doublets,linelab=linelabels,$
                          _extra=initdat.argslinelist) $
      else $
         linelist_with_doublets = $
            ifsf_linelist(lines_with_doublets,linelab=linelabels)

   endif
   if tag_exist(initdat,'donad') then $
      if tag_exist(initdat,'argslinelist') then $
         nadlinelist = ifsf_linelist(['NaD1','NaD2','HeI5876'],$
                                     _extra=initdat.argslinelist) $
      else $
         nadlinelist = ifsf_linelist(['NaD1','NaD2','HeI5876'])

;  Get range file
;
;  plot types, in order; used for correlating with input ranges (array 
;  rangequant)
   hasrangefile=0
   if tag_exist(initmaps,'rangefile') then begin
      if file_test(initmaps.rangefile) then begin
         readcol,initmaps.rangefile,rangeline,rangequant,rangelo,rangehi,$
                 rangencbdiv,format='(A,A,D,D,I)',/silent
         hasrangefile=1
      endif else message,'Range file listed in INITMAPS but not found.'
   endif

;  Restore line maps
   if ~ tag_exist(initdat,'noemlinfit') then $
      restore,file=initdat.outdir+initdat.label+'.lin.xdr'
   
;  Restore continuum parameters
;   if tag_exist(initdat,'decompose_ppxf_fit') OR $
;      tag_exist(initdat,'decompose_qso_fit') then $
   restore,file=initdat.outdir+initdat.label+'.cont.xdr'

;  Get NaD parameters
   if tag_exist(initdat,'donad') then begin
      restore,file=initdat.outdir+initdat.label+'.nadspec.xdr'
      restore,file=initdat.outdir+initdat.label+'.nadfit.xdr'
      if tag_exist(initmaps,'badnademp') then begin
         tagstobad=['WEQ','IWEQ','EMFLUX','EMUL','VEL']
         tagnames=tag_names(nadcube)
         ibad = where(initmaps.badnademp eq 1b,ctbad)
         if ctbad gt 0 then begin
            for i=0,n_elements(tagstobad)-1 do begin
               itag = where(tagnames eq tagstobad[i])
               sizetag = size(nadcube.(itag))
               for j=0,sizetag[3]-1 do begin
                  tmp=nadcube.(itag)[*,*,j]
                  tmp[ibad]=bad
                  nadcube.(itag)[*,*,j]=tmp
               endfor
            endfor
         endif
      endif
   endif
   
   if ~ tag_exist(initdat,'donad') AND $
      tag_exist(initmaps,'noemlinfit') then begin
      message,'No emission line or absorption line data specified.'
   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Compute some things
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  Luminosity and angular size distances
   if tag_exist(initdat,'distance') then begin
      ldist = initdat.distance
      asdist = ldist/(1d + initdat.zsys_gas)^2d
   endif else begin
; Planck 2018 parameters: https://ui.adsabs.harvard.edu/#abs/arXiv:1807.06209
      ldist = lumdist(initdat.zsys_gas,H0=67.4d,Omega_m=0.315d,Lambda0=0.685d,/silent)
      asdist = ldist/(1d + initdat.zsys_gas)^2d
   endelse
   kpc_per_as = asdist*1000d/206265d
   kpc_per_pix = initdat.platescale * kpc_per_as


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Load and process continuum data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  Data cube
   if not tag_exist(initdat,'datext') then datext=1 else datext=initdat.datext
   if not tag_exist(initdat,'varext') then varext=2 else varext=initdat.varext
   if not tag_exist(initdat,'dqext') then dqext=3 else dqext=initdat.dqext
   header=1
   datacube = ifsf_readcube(initdat.infile,/quiet,oned=oned,header=header,$
                            datext=datext,varext=varext,dqext=dqext)
   if tag_exist(initmaps,'fluxfactor') then begin
      datacube.dat *= initmaps.fluxfactor
      datacube.var *= (initmaps.fluxfactor)^2d
   endif
   size_tmp = size(datacube.dat)
   dx = size_tmp[1]
   dy = size_tmp[2]
   dz = size_tmp[3]
;  defined so that center of spaxel at bottom left has coordinates [1,1] 
   if center_axes[0] eq -1 then center_axes = [double(dx)/2d,double(dy)/2d]+0.5d
   if center_nuclei[0] eq -1 then center_nuclei = center_axes
   if tag_exist(initmaps,'vornorm') then begin
      datacube.dat /= rebin(initmaps.vornorm,dx,dy,dz)
   endif

;  Image window
   if tag_exist(initmaps,'plotwin') then begin
      plotwin = initmaps.plotwin
      dxwin = plotwin[2]-plotwin[0]+1
      dywin = plotwin[3]-plotwin[1]+1
   endif else begin
      plotwin = [1,1,dx,dy]
      dxwin = dx
      dywin = dy
   endelse
   
;  Figure aspect ratio multiplier
   if tag_exist(initmaps,'aspectrat') then aspectrat = initmaps.aspectrat $
   else aspectrat = 1d

;  HST data
   dohst=0
   dohstbl=0
   dohstrd=0
   dohstsm=0
   dohstcol=0
   dohstcolsm=0
   if tag_exist(initmaps,'hst') AND tag_exist(initmaps,'hstbl') then begin
      dohstbl=1
      if tag_exist(initmaps.hstbl,'ext') then hstblext=initmaps.hstbl.ext $
      else hstblext = 1
      hstbl = readfits(initmaps.hstbl.file,hstblhead,/silent,$
                       ext=hstblext)
      hst_big_ifsfov = dblarr(4,2)
      if tag_exist(initmaps.hstbl,'platescale') then $
         hstpsbl = initmaps.hstbl.platescale $
      else hstpsbl = 0.05d
      if tag_exist(initmaps.hstbl,'refcoords') then $
         hstrefcoords = initmaps.hstbl.refcoords $
      else hstrefcoords = initmaps.hst.refcoords
      if tag_exist(initmaps.hstbl,'buffac') then $
         hstbl_buffac = initmaps.hstbl.buffac $
      else hstbl_buffac = 2d
      if tag_exist(initmaps.hst,'subim_sm') AND $
         tag_exist(initmaps.hstbl,'sclargs_sm') then begin
         hst_sm_ifsfov = dblarr(4,2)
         bhst_sm = ifsf_hstsubim(hstbl,[initmaps.hst.subim_sm,$
                                 initmaps.hst.subim_sm],$
                                 [dx,dy],initdat.platescale,$
                                 initdat.positionangle,center_nuclei,$
                                 hstrefcoords,$
                                 initmaps.hstbl.scllim,$
                                 sclargs=initmaps.hstbl.sclargs_sm,$
                                 ifsbounds=hst_sm_ifsfov,hstps=hstpsbl)
      endif
      bhst_big = ifsf_hstsubim(hstbl,[initmaps.hst.subim_big,$
                               initmaps.hst.subim_big],$
                               [dx,dy],initdat.platescale,$
                               initdat.positionangle,center_nuclei,$
                               hstrefcoords,$
                               initmaps.hstbl.scllim,$
                               sclargs=initmaps.hstbl.sclargs_big,$
                               ifsbounds=hst_big_ifsfov,hstps=hstpsbl)
      bhst_fov = ifsf_hstsubim(hstbl,[0,0],[dx,dy],initdat.platescale,$
                               initdat.positionangle,center_nuclei,$
                               hstrefcoords,$
                               initmaps.hstbl.scllim,$
                               sclargs=initmaps.hstbl.sclargs_fov,$
                               /fov,hstps=hstpsbl,buffac=hstbl_buffac)
      bhst_fov_ns = ifsf_hstsubim(hstbl,[0,0],[dx,dy],initdat.platescale,$
                                  initdat.positionangle,center_nuclei,$
                                  hstrefcoords,[0,0],/noscl,/fov,$
                                  hstps=hstpsbl,buffac=hstbl_buffac)
      if tag_exist(initmaps,'hstblsm') then begin
         dohstsm=1

;        For F05189, mask central pixels before smoothing
         if initdat.label eq 'f05189' then begin
            size_tmp = size(hstbl)
            map_x_tmp = rebin(dindgen(size_tmp[1]),size_tmp[1],size_tmp[2])
            map_y_tmp = rebin(transpose(dindgen(size_tmp[2])),$
                              size_tmp[1],size_tmp[2])
            map_rkpc_tmp = sqrt((map_x_tmp - (hstrefcoords[0]+$
                                 initmaps.hstbl.nucoffset[0]-1))^2d + $
                                (map_y_tmp - (hstrefcoords[1]+$
                                 initmaps.hstbl.nucoffset[1]-1))^2d) $
                           * initmaps.hstbl.platescale * kpc_per_as
            ipsf = where(map_rkpc_tmp le 0.15d)
            ipsf_bkgd = where(map_rkpc_tmp gt 0.15d AND map_rkpc_tmp le 0.25d)
            hstbl_tmp = hstbl
            hstbl_tmp[ipsf] = median(hstbl[ipsf_bkgd])
            hstblsm = filter_image(hstbl_tmp,fwhm=initmaps.hst.smoothfwhm,/all)
         endif else begin
            hstbltmp = hstbl
            ibadhst = where(hstbl eq 0d,ctbadhst)
            if ctbadhst gt 0 then hstbltmp[ibadhst] = !values.d_nan
            fwhm = initmaps.hst.smoothfwhm
            boxwidth = round((fwhm^2d)/2d + 1d)
            if not boxwidth then boxwidth++
;            hstblsm = filter_image(hstbl,fwhm=initmaps.hst.smoothfwhm,/all)
            hstblsm = filter_image(hstbltmp,smooth=boxwidth,/iter,/all)
            ibadhst = where(finite(hstblsm,/nan),ctbadhst)
            if ctbadhst gt 0 then hstblsm[ibadhst] = bad
            hstbltmp = 0
         endelse
                  
;         bhst_fov_sm = ifsf_hstsubim(hstblsm,[0,0],[dx,dy],$
;                                     initdat.platescale,$
;                                     initdat.positionangle,center_nuclei,$
;                                     hstrefcoords,$
;                                     initmaps.hstblsm.scllim,$
;                                     sclargs=initmaps.hstblsm.sclargs,$
;                                     /fov,hstps=hstpsbl,buffac=hstbl_buffac)
         bhst_fov_sm_ns= ifsf_hstsubim(hstblsm,[0,0],[dx,dy],$
                                       initdat.platescale,$
                                       initdat.positionangle,center_nuclei,$
                                       hstrefcoords,[0,0],/noscl,$
                                       /fov,hstps=hstpsbl,buffac=hstbl_buffac)
         bhst_fov_sm_ns_rb = congrid(bhst_fov_sm_ns,dx,dy,/interp,/center)
         if tag_exist(initdat,'vormap') then begin
            tmp_rb = bhst_fov_sm_ns_rb
            nvor = max(initdat.vormap)
            badvor = where(~ finite(initdat.vormap),ctbad)
            if ctbad gt 0 then tmp_rb[badvor] = 0d
            for i=1,nvor do begin
               ivor = where(initdat.vormap eq i,ctbins)
               tmp_rb[ivor] = total(bhst_fov_sm_ns_rb[ivor])/ctbins
            endfor
            bhst_fov_sm_ns_rb = tmp_rb
         endif
         if tag_exist(initmaps.hstblsm,'bgsub') then $
            bhst_fov_sm_ns_rb -= initmaps.hstblsm.bgsub
      endif
   endif      
   if tag_exist(initmaps,'hst') AND tag_exist(initmaps,'hstrd') then begin
      dohstrd=1
      if tag_exist(initmaps.hstrd,'ext') then hstrdext=initmaps.hstrd.ext $
      else hstrdext = 1
      hstrd = readfits(initmaps.hstrd.file,hstrdhead,/silent,$
                       ext=hstrdext)
      hst_big_ifsfov = dblarr(4,2)
      if tag_exist(initmaps.hstrd,'refcoords') then $
         hstrefcoords = initmaps.hstrd.refcoords $
      else hstrefcoords = initmaps.hst.refcoords
      if tag_exist(initmaps.hstrd,'platescale') then $
         hstpsrd = initmaps.hstrd.platescale $
      else hstpsrd = 0.05d
      if tag_exist(initmaps.hstrd,'buffac') then $
         hstrd_buffac = initmaps.hstrd.buffac $
      else hstrd_buffac = 2d
      if tag_exist(initmaps.hst,'subim_sm') AND $
         tag_exist(initmaps.hstrd,'sclargs_sm') then begin
         hst_sm_ifsfov = dblarr(4,2)
         rhst_sm = ifsf_hstsubim(hstrd,[initmaps.hst.subim_sm,$
                                        initmaps.hst.subim_sm],$
                                 [dx,dy],initdat.platescale,$
                                 initdat.positionangle,center_nuclei,$
                                 hstrefcoords,$
                                 initmaps.hstrd.scllim,$
                                 sclargs=initmaps.hstrd.sclargs_sm,$
                                 ifsbounds=hst_sm_ifsfov,hstps=hstpsrd)
      endif
      rhst_big = ifsf_hstsubim(hstrd,[initmaps.hst.subim_big,$
                               initmaps.hst.subim_big],$
                               [dx,dy],initdat.platescale,$
                               initdat.positionangle,center_nuclei,$
                                 hstrefcoords,$
                               initmaps.hstrd.scllim,$
                               sclargs=initmaps.hstrd.sclargs_big,$
                               ifsbounds=hst_big_ifsfov,hstps=hstpsrd)
      rhst_fov_sc = ifsf_hstsubim(hstrd,[0,0],[dx,dy],initdat.platescale,$
                               initdat.positionangle,center_nuclei,$
;                               0,center_nuclei,$
                                 hstrefcoords,$
                               initmaps.hstrd.scllim,$
                               sclargs=initmaps.hstrd.sclargs_fov,$
                               /fov,hstps=hstpsrd,buffac=hstrd_buffac)
      rhst_fov_ns = ifsf_hstsubim(hstrd,[0,0],[dx,dy],initdat.platescale,$
                                  initdat.positionangle,center_nuclei,$
;                                  0,center_nuclei,$
                                  hstrefcoords,[0,0],/noscl,/fov,$
                                  hstps=hstpsrd,buffac=hstrd_buffac)
      if tag_exist(initmaps,'hstrdsm') then begin
         dohstsm=1

;        For F05189, mask central pixels before smoothing
         if initdat.label eq 'f05189' then begin
            size_tmp = size(hstrd)
            map_x_tmp = rebin(dindgen(size_tmp[1]),size_tmp[1],size_tmp[2])
            map_y_tmp = rebin(transpose(dindgen(size_tmp[2])),$
                              size_tmp[1],size_tmp[2])
            map_rkpc_tmp = sqrt((map_x_tmp - (hstrefcoords[0]+$
                                 initmaps.hstrd.nucoffset[0]-1))^2d + $
                                (map_y_tmp - (hstrefcoords[1]+$
                                 initmaps.hstrd.nucoffset[1]-1))^2d) $
                           * initmaps.hstrd.platescale * kpc_per_as
            ipsf = where(map_rkpc_tmp le 0.15d)
            ipsf_bkgd = where(map_rkpc_tmp gt 0.15d AND map_rkpc_tmp le 0.25d)
            hstrd_tmp = hstrd
            hstrd_tmp[ipsf] = median(hstrd[ipsf_bkgd])
            hstrdsm = filter_image(hstrd_tmp,fwhm=initmaps.hst.smoothfwhm,/all)
         endif else begin
            hstrdtmp = hstrd
            ibadhst = where(hstrd eq 0d,ctbadhst)
            if ctbadhst gt 0 then hstrdtmp[ibadhst] = !values.d_nan
            fwhm = initmaps.hst.smoothfwhm
            boxwidth = round((fwhm^2d)/2d + 1d)
            if not boxwidth then boxwidth++
;            hstrdsm = filter_image(hstrd,fwhm=initmaps.hst.smoothfwhm,/all)
            hstrdsm = filter_image(hstrdtmp,smooth=boxwidth,/iter,/all)
            ibadhst = where(finite(hstrdsm,/nan),ctbadhst)
            if ctbadhst gt 0 then hstrdsm[ibadhst] = bad
            hstrdtmp = 0
         endelse

;         rhst_fov_sm = ifsf_hstsubim(hstrdsm,[0,0],[dx,dy],$
;                                     initdat.platescale,$
;                                     initdat.positionangle,center_nuclei,$
;                                     hstrefcoords,$
;                                     initmaps.hstrdsm.scllim,$
;                                     sclargs=initmaps.hstrdsm.sclargs,$
;                                     /fov,hstps=hstpsrd,buffac=hstrd_buffac)
         rhst_fov_sm_ns= ifsf_hstsubim(hstrdsm,[0,0],[dx,dy],$
                                       initdat.platescale,$
                                       initdat.positionangle,center_nuclei,$
                                       hstrefcoords,[0,0],/noscl,$
                                       /fov,hstps=hstpsrd,buffac=hstrd_buffac)
         rhst_fov_sm_ns_rb = congrid(rhst_fov_sm_ns,dx,dy,/interp,/center)
         if tag_exist(initdat,'vormap') then begin
            tmp_rb = rhst_fov_sm_ns_rb
            nvor = max(initdat.vormap)
            badvor = where(~ finite(initdat.vormap),ctbad)
            if ctbad gt 0 then tmp_rb[badvor] = 0d
            for i=1,nvor do begin
               ivor = where(initdat.vormap eq i,ctbins)
               tmp_rb[ivor] = total(rhst_fov_sm_ns_rb[ivor])/ctbins
            endfor
            rhst_fov_sm_ns_rb = tmp_rb
         endif
         if tag_exist(initmaps.hstrdsm,'bgsub') then $
            rhst_fov_sm_ns_rb -= initmaps.hstrdsm.bgsub
      endif
   endif
   if dohstbl OR dohstrd then dohst=1b
   if dohst AND tag_exist(initmaps,'hstcol') then begin
      dohstcol=1
      if initmaps.hstbl.platescale ne initmaps.hstrd.platescale then begin
         print,'IFSF_MAKEMAPS: EROR: HST blue and red plate scales differ.'
         stop
      endif
;     HST ABmag computations. See
;     http://hla.stsci.edu/hla_faq.html#Source11
;     http://www.stsci.edu/hst/acs/documents/handbooks/currentDHB/acs_Ch52.html#102632
;     http://www.stsci.edu/hst/acs/analysis/zeropoints
;     Convert from e-/s to ABmags in the filters using:
;        ABmag = -2.5 log(e-/s) + ZP
;        ZP(ABmag) = -2.5 log(PHOTFLAM) - 2.408 - 5 log(PHOTPLAM)
;     The result will be mags per arcsec^2.       
      abortrd = 'No PHOTFLAM/PLAM keyword in red continuum image.'
      abortbl = 'No PHOTFLAM/PLAM keyword in blue continuum image.'
      if tag_exist(initmaps.hstbl,'photplam') then $
         pivotbl=initmaps.hstbl.photplam $
      else $
         pivotbl = sxpar(hstblhead,'PHOTPLAM',abortbl)
      if tag_exist(initmaps.hstrd,'photplam') then $
         pivotrd=initmaps.hstrd.photplam $
      else $
         pivotrd = sxpar(hstrdhead,'PHOTPLAM',abortrd)
      if tag_exist(initmaps.hstbl,'zp') then zpbl=initmaps.hstbl.zp $
      else $
         zpbl = -2.5d*alog10(sxpar(hstblhead,'PHOTFLAM',abortbl)) - 2.408d - $
                5d*alog10(sxpar(hstblhead,'PHOTPLAM',abortbl))
      if tag_exist(initmaps.hstrd,'zp') then zprd=initmaps.hstrd.zp $
      else $
         zprd = -2.5d*alog10(sxpar(hstrdhead,'PHOTFLAM',abortrd)) - 2.408d - $
                5d*alog10(sxpar(hstrdhead,'PHOTPLAM',abortrd))
;     Shift one image w.r.t. the other if necessary
;     Use red image as reference by default
      if tag_exist(initmaps.hstbl,'refcoords') AND $
         tag_exist(initmaps.hstrd,'refcoords') then begin
         if initmaps.hstbl.refcoords[0] ne initmaps.hstrd.refcoords[0] OR $
            initmaps.hstbl.refcoords[1] ne initmaps.hstrd.refcoords[1] then begin
               idiff = initmaps.hstrd.refcoords - initmaps.hstbl.refcoords
               hstbl = shift(hstbl,fix(idiff[0]),fix(idiff[1]))
         endif else begin
            hstrefcoords = initmaps.hstrd.refcoords
         endelse
      endif else begin
         hstrefcoords = initmaps.hst.refcoords
      endelse
;     Resize images to same size if necessary
      sizebl = size(hstbl)
      sizerd = size(hstrd)
      if sizebl[1] ne sizerd[1] OR sizebl[2] ne sizerd[2] then begin
         if sizebl[1] ne sizerd[1] then begin
            if sizebl[1] gt sizerd[1] then hstbl = hstbl[0:sizerd[1]-1,*] $
            else hstrd = hstrd[0:sizebl[1]-1,*]
         endif
         if sizebl[2] ne sizerd[2] then begin
            if sizebl[2] gt sizerd[2] then hstbl = hstbl[*,0:sizerd[2]-1] $
            else hstrd = hstrd[*,0:sizebl[2]-1]
         endif
     endif
;     Take a bunch of random samples of HST image
      if tag_exist(initmaps.hst,'sdevuplim') then uplim=initmaps.hst.sdevuplim $
      else uplim = 0.1 ; this gets rid of cosmic rays and stars ...
      if tag_exist(initmaps.hst,'sdevreg') then begin
         sdevreg = initmaps.hst.sdevreg
      endif else begin
         size_hst = size(hstrd)
         pxhst = round(size_hst[1]/10)
         pyhst = round(size_hst[2]/10)
         sdevreg[*,0] = [3*pxhst,4*pxhst,3*pyhst,4*pyhst]
         sdevreg[*,1] = [3*pxhst,4*pxhst,6*pyhst,7*pyhst]
         sdevreg[*,2] = [6*pxhst,7*pxhst,3*pyhst,4*pyhst]
         sdevreg[*,3] = [6*pxhst,7*pxhst,6*pyhst,7*pyhst]
      endelse
      size_sdevreg = size(sdevreg)
      nsdevreg = size_sdevreg[2]
      if size_sdevreg[0] eq 1 then begin
         sdevreg = reform(sdevreg,4,1)
         nsdevreg = 1
      endif
      sdev = dblarr(nsdevreg)
      for i=0,nsdevreg-1 do begin
         hsttmp = hstblsm[sdevreg[0,i]:sdevreg[1,i],sdevreg[2,i]:sdevreg[3,i]]
         sdev[i] = stddev(hsttmp[where(hsttmp ne 0 AND hsttmp le uplim)])
      endfor
      sdevbl = median(sdev)
      for i=0,nsdevreg-1 do begin
         hsttmp = hstrdsm[sdevreg[0,i]:sdevreg[1,i],sdevreg[2,i]:sdevreg[3,i]]
         sdev[i] = stddev(hsttmp[where(hsttmp ne 0 AND hsttmp le uplim)])
      endfor
      sdevrd = median(sdev)
;     Find bad pixels
      if tag_exist(initmaps.hstcol,'sigthr') then colsigthr = initmaps.hstcol.sigthr $
      else colsigthr = 3d
      ibdcol = where(hstrd le colsigthr*sdevrd OR $
                     hstbl le colsigthr*sdevbl,ctbdcol)
      hstcol = -2.5d*alog10(hstbl/hstrd) + zpbl - zprd
;     Galactic extinction correction:
;       (x - y)_intrinsic = (x - y)_obs - (A_x - A_y)
;       galextcor = A_x - A_y
      if tag_exist(initmaps.hstcol,'galextcor') then $
         hstcol -= initmaps.hstcol.galextcor
      hstcol[ibdcol] = bad
      hstcol_buffac = hstbl_buffac ge hstrd_buffac ? hstrd_buffac : hstbl_buffac
;     Extract and scale
      if tag_exist(initmaps.hst,'subim_sm') then begin
         chst_sm = ifsf_hstsubim(hstcol,[initmaps.hst.subim_sm,$
                                         initmaps.hst.subim_sm],$
                                 [dx,dy],initdat.platescale,$
                                 initdat.positionangle,center_nuclei,$
                                 hstrefcoords,$
                                 initmaps.hstcol.scllim,$
                                 sclargs=initmaps.hstcol.sclargs,hstps=hstpsbl)
      endif
      chst_big = ifsf_hstsubim(hstcol,[initmaps.hst.subim_big,$
                               initmaps.hst.subim_big],$
                               [dx,dy],initdat.platescale,$
                               initdat.positionangle,center_nuclei,$
                               hstrefcoords,$
                               initmaps.hstcol.scllim,$
                               sclargs=initmaps.hstcol.sclargs,hstps=hstpsbl)
      chst_fov = ifsf_hstsubim(hstcol,[0,0],[dx,dy],initdat.platescale,$
                               initdat.positionangle,center_nuclei,$
                               hstrefcoords,$
                               initmaps.hstcol.scllim,$
                               sclargs=initmaps.hstcol.sclargs,$
                               /fov,hstps=hstpsbl,buffac=hstcol_buffac)
      if ctbdcol gt 0 then begin
         hstcol_bad = hstcol*0d
         hstcol_bad[ibdcol] = 1d
         chst_fov_bad = ifsf_hstsubim(hstcol_bad,[0,0],[dx,dy],initdat.platescale,$
                                      initdat.positionangle,center_nuclei,$
                                      hstrefcoords,$
                                      initmaps.hstcol.scllim,/noscl,$
                                      /fov,hstps=hstpsbl,/badmask,buffac=hstcol_buffac)
         ibdcol_chst_fov = where(chst_fov_bad eq 1d,ctbdcol_chst_fov)
         if ctbdcol_chst_fov gt 0 then chst_fov[ibdcol_chst_fov] = 255b
         hstcol_bad = 0
      endif
;     Extract unscaled color image
      chst_fov_ns = ifsf_hstsubim(hstcol,[0,0],[dx,dy],initdat.platescale,$
                                  initdat.positionangle,center_nuclei,$
                                  hstrefcoords,$
                                  initmaps.hstcol.scllim,/noscl,/fov,$
                                  hstps=hstpsbl,buffac=hstcol_buffac)
   endif
   if tag_exist(initmaps,'hst') AND tag_exist(initmaps,'hstcolsm') then begin
      dohstcolsm=1
;     Shift one image w.r.t. the other if necessary
;     Use red image as reference by default
      if tag_exist(initmaps.hstbl,'refcoords') AND $
         tag_exist(initmaps.hstrd,'refcoords') then begin
         if initmaps.hstbl.refcoords[0] ne initmaps.hstrd.refcoords[0] AND $
            initmaps.hstbl.refcoords[1] ne initmaps.hstrd.refcoords[1] then begin
               idiff = initmaps.hstrd.refcoords - initmaps.hstbl.refcoords
               hstblsm = shift(hstblsm,fix(idiff[0]),fix(idiff[1]))
         endif else begin
            hstrefcoords = initmaps.hstrd.refcoords
         endelse
      endif else begin
         hstrefcoords = initmaps.hst.refcoords
      endelse
;     Resize images to same size if necessary
      sizebl = size(hstblsm)
      sizerd = size(hstrdsm)
      if sizebl[1] ne sizerd[1] OR sizebl[2] ne sizerd[2] then begin
         if sizebl[1] ne sizerd[1] then begin
            if sizebl[1] gt sizerd[1] then hstblsm = hstblsm[0:sizerd[1]-1,*] $
            else hstrdsm = hstrdsm[0:sizebl[1]-1,*]
         endif
         if sizebl[2] ne sizerd[2] then begin
            if sizebl[2] gt sizerd[2] then hstblsm = hstblsm[*,0:sizerd[2]-1] $
            else hstrdsm = hstrdsm[*,0:sizebl[2]-1]
         endif
      endif
;     Take a bunch of random samples of HST image
      if tag_exist(initmaps.hst,'sdevuplim') then uplim=initmaps.hst.sdevuplim $
      else uplim = 0.1 ; this gets rid of cosmic rays and stars ...
      if tag_exist(initmaps.hst,'sdevreg') then begin
         sdevreg = initmaps.hst.sdevreg
      endif else begin
         size_hst = size(hstrd)
         pxhst = round(size_hst[1]/10)
         pyhst = round(size_hst[2]/10)
         sdevreg[*,0] = [3*pxhst,4*pxhst,3*pyhst,4*pyhst]
         sdevreg[*,1] = [3*pxhst,4*pxhst,6*pyhst,7*pyhst]
         sdevreg[*,2] = [6*pxhst,7*pxhst,3*pyhst,4*pyhst]
         sdevreg[*,3] = [6*pxhst,7*pxhst,6*pyhst,7*pyhst]
      endelse
      size_sdevreg = size(sdevreg)
      nsdevreg = size_sdevreg[2]
      if size_sdevreg[0] eq 1 then begin
         sdevreg = reform(sdevreg,4,1)
         nsdevreg = 1
      endif
      sdev = dblarr(nsdevreg)
      for i=0,nsdevreg-1 do begin
         hsttmp = hstblsm[sdevreg[0,i]:sdevreg[1,i],sdevreg[2,i]:sdevreg[3,i]]
         sdev[i] = stddev(hsttmp[where(hsttmp ne 0 AND hsttmp le uplim)])
      endfor
      sdevbl = median(sdev)
      for i=0,nsdevreg-1 do begin
         hsttmp = hstrdsm[sdevreg[0,i]:sdevreg[1,i],sdevreg[2,i]:sdevreg[3,i]]
         sdev[i] = stddev(hsttmp[where(hsttmp ne 0 AND hsttmp le uplim)])
      endfor
      sdevrd = median(sdev)
;     Find bad pixels
      if tag_exist(initmaps.hst,'colsigthr') then colsigthr = initmaps.hst.colsigthr $
      else colsigthr = 3d
;     Color maps
      hstcolsm = -2.5d*alog10(hstblsm/hstrdsm) + zpbl - zprd
      if tag_exist(initmaps.hstcol,'galextcor') then $
         hstcolsm -= initmaps.hstcol.galextcor
      ibdcol = where(hstrdsm le colsigthr*sdevrd OR $
                     hstblsm le colsigthr*sdevbl)
      hstcolsm[ibdcol] = bad
;     Extract and scale
      if tag_exist(initmaps.hst,'subim_sm') then begin
         cshst_sm = ifsf_hstsubim(hstcolsm,[initmaps.hst.subim_sm,$
                                            initmaps.hst.subim_sm],$
                                  [dx,dy],initdat.platescale,$
                                  initdat.positionangle,center_nuclei,$
                                  hstrefcoords,$
                                  initmaps.hstcolsm.scllim,$
                                  sclargs=initmaps.hstcolsm.sclargs,hstps=hstpsbl)
      endif
      cshst_big = ifsf_hstsubim(hstcolsm,[initmaps.hst.subim_big,$
                                initmaps.hst.subim_big],$
                                [dx,dy],initdat.platescale,$
                                initdat.positionangle,center_nuclei,$
                                hstrefcoords,$
                                initmaps.hstcolsm.scllim,$
                                sclargs=initmaps.hstcolsm.sclargs,hstps=hstpsbl)
      cshst_fov_s = ifsf_hstsubim(hstcolsm,[0,0],[dx,dy],initdat.platescale,$
                                initdat.positionangle,center_nuclei,$
                                hstrefcoords,$
                                initmaps.hstcolsm.scllim,$
                                sclargs=initmaps.hstcolsm.sclargs,$
                                /fov,hstps=hstpsbl,buffac=hstcol_buffac)
;     Extract unscaled color image and convert to same pixel scale as IFS data
      cshst_fov_ns = ifsf_hstsubim(hstcolsm,[0,0],[dx,dy],initdat.platescale,$
                                   initdat.positionangle,center_nuclei,$
                                   hstrefcoords,[0,0],/noscl,/fov,$
                                   hstps=hstpsbl,buffac=hstcol_buffac)
      cshst_fov_rb = congrid(cshst_fov_ns,dx,dy,/interp,/center)
   endif
   if dohst AND dohstcol AND tag_exist(initdat,'vormap') then begin
      dohstcolvor=1
      cshst_fov_rb = -2.5d*alog10(bhst_fov_sm_ns_rb/rhst_fov_sm_ns_rb) + zpbl - zprd
      inan = where(finite(cshst_fov_rb,/nan),ctnan)
      if tag_exist(initmaps.hstcol,'galextcor') then $
         cshst_fov_rb -= initmaps.hstcol.galextcor
      if ctnan gt 0 then cshst_fov_rb[inan] = bad
   endif
   hstrd=0
   hstbl=0
   hstcol=0
   hstrdsm=0
   hstblsm=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Fit QSO PSF
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  Radii in kpc
;  GMOS FOV
   map_x = rebin(dindgen(dx)+1,dx,dy)
   map_y = rebin(transpose(dindgen(dy)+1),dx,dy)
   map_r = sqrt((map_x - center_axes[0])^2d + (map_y - center_axes[1])^2d)
   map_rkpc_ifs = map_r * kpc_per_pix

;  PA E of N, in degrees; spaxel with [0,0] has PA = bad
   map_pa = dblarr(dx,dy)
   map_xaxis = map_x - center_axes[0]
   map_yaxis = map_y - center_axes[1]
   map_pa = ifsf_pa(map_xaxis,map_yaxis)
   map_pa += initdat.positionangle
   iphase = where(map_pa ge 360d,ctphase)
   if ctphase gt 0 then map_pa[iphase] -= 360d
   inuc = where(map_r eq 0d,ctnuc)
   if ctnuc gt 0 then map_pa[inuc] = bad

   if tag_exist(initdat,'decompose_qso_fit') then begin
      inan = where(finite(contcube.qso_mod,/nan),ctnan)
      if ctnan gt 0 then contcube.qso_mod[inan] = 0d
      qso_map = total(contcube.qso_mod,3) / contcube.npts
      maxqso_map = max(qso_map)
;      qso_err = stddev(contcube.qso,dim=3,/double)
;      qso_err = sqrt(total(datacube.var,3))
      qso_err = sqrt(median(datacube.var,dim=3,/double))
      ibd = where(~ finite(qso_err),ctbd)
      if ctbd gt 0 then begin
         qso_err[ibd] = bad
         qso_map[ibd] = 0d
      endif
      qso_map /= maxqso_map
      qso_err /= maxqso_map
;      izero = where(qso_map eq 0d,ctzero)
;      maxerr = max(qso_err)
;      lowthresh=1d-4
;      if ctzero gt 0 then begin
;         qso_map[izero] = lowthresh
;         qso_err[izero] = maxerr*100d
;      endif
;      ilow = where(qso_map lt lowthresh,ctlow)
;      if ctlow gt 0 then begin
;         qso_map[ilow] = lowthresh
;         qso_err[ilow] = maxerr*100d
;      endif
      
;      qso_err /= max(median(datacube.dat,dim=3,/double))
      
;     2D Moffat fit to continuum flux vs. radius
      parinit = REPLICATE({value:0d, fixed:0b, limited:[0B,0B], tied:'', $
                           limits:[0d,0d]},8)
      est=[0d,1d,1d,1d,center_nuclei[0]-1d,center_nuclei[1]-1d,0d,3d]
      parinit.value = est
      parinit[7].limited = [1b,1b]
      parinit[7].limits = [0d,5d]
      qso_fit = mpfit2dpeak(qso_map,qso_fitpar,/circular,/moffat,est=est,$
                            error=qso_err,parinfo=parinit)
      map_rnuc = sqrt((map_x - qso_fitpar[4]+1)^2d + $
                      (map_y - qso_fitpar[5]+1)^2d)
      map_rnuckpc_ifs = map_rnuc * kpc_per_pix
      psf1d_x = dindgen(101)/100d*max(map_rnuckpc_ifs)
      psf1d_y = alog10(moffat(psf1d_x,[qso_fitpar[1],0d,$
                              qso_fitpar[2]*kpc_per_pix,$
                              qso_fitpar[7]]))

   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Compute radii and centers
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  coordinates in kpc
   xran_kpc = double([-(center_axes[0]-0.5),dx-(center_axes[0]-0.5)]) $
              * kpc_per_pix
   yran_kpc = double([-(center_axes[1]-0.5),dy-(center_axes[1]-0.5)]) $
              * kpc_per_pix
   center_nuclei_kpc_x = (center_nuclei[0,*] - center_axes[0]) $
                         * kpc_per_pix
   center_nuclei_kpc_y = (center_nuclei[1,*] - center_axes[1]) $
                         * kpc_per_pix
;  in image window
   xwinran_kpc = double([-(center_axes[0]-0.5-(plotwin[0]-1)),$
                         dxwin-(center_axes[0]-0.5-(plotwin[0]-1))]) $
                        * kpc_per_pix
   ywinran_kpc = double([-(center_axes[1]-0.5-(plotwin[1]-1)),$
                         dywin-(center_axes[1]-0.5-(plotwin[1]-1))]) $
                        * kpc_per_pix
;   center_nuclei_kpc_xwin = (center_nuclei[0,*]-center_axes[0]-(plotwin[0]-1)) $
;                            * kpc_per_pix
;   center_nuclei_kpc_ywin = (center_nuclei[1,*]-center_axes[1]-(plotwin[1]-1)) $
;                            * kpc_per_pix
   center_nuclei_kpc_xwin = center_nuclei_kpc_x
   center_nuclei_kpc_ywin = center_nuclei_kpc_y
   
;  HST FOV
   if (dohstrd OR dohstbl) then begin
      if dohstbl then size_subim = size(bhst_fov) $
      else size_subim = size(rhst_fov_sc)
      map_x_hst = rebin(dindgen(size_subim[1])+1,size_subim[1],size_subim[2])
      map_y_hst = rebin(transpose(dindgen(size_subim[2])+1),$
                        size_subim[1],size_subim[2])
;     Locations of [0,0] point on axes and nuclei, in HST pixels (single-offset
;     indices).
      center_axes_hst = (center_axes-0.5d) * double(size_subim[1]/dx) + 0.5d
      center_nuclei_hst = (center_nuclei-0.5d) * double(size_subim[1]/dx) + 0.5d
;     Radius of each HST pixel from axis [0,0] point, in HST pixels
      map_r_hst = sqrt((map_x_hst - center_axes_hst[0])^2d + $
                       (map_y_hst - center_axes_hst[1])^2d)
      if dohstbl then begin
;        Radius of each HST pixel from axis [0,0] point, in kpc
         map_rkpc_hst = map_r_hst * initmaps.hstbl.platescale * kpc_per_as
         hstplatescale = initmaps.hstbl.platescale
         kpc_per_hstpix = hstplatescale * kpc_per_as
         if dohstrd then begin
            if initmaps.hstbl.platescale ne initmaps.hstrd.platescale then begin
               print,'WARNING: HST blue and red plate scales differ;'
               print,'         using blue platescale for radius calculations.'
            endif
         endif
      endif else begin
         map_rkpc_hst = map_r_hst * initmaps.hstrd.platescale * kpc_per_as      
         hstplatescale = initmaps.hstrd.platescale
         kpc_per_hstpix = hstplatescale * kpc_per_as
      endelse
      if dohstbl then $
         if tag_exist(initmaps.hstbl,'nucoffset') then begin
;           Radius of each blue HST pixel from axis [0,0] point, in HST pixels, 
;           with by-hand offset applied
            map_r_bhst = $
               sqrt((map_x_hst - $
                     (center_axes_hst[0]+initmaps.hstbl.nucoffset[0]))^2d + $
                    (map_y_hst - $
                     (center_axes_hst[1]+initmaps.hstbl.nucoffset[1]))^2d)
;           ... and now in kpc
            map_rkpc_bhst = map_r_bhst * initmaps.hstbl.platescale * kpc_per_as
         endif else map_rkpc_bhst = map_rkpc_hst
      if dohstrd then $
         if tag_exist(initmaps.hstrd,'nucoffset') then begin
            map_r_rhst = $
               sqrt((map_x_hst - $
                     (center_axes_hst[0]+initmaps.hstrd.nucoffset[0]))^2d + $
                    (map_y_hst - $
                     (center_axes_hst[1]+initmaps.hstrd.nucoffset[1]))^2d)
            map_rkpc_rhst = map_r_rhst * initmaps.hstrd.platescale * kpc_per_as
         endif else map_rkpc_rhst = map_rkpc_hst 
   endif else map_rkpc_hst = 0d

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Process emission lines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   emlvel = 0d
   emlflxcor_pp = 0d
   emlflxcor_med = 0d
   ebv = 0d
   ebvmed = 0d
   errebv = 0d
   lr = hash()
   lrerrlo = hash()
   lrerrhi = hash()
   elecdenmap = hash()
   elecdenmap_errlo = hash()
   elecdenmap_errhi = hash()
   if ~ tag_exist(initdat,'noemlinfit') then begin

;;     Sort emission line components if requested
;      if tag_exist(initmaps,'fcnsortcomp') AND $
;         tag_exist(initmaps,'sortlines') AND $
;         tag_exist(initmaps,'sorttype') then begin
;         if tag_exist(initmaps,'argssortcomp') then $
;            linmaps = call_function(initmaps.fcnsortcomp,dx,dy,linmaps,$
;                                    initdat.linetie,initmaps.sortlines,$
;                                    initmaps.sorttype,$
;                                    _extra=initmaps.argssortcomp) $
;         else $
;            linmaps = call_function(initmaps.fcnsortcomp,dx,dy,linmaps,$
;                                    initdat.linetie,initmaps.sortlines,$
;                                    initmaps.sorttype)
;      endif

;;     Set reference line for rotation curve, for defining outflows, by resorting
;;     OUTLINES so that reference line comes first.
;      if ~ tag_exist(initmaps,'diskline') then diskline='Halpha' $
;      else diskline = initmaps.diskline
;      idiskline = outlines.where(diskline,count=ctdiskline)
;      if ctdiskline gt 0 then outlines.move,idiskline,0 $
;      else message,$
;         'Disk line not found; using first line in OUTLINES hash.',/cont
                                  
;      ofpars_line=0b
;      sigthresh=0b
;      ofthresh=0b ; threshold for defining outflow
;      ofignore=0b
;      diffthresh=0b
;      if tag_exist(initmaps,'compof') then begin
;         ofpars_line=1b
;         if tag_exist(initmaps,'compof_sigthresh') then $
;            sigthresh=initmaps.compof_sigthresh
;         if tag_exist(initmaps,'compof_ofthresh') then $
;            ofthresh=initmaps.compof_ofthresh
;         if tag_exist(initmaps,'compof_diffthresh') then $
;            diffthresh=initmaps.compof_diffthresh
;         if tag_exist(initmaps,'compof_ignore') then $
;            ofignore=initmaps.compof_ignore
;      endif
;      if line eq diskline then diskrot=0b $
;      else diskrot=linspecpars[diskline].vpk
;      if tag_exist(initmaps,'compof') then ofpars[line] = ofpars_line
;      endforeach
;      linspecpars_tags = tag_names(linspecpars[outlines[0]])
;      if tag_exist(initmaps,'compof') then $
;         ofpars_tags = tag_names(ofpars[outlines[0]])

;     Compute CVDF velocities and fluxes
      if tag_exist(initmaps,'cvdf') then begin
         if tag_exist(initmaps.cvdf,'flux_maps') then $
            fluxvels=initmaps.cvdf.flux_maps $
         else fluxvels=0b
         if tag_exist(initmaps.cvdf,'sigcut') then $
            sigcut=initmaps.cvdf.sigcut $
         else sigcut=0b
      endif
      emlcompvel = ifsf_cmpcompvals(emlwav,emlsig,emlflx,initdat.zsys_gas,$
                                    emlwaverr=emlwaverr,emlsigerr=emlsigerr)
      emlvel = emlcompvel
      if not tag_exist(initdat,'nocvdf') then begin
         emlcvdfvel = ifsf_cmpcvdfvals(emlcvdf,emlflx,emlflxerr,$
                                       fluxvels=fluxvels,sigcut=sigcut)
         emlvel = emlcvdfvel + emlcompvel
      endif
        
;     Extinction maps
;     flux summed over components
      if tag_exist(initmaps,'ebv') then begin
;        Calculate ...
         if tag_exist(initmaps.ebv,'calc') then begin
            ebv = hash()
            errebv = hash()
            ebvmed = hash()
            foreach key,initmaps.ebv.calc do begin
               if tag_exist(initmaps,'argslineratios') then $
                  ebvtmp = $
                     ifsf_lineratios(emlflx[key],emlflxerr[key],linelist,/ebvonly,$
                                     errlo=errtmp,_extra=initmaps.argslineratios) $
               else ebvtmp = $
                  ifsf_lineratios(emlflx[key],emlflxerr[key],linelist,/ebvonly,$
                                  errlo=errtmp)
               ebv[key] = ebvtmp['ebv']
               errebv[key] = errtmp['ebv']
               igdebv = where(ebv[key] ne bad,ctgdebv)
               if ctgdebv gt 0 then ebvmed[key] = median(ebv[key,igdebv])
            endforeach
;           ... and apply.
            if tag_exist(initmaps.ebv,'apply') then begin
               ebvkey = 'ftot'
               emlflxcor_pp = hash()
               emlflxcor_med = hash()
               emlflxerrcor_pp = hash()
               emlflxerrcor_med = hash()
               emlflxcor_pp['ftot'] = hash()
               emlflxcor_med['ftot'] = hash()
               emlflxerrcor_pp['ftot'] = hash()
               emlflxerrcor_med['ftot'] = hash()
               for icomp=1,initdat.maxncomp do begin
                  stric = string(icomp,format='(I0)')
                  emlflxcor_pp['fc'+stric] = hash()
                  emlflxcor_med['fc'+stric] = hash()
                  emlflxerrcor_pp['fc'+stric] = hash()
                  emlflxerrcor_med['fc'+stric] = hash()
                  foreach line,linelist.keys() do begin
                     emlflxcor_pp['fc'+stric,line] = emlflx['fc'+stric,line]
                     emlflxcor_med['fc'+stric,line] = emlflx['fc'+stric,line]
                     emlflxerrcor_pp['fc'+stric,line] = emlflx['fc'+stric,line]
                     emlflxerrcor_med['fc'+stric,line] = emlflx['fc'+stric,line]
                     flx = emlflx['fc'+stric,line]
                     flxerr = emlflxerr['fc'+stric,line]
                     igdflx = where(flx gt 0 AND flx ne bad,ctgdflx)
                     igdebv = where(ebv[ebvkey] ne bad,ctgdebv)
                     if ctgdflx gt 0 AND ctgdebv gt 0 then begin
                        ebvuse = ebv[ebvkey]
                        ibadebv = cgsetdifference(igdflx,igdebv,count=ctbadebv)
                        if ctbadebv gt 0 then $
                           ebvuse[ibadebv] = ebvmed[ebvkey]
                        emlflxcor_pp['fc'+stric,line,igdflx] = $
                           ifsf_dustcor_ccm(linelist[line],flx[igdflx],$
                                            ebvuse[igdflx])
                        emlflxcor_med['fc'+stric,line,igdflx] = $
                           ifsf_dustcor_ccm(linelist[line],flx[igdflx],$
                                            ebvmed[ebvkey])
; Should really propagate E(B-V) error as well ...
                        emlflxerrcor_pp['fc'+stric,line,igdflx] = $
                           ifsf_dustcor_ccm(linelist[line],flxerr[igdflx],$
                                            ebvuse[igdflx])
                        emlflxerrcor_med['fc'+stric,line,igdflx] = $
                           ifsf_dustcor_ccm(linelist[line],flxerr[igdflx],$
                                            ebvmed[ebvkey])
                     endif
                  endforeach
               endfor
            endif

         endif
      endif else begin
         ebv = 0b
         errebv = 0b
      endelse

;     Line ratios
      if tag_exist(initmaps,'lr') then begin
         if tag_exist(initmaps.lr,'calc') then begin
            foreach key,initmaps.lr.calc do begin
               lrtmp = $
                  ifsf_lineratios(emlflx[key],emlflxerr[key],linelist,/lronly,$
                                  errlo=errlotmp,errhi=errhitmp)
               lr[key] = hash()
               lrerrlo[key] = hash()
               lrerrhi[key] = hash()
               foreach lrloop,lrtmp.keys() do begin 
                  lr[key,lrloop] = lrtmp[lrloop]
                  lrerrlo[key,lrloop] = errlotmp[lrloop]
                  lrerrhi[key,lrloop] = errhitmp[lrloop]
               endforeach


;              Compute electron densities
;              Use #s from Sanders, Shapley, et al. 2015
               if lrtmp.haskey('s2') then begin
;                 Densities
                  tmps2map = lrtmp['s2']
                  igds2 = where(tmps2map ne bad AND finite(tmps2map),ctgd)
;                  ibds2 = where(tmps2map eq bad OR ~ finite(tmps2map),ctgd)
                  tmps2mapgd = 10d^tmps2map[igds2]
                  tmpdenmap = dblarr(dx,dy)+bad
                  igdden = where(tmps2mapgd gt s2_minrat AND $
                                 tmps2mapgd lt s2_maxrat)
                  tmpdenmapgd = tmpdenmap[igds2]
                  tmpdenmapgd[igdden] = alog10((s2_c*tmps2mapgd[igdden] - s2_a*s2_b)/$
                                               (s2_a - tmps2mapgd[igdden]))
                  ilo = where(tmps2mapgd ge s2_maxrat OR $
                              tmpdenmapgd lt alog10(s2_minden),ctlo)
                  ihi = where(tmps2mapgd le s2_minrat OR $
                              (tmpdenmapgd gt alog10(s2_maxden) AND $
                               tmpdenmapgd ne bad),cthi)
                  if ctlo gt 0 then tmpdenmapgd[ilo] = alog10(s2_minden)
                  if cthi gt 0 then tmpdenmapgd[ihi] = alog10(s2_maxden)
                  tmpdenmap[igds2] = tmpdenmapgd

;                 Density upper limits
                  tmps2maperrlo = errlotmp['s2']
                  tmps2mapgderrlo = tmps2mapgd - $
                                    10d^(tmps2map[igds2]-tmps2maperrlo[igds2])
                  tmpdenmaperrhi = dblarr(dx,dy) + bad
                  tmpdenmaperrhigd = tmpdenmaperrhi[igds2]
                  tmpdenmaphi = dblarr(dx,dy) + bad
                  tmpdenmaphigd = tmpdenmaphi[igds2]
                  igdden = where(tmps2mapgd - tmps2mapgderrlo gt s2_minrat AND $
                                 tmps2mapgd - tmps2mapgderrlo lt s2_maxrat)
                  tmpdenmaphigd[igdden] = $
                     alog10((s2_c*(tmps2mapgd[igdden]-tmps2mapgderrlo[igdden]) - s2_a*s2_b)/$
                            (s2_a - (tmps2mapgd[igdden]-tmps2mapgderrlo[igdden])))
                  tmpdenmaperrhigd[igdden] = tmpdenmaphigd[igdden] - $
                                             tmpdenmapgd[igdden]
                  ilo = where(tmps2mapgd - tmps2mapgderrlo ge s2_maxrat OR $
                              tmpdenmaphigd le alog10(s2_minden),ctlo)
                  ihi = where((tmps2mapgd - tmps2mapgderrlo le s2_minrat AND $
                              tmps2mapgd - tmps2mapgderrlo gt 0d) OR $
                              (tmpdenmaphigd ge alog10(s2_maxden) AND $
                               tmpdenmaphigd ne bad),cthi)
                  if ctlo gt 0 then tmpdenmaperrhigd[ilo] = 0d
                  if cthi gt 0 then tmpdenmaperrhigd[ihi] = $
                     alog10(s2_maxden) - tmpdenmapgd[ihi]
                  tmpdenmaperrhi[igds2] = tmpdenmaperrhigd

;                 Density lower limits
                  tmps2maperrhi = errhitmp['s2']
                  tmps2mapgderrhi = 10d^(tmps2map[igds2]+tmps2maperrhi[igds2]) - $
                                         tmps2mapgd
                  tmpdenmaperrlo = dblarr(dx,dy)+bad
                  tmpdenmaperrlogd = tmpdenmaperrlo[igds2]
                  tmpdenmaplo = dblarr(dx,dy)+bad
                  tmpdenmaplogd = tmpdenmaplo[igds2]
                  igdden = where(tmps2mapgd + tmps2mapgderrhi gt s2_minrat AND $
                                 tmps2mapgd + tmps2mapgderrhi lt s2_maxrat)
                  tmpdenmaplogd[igdden] = $
                     alog10((s2_c*(tmps2mapgd[igdden]+tmps2mapgderrhi[igdden]) - s2_a*s2_b)/$
                            (s2_a - (tmps2mapgd[igdden]+tmps2mapgderrhi[igdden])))
                  tmpdenmaperrlogd[igdden] = tmpdenmapgd[igdden] - $
                                             tmpdenmaplogd[igdden]
                  ilo = where(tmps2mapgd + tmps2mapgderrhi ge s2_maxrat OR $
                              tmpdenmaplogd le alog10(s2_minden),ctlo)
                  ihi = where((tmps2mapgd + tmps2mapgderrhi le s2_minrat AND $
                              tmps2mapgd + tmps2mapgderrhi gt 0d) OR $
                              (tmpdenmaplogd ge alog10(s2_maxden) AND $
                               tmpdenmaplogd ne bad),cthi)
                  if ctlo gt 0 then tmpdenmaperrlogd[ilo] = $
                     tmpdenmapgd[ilo] - alog10(s2_minden)
                  if cthi gt 0 then tmpdenmaperrlogd[ihi] = 0d
                  tmpdenmaperrlo[igds2] = tmpdenmaperrlogd

                  elecdenmap[key] = tmpdenmap
                  elecdenmap_errlo[key] = tmpdenmaperrlo
                  elecdenmap_errhi[key] = tmpdenmaperrhi
               endif

            endforeach
         endif
      endif

   endif else begin
      emlflx = 0d
   endelse


   if tag_exist(initdat,'donad') then begin

;     Apply S/N cut. Funny logic is to avoid loops.
     map = nadfit.weqabs[*,*,0]
     maperravg = (nadfit.weqabserr[*,*,0] + nadfit.weqabserr[*,*,1]) / 2d
     mask1 = dblarr(dx,dy)
     mask2 = dblarr(dx,dy)
     igd = where(map gt 0d AND $
                 map ne bad AND $
                 map ge initmaps.nadabsweq_snrthresh*maperravg)
     ibd = where(map eq 0d OR $
                 map eq bad OR $
                 map lt initmaps.nadabsweq_snrthresh*maperravg)
     mask1[igd]=1d
     mask2[ibd]=bad
     rbmask1 = rebin(mask1,dx,dy,initnad.maxncomp)
     rbmask2 = rebin(mask2,dx,dy,initnad.maxncomp)
     nadfit.waveabs*=rbmask1
     nadfit.waveabserr*=rbmask1
     nadfit.sigmaabs*=rbmask1
     nadfit.sigmaabserr*=rbmask1
     nadfit.tau*=rbmask1
     nadfit.tauerr*=rbmask1
     nadfit.waveabs+=rbmask2
     nadfit.waveabserr+=rbmask2
     nadfit.sigmaabs+=rbmask2
     nadfit.sigmaabserr+=rbmask2
     nadfit.tau+=rbmask2
     nadfit.tauerr+=rbmask2


;     Compute velocities and column densities of NaD model fits
      nadabscftau = dblarr(dx,dy,maxnadabscomp+1)+bad
      errnadabscftau = dblarr(dx,dy,maxnadabscomp+1,2)+bad
      nadabstau = dblarr(dx,dy,maxnadabscomp+1)+bad
      errnadabstau = dblarr(dx,dy,maxnadabscomp+1,2)+bad
      nadabscf = dblarr(dx,dy,maxnadabscomp+1)+bad
      nadabsvel = dblarr(dx,dy,maxnadabscomp+1)+bad
      errnadabsvel = dblarr(dx,dy,maxnadabscomp+1,2)+bad
      nademvel = dblarr(dx,dy,maxnademcomp+1)+bad
      errnademvel = dblarr(dx,dy,maxnademcomp+1)+bad
      nadabssig = dblarr(dx,dy,maxnadabscomp+1)+bad
      errnadabssig = dblarr(dx,dy,maxnadabscomp+1,2)+bad
      nademsig = dblarr(dx,dy,maxnademcomp+1)+bad
      errnademsig = dblarr(dx,dy,maxnademcomp+1)+bad
      nadabsv98 = dblarr(dx,dy,maxnadabscomp+1)+bad
;      errnadabsv98 = dblarr(dx,dy,maxnadabscomp+1,2)+bad
      nademv98 = dblarr(dx,dy,maxnademcomp+1)+bad
;      errnademv98 = dblarr(dx,dy,maxnademcomp+1,2)+bad
      nadabsnh = dblarr(dx,dy)+bad
      llnadabsnh = bytarr(dx,dy)
      errnadabsnh = dblarr(dx,dy,2)+bad
      nadabslnnai = dblarr(dx,dy)+bad
      llnadabslnnai = bytarr(dx,dy)
      errnadabslnnai = dblarr(dx,dy,2)+bad
      nadabsnhcf = dblarr(dx,dy)+bad
      errnadabsnhcf = dblarr(dx,dy,2)+bad
      nadabscnhcf = dblarr(dx,dy,maxnadabscomp+1)+bad
      errnadabscnhcf = dblarr(dx,dy,maxnadabscomp+1,2)+bad
      nadabsncomp = intarr(dx,dy)+bad
      nademncomp = intarr(dx,dy)+bad
      for i=0,dx-1 do begin
         for j=0,dy-1 do begin
            igd = where(nadfit.waveabs[i,j,*] ne bad AND $
                        nadfit.waveabs[i,j,*] ne 0,ctgd)
            if ctgd gt 0 then begin
               tmpcftau=nadfit.cf[i,j,igd]*nadfit.tau[i,j,igd]
               tmpcf=nadfit.cf[i,j,igd]
               tmptau=nadfit.tau[i,j,igd]
               tmpcftauerr = nadfit.cferr[i,j,igd,*]
               tmpltauerrlo = nadfit.tauerr[i,j,igd,0] ; errors are in log(tau) space
               tmpltauerrhi = nadfit.tauerr[i,j,igd,1] ; errors are in log(tau) space
               tmpwaveabs = nadfit.waveabs[i,j,igd]
               tmpwaveabserr = nadfit.waveabserr[i,j,igd,*]
               tmpsigabs=nadfit.sigmaabs[i,j,igd]
               tmpsigabserr = nadfit.sigmaabserr[i,j,igd,*]
               ineg = where(tmpcftau[0,0,igd] gt 0d AND $
                            tmpcftau[0,0,igd] lt tmpcftauerr[0,0,igd,0],ctneg)
               if ctneg gt 0 then tmpcftau[0,0,igd,0] = 0d
               nnai = total(tmptau*nadfit.sigmaabs[i,j,igd]) / $
                            (1.497d-15/sqrt(2d)*nadlinelist['NaD1']*0.3180d)
               nadabslnnai[i,j] = alog10(nnai)
               nnaicf = total(tmpcftau*nadfit.sigmaabs[i,j,igd]*$
                        nadfit.cf[i,j,igd]) / $
                        (1.497d-15/sqrt(2d)*nadlinelist['NaD1']*0.3180d)
               tmptauerrlo = tmptau - 10d^(alog10(tmptau) - tmpltauerrlo)
               tmptauerrhi = 10d^(alog10(tmptau) + tmpltauerrhi) - tmptau
;              get saturated points
               isat = where(tmptau eq taumax,ctsat)
               if ctsat gt 0 then tmptauerrhi[isat] = 0d
               nnaierrlo = $
                  nnai*sqrt(total((tmptauerrlo/tmptau)^2d + $
                                  (nadfit.sigmaabserr[i,j,igd,0]/$
                                   nadfit.sigmaabs[i,j,igd])^2d))
               nnaierrhi = $
                  nnai*sqrt(total((tmptauerrhi/tmptau)^2d + $
                                  (nadfit.sigmaabserr[i,j,igd,1]/$
                                   nadfit.sigmaabs[i,j,igd])^2d))
               nnaicferr = $
                  [nnaicf,nnaicf]*sqrt(total((nadfit.cferr[i,j,igd,*]/$
                  rebin(tmpcftau,ctgd,2))^2d + $
                  (nadfit.sigmaabserr[i,j,igd,*]/$
                  rebin(nadfit.sigmaabs[i,j,igd],ctgd,2))^2d,3))
               nadabsnh[i,j] = nnai/(1-ionfrac)/10^(naabund - nadep)
               nadabsnhcf[i,j] = nnaicf/(1-ionfrac)/10^(naabund - nadep)
               errnadabsnh[i,j,*] = [nnaierrlo,nnaierrhi]/$
                                    (1-ionfrac)/10^(naabund - nadep)
               errnadabsnhcf[i,j,*] = $
                  nadabsnhcf[i,j]*sqrt((nnaicferr/nnaicf)^2d + $
                                       (oneminusionfrac_relerr)^2d)
               nadabsncomp[i,j] = ctgd
               if ctsat gt 0 then begin
                  llnadabsnh[i,j] = 1b
                  errnadabsnh[i,j,1] = 0d
                  llnadabslnnai[i,j] = 1b
                  errnadabslnnai[i,j,1] = 0d
               endif
            endif
;           Sort absorption line wavelengths. In output arrays,
;           first element of third dimension holds data for spaxels with only
;           1 component. Next elements hold velocities for spaxels with more than
;           1 comp, in order of increasing blueshift. Formula for computing error
;           in velocity results from computing derivative in Wolfram Alpha w.r.t.
;           lambda and rearranging on paper.
            if ctgd eq 1 then begin
;              Set low sigma error equal to best-fit minus 5 km/s 
;              (unlikely that it's actually lower than this!) 
;              if low is greater than best-fit value
               if tmpsigabs lt tmpsigabserr[0] then tmpsigabserr[0] = tmpsigabs-5d
               errnadabslnnai[i,j,0] = $
                  sqrt(tmpltauerrlo^2d + $
                       (alog10(tmpsigabs)-alog10(tmpsigabs-tmpsigabserr[0]))^2d)
               errnadabslnnai[i,j,1] = $
                  sqrt(tmpltauerrhi^2d + $
                       (alog10(tmpsigabs+tmpsigabserr[1])-alog10(tmpsigabs))^2d)
               nadabscnhcf[i,j,0]=nadabsnhcf[i,j]
               errnadabscnhcf[i,j,0,*]=errnadabsnhcf[i,j,*]
               zdiff = $
                  tmpwaveabs/(nadlinelist['NaD1']*(1d + initdat.zsys_gas)) - 1d
               nadabsvel[i,j,0] = c_kms * ((zdiff+1d)^2d - 1d) / $
                                          ((zdiff+1d)^2d + 1d)
               errnadabsvel[i,j,0,*] = $
                  c_kms * (4d/(nadlinelist['NaD1']*(1d + initdat.zsys_gas))*$
                           ([zdiff,zdiff]+1d)/(([zdiff,zdiff]+1d)^2d + 1d)^2d) * $
                           tmpwaveabserr
               nadabssig[i,j,0] = tmpsigabs
               errnadabssig[i,j,0,*] = tmpsigabserr
               nadabsv98[i,j,0] = nadabsvel[i,j,0]-2d*nadabssig[i,j,0]
;               errnadabsv98[i,j,0] = $
;                  sqrt(errnadabsvel[i,j,0]^2d + 4d*errnadabssig[i,j,0]^2d)
               nadabscftau[i,j,0] = tmpcftau
               errnadabscftau[i,j,0,*] = tmpcftauerr
               nadabstau[i,j,0] = tmptau
               nadabscf[i,j,0] = tmpcf
               errnadabstau[i,j,0,*] = [tmptauerrlo,tmptauerrhi]
            endif else if ctgd gt 1 then begin
               errnadabslnnai[i,j,0] = 0d
               errnadabslnnai[i,j,1] = 0d
               for k=0,ctgd-1 do begin
                  if tmpsigabs[0,0,k] lt tmpsigabserr[0,0,k,0] then $
                     tmpsigabserr[0,0,k,0] = tmpsigabs[0,0,k] - 5d
                  errnadabslnnai[i,j,0] += $
                     tmpltauerrlo[0,0,k]^2d + $
                     (alog10(tmpsigabs[0,0,k])-$
                      alog10(tmpsigabs[0,0,k]-tmpsigabserr[0,0,k,0]))^2d
                  errnadabslnnai[i,j,1] += $
                     tmpltauerrhi[0,0,k]^2d + $
                     (alog10(tmpsigabs[0,0,k]+tmpsigabserr[0,0,k,1])-$
                      alog10(tmpsigabs[0,0,k]))^2d
               endfor
               errnadabslnnai[i,j,*] = sqrt(errnadabslnnai[i,j,*])
               sortgd = sort(tmpwaveabs)
;               nnai = reverse(tmptau[sortgd]*tmpsigabs[sortgd]) / $
;                      (1.497d-15/sqrt(2d)*nadlinelist['NaD1']*0.3180d)
               nnaicf = reverse(tmpcftau[sortgd]*tmpsigabs[sortgd]) / $
                        (1.497d-15/sqrt(2d)*nadlinelist['NaD1']*0.3180d)
;               nadabscnh[i,j,1:ctgd] = nnai/(1-ionfrac)/10^(naabund - nadep)
               nadabscnhcf[i,j,1:ctgd] = nnaicf/(1-ionfrac)/10^(naabund - nadep)
;               nnaierrlo = $
;                  nnai*sqrt(reverse((tmptauerrlo[sortgd]/tmptau[sortgd])^2d + $
;                                    (tmpsigabserr[sortgd]/tmpsigabs[sortgd])^2d))
;               nnaierrhi = $
;                  nnai*sqrt(reverse((tmptauerrhi[sortgd]/tmptau[sortgd])^2d + $
;                                    (tmpsigabserr[sortgd]/tmpsigabs[sortgd])^2d))
               nnaicferr = $
                  rebin(nnaicf,1,1,ctgd,2)*$
                  sqrt(reverse((tmpcftauerr[0,0,sortgd,*]/$
                  rebin(tmpcftau[sortgd],1,1,ctgd,2))^2d + $
                  (tmpsigabserr[0,0,sortgd,*]/$
                  rebin(tmpsigabs[sortgd],1,1,ctgd,2))^2d,3))
;               nnaicferrhi = $
;                  nnaicf*sqrt(reverse((tmpcftauerrhi[sortgd]/tmpcftau[sortgd])^2d + $
;                                    (tmpsigabserr[sortgd]/tmpsigabs[sortgd])^2d))
;               errnadabscnh[i,j,1:ctgd,0] = nnaierrlo / $
;                                            (1-ionfrac)/10^(naabund - nadep)
;               errnadabscnh[i,j,1:ctgd,1] = nnaierrhi / $
;                                            (1-ionfrac)/10^(naabund - nadep)
;               errnadabscnhcf[i,j,1:ctgd,*] = nnaicferr / $
;                  (1-ionfrac)/10^(naabund - nadep)
               errnadabscnhcf[i,j,1:ctgd,*] = $
                  rebin(nadabscnhcf[i,j,1:ctgd],1,1,ctgd,2)*$
                  sqrt((nnaicferr/rebin(nnaicf,1,1,ctgd,2))^2d + $
                     (rebin(oneminusionfrac_relerr,1,1,ctgd,2))^2d)
               zdiff = reverse(tmpwaveabs[sortgd])/$
                       (nadlinelist['NaD1']*(1d + initdat.zsys_gas)) - 1d
               nadabsvel[i,j,1:ctgd] = c_kms * ((zdiff+1d)^2d - 1d) / $
                                               ((zdiff+1d)^2d + 1d)
               errnadabsvel[i,j,1:ctgd,*] = $
                  c_kms * (4d/(nadlinelist['NaD1']*(1d + initdat.zsys_gas))*$
                           (rebin(zdiff,1,1,ctgd,2)+1d)/$
                           ((rebin(zdiff,1,1,ctgd,2)+1d)^2d + 1d)^2d) * $
                           reverse(tmpwaveabserr[0,0,sortgd,*],3)
               nadabssig[i,j,1:ctgd] = reverse(tmpsigabs[sortgd])
               errnadabssig[i,j,1:ctgd,*] = reverse(tmpsigabserr[0,0,sortgd,*],3)
               nadabsv98[i,j,1:ctgd] = nadabsvel[i,j,1:ctgd]-$
                                       2d*nadabssig[i,j,1:ctgd]
;               errnadabsv98[i,j,1:ctgd] = $
;                  sqrt(errnadabsvel[i,j,1:ctgd]^2d + $
;                  4d*errnadabssig[i,j,1:ctgd]^2d)
               nadabscftau[i,j,1:ctgd] = reverse(tmpcftau[sortgd])
               errnadabscftau[i,j,1:ctgd,*] = reverse(tmpcftauerr[0,0,sortgd,*],3)
               nadabstau[i,j,1:ctgd] = reverse(tmptau[sortgd])
               errnadabstau[i,j,1:ctgd,0] = reverse(tmptauerrlo[sortgd])
               errnadabstau[i,j,1:ctgd,1] = reverse(tmptauerrhi[sortgd])
               nadabscf[i,j,1:ctgd] = reverse(tmpcf[sortgd])
            endif
;           Sort emission line wavelengths. In output velocity array,
;           first element of third dimension holds data for spaxels with only
;           1 component. Next elements hold velocities for spaxels with more than
;           1 comp, in order of increasing redshift.
            igd = where(nadfit.waveem[i,j,*] ne bad AND $
                        nadfit.waveem[i,j,*] ne 0,ctgd)
            if ctgd gt 0 then begin
               nademncomp[i,j] = ctgd
               tmpwaveem = nadfit.waveem[i,j,igd]
               tmpwaveemerr = mean(nadfit.waveemerr[i,j,*,*],dim=4)
               tmpwaveemerr = tmpwaveemerr[igd]
               tmpsigem = nadfit.sigmaem[i,j,igd]
               tmpsigemerr = mean(nadfit.sigmaemerr[i,j,*,*],dim=4)
               tmpsigemerr = tmpsigemerr[igd]
            endif
            if ctgd eq 1 then begin
               zdiff = tmpwaveem/$
                       (nadlinelist['NaD1']*(1d + initdat.zsys_gas)) - 1d
               nademvel[i,j,0] = c_kms * ((zdiff+1d)^2d - 1d) / $
                                 ((zdiff+1d)^2d + 1d)
               errnademvel[i,j,0] = $
                  c_kms * (4d/(nadlinelist['NaD1']*(1d + initdat.zsys_gas))*$
                           (zdiff+1d)/((zdiff+1d)^2d + 1d)^2d) * tmpwaveemerr
               nademsig[i,j,0] = tmpsigem
               errnadabssig[i,j,0] = tmpsigemerr
               nademv98[i,j,0] = nademvel[i,j,0]+2d*nademsig[i,j,0]
;               errnademv98[i,j,0] = $
;                  sqrt(errnademvel[i,j,0]^2d + 4d*errnademsig[i,j,0]^2d)
            endif else if ctgd gt 1 then begin
               sortgd = sort(tmpwaveem)
               zdiff = tmpwaveem[sortgd]/$
                       (nadlinelist['NaD1']*(1d + initdat.zsys_gas)) - 1d
               nademvel[i,j,1:ctgd] = c_kms * ((zdiff+1d)^2d - 1d) / $
                                      ((zdiff+1d)^2d + 1d)               
               errnademvel[i,j,1:ctgd] = $
                  c_kms * (4d/(nadlinelist['NaD1']*(1d + initdat.zsys_gas))*$
                           (zdiff+1d)/((zdiff+1d)^2d + 1d)^2d) * $
                  tmpwaveemerr[sortgd]
               nademsig[i,j,1:ctgd] = tmpsigem[sortgd]
               errnademsig[i,j,1:ctgd] = tmpsigemerr[sortgd]
               nademv98[i,j,1:ctgd] = nademvel[i,j,1:ctgd]+2d*nademsig[i,j,1:ctgd]
;               errnademv98[i,j,1:ctgd] = $
;                  sqrt(errnademvel[i,j,1:ctgd]^2d + $
;                  4d*errnademsig[i,j,1:ctgd]^2d)                                       
            endif
         endfor
      endfor
      
;     Parse actual numbers of NaD components
      ionecomp = where(nadabsncomp eq 1,ctonecomp)
      if ctonecomp gt 0 then donadabsonecomp=1b else donadabsonecomp = 0b      
      imulticomp = where(nadabsncomp gt 1 AND nadabsncomp ne bad,ctmulticomp)
      if ctmulticomp gt 0 then donadabsmulticomp=1b else donadabsmulticomp = 0b
      igd_tmp = where(nadabsncomp ne bad,ctgd_tmp)
      if ctgd_tmp gt 0 then maxnadabsncomp_act=max(nadabsncomp[igd_tmp]) $
      else maxnadabsncomp_act = 0

      ionecomp = where(nademncomp eq 1,ctonecomp)
      if ctonecomp gt 0 then donademonecomp=1b else donademonecomp = 0b      
      imulticomp = where(nademncomp gt 1 AND nademncomp ne bad,ctmulticomp)
      if ctmulticomp gt 0 then donademmulticomp=1b else donademmulticomp = 0b      
      igd_tmp = where(nademncomp ne bad,ctgd_tmp)
      if ctgd_tmp gt 0 then maxnademncomp_act=max(nademncomp[igd_tmp]) $
      else maxnademncomp_act = 0

;     Absorption lines: Cumulative velocity distribution functions
      nadabscvdf = $
         ifsf_cmpcvdf_abs(nadfit.waveabs,$
                          mean(nadfit.waveabserr,dim=4),$
                          nadfit.sigmaabs,$
                          mean(nadfit.sigmaabserr,dim=4),$
                          nadfit.tau,$
                          mean(nadfit.tauerr,dim=4),$
                          initnad.maxncomp,nadlinelist['NaD1'],$
                          initdat.zsys_gas)
      nadabscvdfvals = ifsf_cmpcvdfvals_abs(nadabscvdf)

;     Emission lines: Cumulative velocity distribution functions
      nademcvdf = $
         ifsf_cmpcvdf_abs(nadfit.waveem,$
                          mean(nadfit.waveemerr,dim=4),$
                          nadfit.sigmaem,$
                          mean(nadfit.sigmaemerr,dim=4),$
                          nadfit.flux,$
                          mean(nadfit.fluxerr,dim=4),$
                          initnad.maxncomp,nadlinelist['NaD1'],$
                          initdat.zsys_gas)
      nademcvdfvals = ifsf_cmpcvdfvals_abs(nademcvdf)


   endif else begin
      nadabsvel=0b
      nadabsv98=0b
      ibd_nadabs_fitweq=0b
      igd_nadabs_fitweq=0b
   endelse

;  Compass rose
   angarr_rad = initdat.positionangle*!DPi/180d
   sinangarr = sin(angarr_rad)
   cosangarr = cos(angarr_rad)
;  starting point and length of compass rose normalized to plot panel
   xarr0_norm = 0.95
   yarr0_norm = 0.05
   rarr_norm = 0.2
   rlaboff_norm = 0.05
   laboff_norm = 0.0
   carr = 'White'
;  Coordinates in kpc for arrow coordinates:
;  Element 1: starting point
;  2: end of N arrow
;  3: end of E arrow
;  4: N label
;  5: E label
   xarr_kpc = dblarr(5)
   yarr_kpc = dblarr(5)
;  average panel dimension
   pdim = ((xran_kpc[1]-xran_kpc[0]) + (yran_kpc[1]-yran_kpc[0]))/2d
   xarr_kpc[0] = xarr0_norm * (xran_kpc[1]-xran_kpc[0]) + xran_kpc[0]
   xarr_kpc[1] = xarr_kpc[0] + rarr_norm*pdim*sinangarr
   xarr_kpc[2] = xarr_kpc[0] - rarr_norm*pdim*cosangarr
   xarr_kpc[3] = xarr_kpc[0] + (rarr_norm+rlaboff_norm)*pdim*sinangarr
   xarr_kpc[4] = xarr_kpc[0] - (rarr_norm+rlaboff_norm)*pdim*cosangarr
   yarr_kpc[0] = yarr0_norm * (yran_kpc[1]-yran_kpc[0]) + yran_kpc[0]
   yarr_kpc[1] = yarr_kpc[0] + rarr_norm*pdim*cosangarr
   yarr_kpc[2] = yarr_kpc[0] + rarr_norm*pdim*sinangarr
   yarr_kpc[3] = yarr_kpc[0] + (rarr_norm+rlaboff_norm)*pdim*cosangarr
   yarr_kpc[4] = yarr_kpc[0] + (rarr_norm+rlaboff_norm)*Pdim*sinangarr

   minyarr_kpc = min(yarr_kpc)
   if minyarr_kpc lt yran_kpc[0] then yarr_kpc -= minyarr_kpc - yran_kpc[0]
   maxxarr_kpc = max(xarr_kpc)
   if maxxarr_kpc gt xran_kpc[1] then xarr_kpc -= maxxarr_kpc - xran_kpc[1]

;  Compute coordinates for disk axis lines
   diskaxes_endpoints = 0b
   if tag_exist(initmaps,'plotdiskaxes') then begin

      diskaxes_endpoints = dblarr(2,2,2)
      for i=0,1 do begin
         halflength=initmaps.plotdiskaxes.length[i]/2d
         sinangle_tmp = sin(initmaps.plotdiskaxes.angle[i]*!DPi/180d)
         cosangle_tmp = cos(initmaps.plotdiskaxes.angle[i]*!DPi/180d)
         xends = initmaps.plotdiskaxes.xcenter[i]+[halflength*sinangle_tmp,$
                 -halflength*sinangle_tmp]
         yends = initmaps.plotdiskaxes.ycenter[i]+[-halflength*cosangle_tmp,$
                 halflength*cosangle_tmp]
         diskaxes_endpoints[*,*,i] = [[xends],[yends]]
      endfor
   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Fit PSF to Emission Line Map
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;   if tag_exist(initmaps,'fit_empsf') then begin
;      linmap_tmp = $
;         emlflx[string(initmaps.fit_empsf.vel,format='(I0)'),$
;                initmaps.fit_empsf.line]
;      vel = initmaps.fit_empsf.vel
;      ivel = value_locate(linmap_tmp.vel,vel)
;      empsf_map = linmap_tmp.flux[*,*,ivel]
;      maxempsf_map = max(empsf_map,imax)
;      empsf_map /= maxempsf_map
;
;      ;     Use error in total flux for error in line
;      empsf_err = tlinmaps[initmaps.fit_empsf.line,*,*,0,1]
;      empsf_err /= empsf_err[imax]
;
;      ;     2D Moffat fit to continuum flux vs. radius
;      parinfo = REPLICATE({fixed:0b},8)
;      parinfo[0].fixed = 1b
;      est=[0d,1d,1d,1d,center_nuclei[0]-1d,center_nuclei[1]-1d,0d,2.5d]
;      empsf_fit = $
;         mpfit2dpeak(empsf_map,empsf_fitpar,/moffat,/circular,est=est,$
;         parinfo=parinfo,error=empsf_err)
;
;      map_rempsf = sqrt((map_x - empsf_fitpar[4]+1)^2d + $
;         (map_y - empsf_fitpar[5]+1)^2d)
;      map_rempsfkpc_ifs = map_rempsf * kpc_per_pix
;      empsf1d_x = dindgen(101)/100d*max(map_rempsfkpc_ifs)
;      empsf1d_y = alog10(moffat(empsf1d_x,[empsf_fitpar[1],0d,$
;         empsf_fitpar[2]*kpc_per_pix,$
;         empsf_fitpar[7]]))
;
;   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Continuum plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if tag_exist(initmaps,'ct') then begin
      if dohst then ctsumrange_tmp = initmaps.ct.sumrange_hstcomp $
      else ctsumrange_tmp = initmaps.ct.sumrange
      capifs = string(ctsumrange_tmp[0],'-',ctsumrange_tmp[1],$
                      format='(I0,A0,I0)')
      if tag_exist(initmaps.ct,'sumrange_lab') then begin
         if initmaps.ct.sumrange_lab eq 'microns' then $
            capifs = string(ctsumrange_tmp[0]/1d4,'-',ctsumrange_tmp[1]/1d4,$
                            format='(D0.2,A0,D0.2)')
      endif
      if tag_exist(initmaps.ct,'charscale') then charscale=initmaps.ct.charscale $
      else charscale = 1d
   endif
   if dohst then begin
      if tag_exist(initmaps.hst,'source') then capsource = initmaps.hst.source $
      else capsource = 'HST'
      if dohstrd AND dohstbl then $
         caphst = textoidl(initmaps.hstbl.label+'+'+initmaps.hstrd.label) $
      else if dohstrd then $
         caphst = textoidl(initmaps.hstrd.label) $
      else $
         caphst = textoidl(initmaps.hstbl.label)
   endif
   ;  arrays for positions for zoom box
   posbox1x = dblarr(2)
   posbox1y = dblarr(2)
   posbox2x = dblarr(2)
   posbox2y = dblarr(2)
   
;  Figure out correct image size in inches
   ysize_in = 2.2d
   aspectrat_fov=double(dx)/double(dy)
   npanels_ifsfov = 0
   if tag_exist(initmaps,'ct') then npanels_ifsfov = 1d
   if dohst then npanels_ifsfov += 1d
   if dohstsm then npanels_ifsfov += 1d
   if npanels_ifsfov eq 0 then begin
      print,'IFSF_MAKEMAPS: Error -- no continuum images to plot.'
      stop
   endif
   imgheight_in = 1.6d
   xmargin_in = 0.4d
   ymargin_in = (ysize_in - imgheight_in)/2d
   ifsimg_width_in = imgheight_in*aspectrat_fov*npanels_ifsfov
   ;  Sizes and positions of image windows in real and normalized coordinates
   if dohst then begin
      xsize_in = imgheight_in + xmargin_in + ifsimg_width_in
      xfrac_margin = xmargin_in / xsize_in
      xfrac_hstbig = imgheight_in / xsize_in
      xfrac_ifsfov_width = imgheight_in*aspectrat_fov / xsize_in
      yfracb = ymargin_in/ysize_in
      yfract = 1d - ymargin_in/ysize_in
      pos_hstbig = [0d,yfracb,xfrac_hstbig,yfract]
      pos_ifsfov = dblarr(4,fix(npanels_ifsfov))
      pos_ifsfov[*,0] = [xfrac_hstbig+xfrac_margin,yfracb,$
                         xfrac_hstbig+xfrac_margin+xfrac_ifsfov_width,yfract]
      for i=1,fix(npanels_ifsfov)-1 do $
         pos_ifsfov[*,i] = pos_ifsfov[*,i-1] + [xfrac_ifsfov_width,0d,$
                                                xfrac_ifsfov_width,0d]
;     Instrument labels
      lineoff = 0.1d*xfrac_hstbig
      xhstline = [pos_hstbig[0]+lineoff,$
                  pos_ifsfov[2,npanels_ifsfov-2]-lineoff]
      yhstline = [yfracb*0.75d,yfracb*0.75d]
      xhstline_tpos = (xhstline[1]+xhstline[0])/2d
      yhstline_tpos = yfracb*0.15d
      xifsline = [pos_ifsfov[0,npanels_ifsfov-1]+lineoff,$
                  pos_ifsfov[2,npanels_ifsfov-1]-lineoff]
      yifsline = [yfracb*0.75d,yfracb*0.75d]
      xifsline_tpos = (xifsline[1]+xifsline[0])/2d
      yifsline_tpos = yfracb*0.15d
   endif else begin
      ysize_in = 2.2d + ymargin_in
      xsize_in = xmargin_in + ifsimg_width_in
      yfracb = ymargin_in/ysize_in
      yfract = 1d - ymargin_in*2d/ysize_in
      xfrac_margin = xmargin_in/xsize_in
      yfrac_margin = ymargin_in/ysize_in
      pos_ifsfov = [xfrac_margin,yfracb,1d,yfract]
;     Instrument labels
      lineoff = 0.1d*ifsimg_width_in/xsize_in
      xifsline = [pos_ifsfov[0]+lineoff,$
                  pos_ifsfov[2]-lineoff]
      yifsline = [yfracb*0.75d,yfracb*0.75d]
      xifsline_tpos = (xifsline[1]+xifsline[0])/2d
      yifsline_tpos = yfracb*0.15d

   endelse

   cgps_open,initdat.mapdir+initdat.label+'cont.eps',$
             charsize=1d*charscale,$
             /encap,/inches,xs=xsize_in,ys=ysize_in,/qui,/nomatch

   if dohst then begin
      
      cgarrow,xhstline[0],yhstline[0],xhstline[1],yhstline[1],$
              /norm,thick=8,hsize=0,color='Red'
      cgtext,capsource+': '+caphst,xhstline_tpos,yhstline_tpos,chars=1d*charscale,$
             color='Red',/norm,align=0.5d
      
;     HST continuum, large scale
      if dohstbl then size_subim = size(bhst_big) $
      else size_subim = size(rhst_big)
      if dohstrd AND dohstbl then begin
          mapscl = bytarr(3,size_subim[1],size_subim[2])
          mapscl[0,*,*] = rhst_big
          mapscl[2,*,*] = bhst_big
          mapscl[1,*,*] = byte((double(rhst_big)+double(bhst_big))/2d)
          if size_subim[1] lt resampthresh OR size_subim[2] lt resampthresh then $
             mapscl = rebin(mapscl,3,size_subim[1]*samplefac,size_subim[2]*samplefac,/sample)
      endif else begin
         mapscl = bytarr(size_subim[1],size_subim[2])
         if dohstrd then mapscl = rhst_big
         if dohstbl then mapscl = bhst_big
         if size_subim[1] lt resampthresh OR size_subim[2] lt resampthresh then $
            mapscl = rebin(mapscl,size_subim[1]*samplefac,size_subim[2]*samplefac,/sample)
      endelse
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos_hstbig,opos=truepos,$
              /noerase,missing_value=bad,missing_index=255,$
              missing_color='white'
;     Assumes IFS FOV coordinates are 0-offset, with [0,0] at a pixel center
      cgplot,[0],xsty=5,ysty=5,xran=[-0.5,size_subim[1]-0.5],$
             yran=[-0.5,size_subim[2]-0.5],position=truepos,$
             /nodata,/noerase,title=initdat.name
      cgoplot,[hst_big_ifsfov[*,0],hst_big_ifsfov[0,0]],$
              [hst_big_ifsfov[*,1],hst_big_ifsfov[0,1]],color='Red'

      imsize = string(fix(initmaps.hst.subim_big*kpc_per_as),format='(I0)')
      cgtext,size_subim[1]*0.05,size_subim[2]*0.9,$
             textoidl(imsize+'\times'+imsize+' kpc'),$
             color='white'
;      cgtext,size_subim[1]*0.05,size_subim[2]*0.05,caphst,color='white'
      posbox1x[0] = truepos[0]+(truepos[2]-truepos[0])*$
                    hst_big_ifsfov[3,0]/size_subim[1]
      posbox1y[0] = truepos[1]+(truepos[3]-truepos[1])*$
                    hst_big_ifsfov[3,1]/size_subim[2]
      posbox2x[0] = truepos[0]+(truepos[2]-truepos[0])*$
                    hst_big_ifsfov[0,0]/size_subim[1]
      posbox2y[0] = truepos[1]+(truepos[3]-truepos[1])*$
                    hst_big_ifsfov[0,1]/size_subim[2]

;     HST continuum, IFS FOV
      if dohstbl then size_subim = size(bhst_fov) $
      else size_subim = size(rhst_fov_sc)
      if dohstbl AND dohstrd then begin
         mapscl = bytarr(3,size_subim[1],size_subim[2])
         mapscl[0,*,*] = rhst_fov_sc
         mapscl[2,*,*] = bhst_fov
         mapscl[1,*,*] = byte((double(rhst_fov_sc)+double(bhst_fov))/2d)
         ctmap = (rhst_fov_ns+bhst_fov_ns)/2d
         if size_subim[1] lt resampthresh OR size_subim[2] lt resampthresh then $
            mapscl = rebin(mapscl,3,size_subim[1]*samplefac,size_subim[2]*samplefac,/sample)
      endif else begin
         mapscl = bytarr(size_subim[1],size_subim[2])
         if dohstrd then begin
            mapscl = rhst_fov_sc
            ctmap = rhst_fov_ns
         endif else if dohstbl then begin
            mapscl = bhst_fov
            ctmap = bhst_fov_ns
         endif
         if size_subim[1] lt resampthresh OR size_subim[2] lt resampthresh then $
            mapscl = rebin(mapscl,size_subim[1]*samplefac,size_subim[2]*samplefac,/sample)
      endelse

      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos_ifsfov[*,0],opos=truepos,$
             /noerase,missing_value=bad,missing_index=255,$
             missing_color='white'
      if tag_exist(initmaps.hst,'fithstpeak') AND $
         tag_exist(initmaps.hst,'fithstpeakwin_kpc') then begin
         nucfit_dwin_kpc = initmaps.hst.fithstpeakwin_kpc
         nucfit_halfdwin_hstpix = round(nucfit_dwin_kpc/kpc_per_hstpix/2d)
;        subsets of images for peak fitting, centered around (first) nucleus
         xhst_sub = round(center_nuclei_hst[0,0]) + $
                    [-nucfit_halfdwin_hstpix,nucfit_halfdwin_hstpix]
         yhst_sub = round(center_nuclei_hst[1,0]) + $
                    [-nucfit_halfdwin_hstpix,nucfit_halfdwin_hstpix]
         ctmap_center = ctmap[xhst_sub[0]:xhst_sub[1],$
                              yhst_sub[0]:yhst_sub[1]]
;        Circular moffat fit
         yfit = mpfit2dpeak(ctmap_center,a,/moffat,/circular)
;        Fitted peak coordinate in HST pixels; single-offset coordinates,
;        [1,1] at a pixel center
         peakfit_hstpix = [a[4]+xhst_sub[0]+1,a[5]+yhst_sub[0]+1]
         peakfit_hst_distance_from_nucleus_hstpix = peakfit_hstpix - $
                                                    center_nuclei_hst[*,0]
         peakfit_hst_distance_from_nucleus_kpc = $
            peakfit_hst_distance_from_nucleus_hstpix * kpc_per_hstpix
         size_hstpix = size(ctmap)
         cgplot,[0],xsty=5,ysty=5,xran=[0.5,size_hstpix[1]+0.5],$
                yran=[0.5,size_hstpix[2]+0.5],position=truepos,$
                /nodata,/noerase
         cgoplot,peakfit_hstpix[0],peakfit_hstpix[1],psym=1,color='Red'
      endif else begin
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase
      endelse
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,$
                       center_nuclei_kpc_y,/toplab
      cgtext,xran_kpc[0]+(xran_kpc[1]-xran_kpc[0])*0.05,$
             yran_kpc[1]-(yran_kpc[1]-yran_kpc[0])*0.1,$
             'IFS FOV',/data,color='white'
      ifsf_plotcompass,xarr_kpc,yarr_kpc,carr=carr,/nolab,hsize=150d,hthick=2d
      posbox1x[1] = truepos[0]
      posbox1y[1] = truepos[3]
      posbox2x[1] = truepos[0]
      posbox2y[1] = truepos[1]

;     smoothed HST continuum, IFS FOV
      if dohstsm then begin
;;        3-color image
;         mapscl = bytarr(3,size_subim[1],size_subim[2])
;         if dohstrd then mapscl[0,*,*] = rhst_fov_sm
;         if dohstbl then mapscl[2,*,*] = bhst_fov_sm
;         if dohstrd AND dohstbl then $
;            mapscl[1,*,*] = byte((double(rhst_fov_sm)+double(bhst_fov_sm))/2d)
;        Flux image
         if dohstbl AND dohstrd then begin
            ctmap = (double(rhst_fov_sm_ns_rb)+double(bhst_fov_sm_ns_rb))/2d
            if tag_exist(initmaps.hstblsm,'beta') then $
               beta=initmaps.hstblsm.beta $
            else if tag_exist(initmaps.hstrdsm,'beta') then $
               beta=initmaps.hstrdsm.beta $
            else beta=1d
            if tag_exist(initmaps.hstblsm,'stretch') then $
               stretch=initmaps.hstblsm.stretch $
            else if tag_exist(initmaps.hstrdsm,'stretch') then $
               stretch=initmaps.hstrdsm.stretch $
            else stretch=1
            if tag_exist(initmaps.hstblsm,'scllim') then $
               scllim=initmaps.hstblsm.scllim $
            else if tag_exist(initmaps.hstrdsm,'scllim') then $
               scllim=initmaps.hstrdsm.scllim $
            else scllim=[min(ctmap),max(ctmap)]               
         endif else if dohstbl then begin
            ctmap = bhst_fov_sm_ns_rb
            if tag_exist(initmaps.hstblsm,'beta') then $
               beta=initmaps.hstblsm.beta $
            else beta=1d
            if tag_exist(initmaps.hstblsm,'stretch') then $
               stretch=initmaps.hstblsm.stretch $
            else stretch=1
            if tag_exist(initmaps.hstblsm,'scllim') then $
               scllim=initmaps.hstblsm.scllim $
            else scllim=[min(ctmap),max(ctmap)]
         endif else begin
            ctmap = rhst_fov_sm_ns_rb
            if tag_exist(initmaps.hstrdsm,'beta') then $
               beta=initmaps.hstrdsm.beta $
            else beta=1d
            if tag_exist(initmaps.hstrdsm,'stretch') then $
               stretch=initmaps.hstrdsm.stretch $
            else stretch=1
            if tag_exist(initmaps.hstrdsm,'scllim') then $
               scllim=initmaps.hstrdsm.scllim $
            else scllim=[min(ctmap),max(ctmap)]
         endelse
         ctmap /= max(ctmap)
         zran = scllim
         dzran = zran[1]-zran[0]
;         if tag_exist(initmaps.ct,'beta') then beta=initmaps.ct.beta else beta=1d
         mapscl = cgimgscl(ctmap,minval=zran[0],max=zran[1],$
                           stretch=stretch,beta=beta)
         if size_subim[1] lt resampthresh OR size_subim[2] lt resampthresh then $
            mapscl = rebin(mapscl,dx*samplefac,dy*samplefac,/sample)
         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos_ifsfov[*,1],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         if tag_exist(initmaps.ct,'fitifspeak') AND $
            tag_exist(initmaps.ct,'fitifspeakwin_kpc') then begin
            nucfit_dwin_kpc = initmaps.ct.fitifspeakwin_kpc
            nucfit_halfdwin_pix = round(nucfit_dwin_kpc/kpc_per_pix/2d)
;           subsets of images for peak fitting, centered around (first) nucleus
            x_sub = round(center_nuclei[0,0]) + $
                    [-nucfit_halfdwin_pix,nucfit_halfdwin_pix]
            y_sub = round(center_nuclei[1,0]) + $
                    [-nucfit_halfdwin_pix,nucfit_halfdwin_pix]
            ctmap_center = ctmap[x_sub[0]:x_sub[1],$
                                 y_sub[0]:y_sub[1]]
;           Circular moffat fit
            yfit = mpfit2dpeak(ctmap_center,a,/moffat,/circular)
;           Fitted peak coordinate in IFS pixels; single-offset coordinates,
;           [1,1] at a pixel center
            peakfit_pix = [a[4]+x_sub[0]+1,a[5]+y_sub[0]+1]
            peakfit_hstconv_distance_from_nucleus_pix = peakfit_pix - $
                                                        center_nuclei[*,0]
            peakfit_hstconv_distance_from_nucleus_kpc = $
               peakfit_hstconv_distance_from_nucleus_pix * kpc_per_pix
            cgplot,[0],xsty=5,ysty=5,xran=[0.5,dx+0.5],$
                   yran=[0.5,dy+0.5],position=truepos,$
                   /nodata,/noerase
            cgoplot,peakfit_pix[0],peakfit_pix[1],psym=1,color='Red'
         endif else begin
            cgplot,[0],xsty=5,ysty=5,position=truepos,$
                   /nodata,/noerase
         endelse
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,$
                          center_nuclei_kpc_x,center_nuclei_kpc_y,/nolab
         cgtext,xran_kpc[0]+(xran_kpc[1]-xran_kpc[0])*0.05,$
                yran_kpc[1]-(yran_kpc[1]-yran_kpc[0])*0.1,$
                'IFS FOV, conv.',/data,color='white'
      endif

      cgplot,posbox1x,posbox1y,color='Red',$
         xsty=5,ysty=5,/noerase,xran=[0,1],yran=[0,1],pos=[0,0,1,1]
      cgoplot,posbox2x,posbox2y,color='Red'

   endif

   if tag_exist(initmaps,'ct') then begin

      cgarrow,xifsline[0],yifsline[0],xifsline[1],yifsline[1],$
         /norm,thick=8,hsize=0,color='Blue'
      cgtext,'IFS',xifsline_tpos,yifsline_tpos,chars=1d*charscale,$
         color='Blue',/norm,align=0.5

      ictlo = value_locate(datacube.wave,ctsumrange_tmp[0])
      icthi = value_locate(datacube.wave,ctsumrange_tmp[1])
      zran = initmaps.ct.scllim
      dzran = zran[1]-zran[0]
      if tag_exist(initmaps.ct,'domedian') then $
         ctmap = median(datacube.dat[*,*,ictlo:icthi],dim=3,/double)*$
                 double(icthi-ictlo+1) $
      else $
         ctmap = total(datacube.dat[*,*,ictlo:icthi],3)
      ctmap /= max(ctmap)
      if tag_exist(initmaps.ct,'beta') then beta=initmaps.ct.beta else beta=1d
      mapscl = cgimgscl(rebin(ctmap,dx*samplefac,dy*samplefac,/sample),$
                        minval=zran[0],max=zran[1],$
                        stretch=initmaps.ct.stretch,beta=beta)                        
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos_ifsfov[*,fix(npanels_ifsfov) - 1],$
              opos=truepos,/noerase,missing_value=bad,missing_index=255,$
              missing_color='white'
      if tag_exist(initmaps.ct,'fitifspeak') AND $
         tag_exist(initmaps.ct,'fitifspeakwin_kpc') then begin
         nucfit_dwin_kpc = initmaps.ct.fitifspeakwin_kpc
         nucfit_halfdwin_pix = round(nucfit_dwin_kpc/kpc_per_pix/2d)
;           subsets of images for peak fitting, centered around (first) nucleus
         x_sub = round(center_nuclei[0,0]) + $
                 [-nucfit_halfdwin_pix,nucfit_halfdwin_pix]
         y_sub = round(center_nuclei[1,0]) + $
                 [-nucfit_halfdwin_pix,nucfit_halfdwin_pix]
         ctmap_center = ctmap[x_sub[0]:x_sub[1],$
                              y_sub[0]:y_sub[1]]
;           Circular moffat fit
         yfit = mpfit2dpeak(ctmap_center,a,/moffat,/circular)
;        Fitted peak coordinate in IFS pixels; single-offset coordinates,
;        [1,1] at a pixel center
         peakfit_pix = [a[4]+x_sub[0]+1,a[5]+y_sub[0]+1]
         peakfit_pix_ifs = peakfit_pix ; save for later
         peakfit_ifs_distance_from_nucleus_pix = peakfit_pix - $
                                                 center_nuclei[*,0]
         peakfit_ifs_distance_from_nucleus_kpc = $
            peakfit_ifs_distance_from_nucleus_pix * kpc_per_pix
         cgplot,[0],xsty=5,ysty=5,xran=[0.5,dx+0.5],$
            yran=[0.5,dy+0.5],position=truepos,$
            /nodata,/noerase
         cgoplot,peakfit_pix[0],peakfit_pix[1],psym=1,color='Red'
      endif else begin
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase,title=title_tmp
      endelse
      if npanels_ifsfov eq 1 then $
         cgtext,initdat.name,0.5,1d - yfrac_margin,$
                /norm,align=0.5,chars=1.25d*charscale
      if npanels_ifsfov gt 1 then begin
         nolab_tmp=1b
         toplab_tmp=0b
      endif else begin
         nolab_tmp=0b
         toplab_tmp=1b
      endelse
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,$
                       center_nuclei_kpc_x,center_nuclei_kpc_y,$
                       nolab=nolab_tmp,toplab=toplab_tmp
      cgtext,xran_kpc[0]+(xran_kpc[1]-xran_kpc[0])*0.05,$
             yran_kpc[1]-(yran_kpc[1]-yran_kpc[0])*0.1,$
             capifs,/data,color='white'
      if npanels_ifsfov eq 1 then $
         ifsf_plotcompass,xarr_kpc,yarr_kpc,carr=carr,/nolab,$
                          hsize=150d,hthick=2d

   endif

   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Continuum color plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if tag_exist(initmaps,'hstcol') OR tag_exist(initmaps,'hstcolsm') then begin

   if tag_exist(initmaps,'ct') then begin
      if dohst then ctsumrange_tmp = initmaps.ct.sumrange_hstcomp $
      else ctsumrange_tmp = initmaps.ct.sumrange
      capifs = string(ctsumrange_tmp[0],'-',ctsumrange_tmp[1],$
                      format='(I0,A0,I0)')
      if tag_exist(initmaps.ct,'sumrange_lab') then begin
         if initmaps.ct.sumrange_lab eq 'microns' then $
            capifs = string(ctsumrange_tmp[0]/1d4,'-',ctsumrange_tmp[1]/1d4,$
                            format='(D0.2,A0,D0.2)')
      endif
      if tag_exist(initmaps.ct,'charscale') then charscale=initmaps.ct.charscale $
      else charscale = 1d
   endif
   caphst = textoidl(initmaps.hstbl.label+'-'+initmaps.hstrd.label)
;  arrays for positions for zoom box
   posbox1x = dblarr(2)
   posbox1y = dblarr(2)
   posbox2x = dblarr(2)
   posbox2y = dblarr(2)
   
;  Figure out correct image size in inches
   ysize_in = 2.2d
   aspectrat_fov=double(dx)/double(dy)
   npanels_ifsfov = 0
   if tag_exist(initmaps,'ct') then npanels_ifsfov = 1d
   npanels_ifsfov += 1d
   if dohstcolsm then npanels_ifsfov += 1d
   if npanels_ifsfov eq 0 then begin
      print,'IFSF_MAKEMAPS: Error -- no continuum images to plot.'
      stop
   endif
   imgheight_in = 1.6d
   xmargin_in = 0.4d
   ymargin_in = (ysize_in - imgheight_in)/2d
   ifsimg_width_in = imgheight_in*aspectrat_fov*npanels_ifsfov
   ;  Sizes and positions of image windows in real and normalized coordinates
   if dohst then begin
      xsize_in = imgheight_in + xmargin_in + ifsimg_width_in
      xfrac_margin = xmargin_in / xsize_in
      xfrac_hstbig = imgheight_in / xsize_in
      xfrac_ifsfov_width = imgheight_in*aspectrat_fov / xsize_in
      yfracb = ymargin_in/ysize_in
      yfract = 1d - ymargin_in/ysize_in
      pos_hstbig = [0d,yfracb,xfrac_hstbig,yfract]
      pos_ifsfov = dblarr(4,fix(npanels_ifsfov))
      pos_ifsfov[*,0] = [xfrac_hstbig+xfrac_margin,yfracb,$
                         xfrac_hstbig+xfrac_margin+xfrac_ifsfov_width,yfract]
      for i=1,fix(npanels_ifsfov)-1 do $
         pos_ifsfov[*,i] = pos_ifsfov[*,i-1] + [xfrac_ifsfov_width,0d,$
                                                xfrac_ifsfov_width,0d]
;     Instrument labels
      lineoff = 0.1d*xfrac_hstbig
      xhstline = [pos_hstbig[0]+lineoff,$
                  pos_ifsfov[2,npanels_ifsfov-2]-lineoff]
      yhstline = [yfracb*0.75d,yfracb*0.75d]
      xhstline_tpos = (xhstline[1]+xhstline[0])/2d
      yhstline_tpos = yfracb*0.15d
      xifsline = [pos_ifsfov[0,npanels_ifsfov-1]+lineoff,$
                  pos_ifsfov[2,npanels_ifsfov-1]-lineoff]
      yifsline = [yfracb*0.75d,yfracb*0.75d]
      xifsline_tpos = (xifsline[1]+xifsline[0])/2d
      yifsline_tpos = yfracb*0.15d
   endif else begin
      ysize_in = 2.2d + ymargin_in
      xsize_in = xmargin_in + ifsimg_width_in
      yfracb = ymargin_in/ysize_in
      yfract = 1d - ymargin_in*2d/ysize_in
      xfrac_margin = xmargin_in/xsize_in
      yfrac_margin = ymargin_in/ysize_in
      pos_ifsfov = [xfrac_margin,yfracb,1d,yfract]
;     Instrument labels
      lineoff = 0.1d*ifsimg_width_in/xsize_in
      xifsline = [pos_ifsfov[0]+lineoff,$
                  pos_ifsfov[2]-lineoff]
      yifsline = [yfracb*0.75d,yfracb*0.75d]
      xifsline_tpos = (xifsline[1]+xifsline[0])/2d
      yifsline_tpos = yfracb*0.15d

   endelse

   cgps_open,initdat.mapdir+initdat.label+'color.eps',$
             charsize=1d*charscale,$
             /encap,/inches,xs=xsize_in,ys=ysize_in,/qui,/nomatch

   if dohst then begin
      
      cgarrow,xhstline[0],yhstline[0],xhstline[1],yhstline[1],$
              /norm,thick=8,hsize=0,color='Red'
      cgtext,capsource+': '+caphst,xhstline_tpos,yhstline_tpos,chars=1d*charscale,$
             color='Red',/norm,align=0.5d
      
;     HST continuum, large scale
      size_subim = size(chst_big)
      mapscl = chst_big
      if size_subim[1] lt resampthresh OR size_subim[2] lt resampthresh then $
         mapscl = rebin(mapscl,size_subim[1]*samplefac,size_subim[2]*samplefac,/sample)
      cgloadct,65
      cgimage,mapscl,/keep,pos=pos_hstbig,opos=truepos,$
              /noerase,missing_value=bad,missing_index=255,$
              missing_color='white'
;     Assumes IFS FOV coordinates are 0-offset, with [0,0] at a pixel center
      cgplot,[0],xsty=5,ysty=5,xran=[-0.5,size_subim[1]-0.5],$
             yran=[-0.5,size_subim[2]-0.5],position=truepos,$
             /nodata,/noerase,title=initdat.name
      cgoplot,[hst_big_ifsfov[*,0],hst_big_ifsfov[0,0]],$
              [hst_big_ifsfov[*,1],hst_big_ifsfov[0,1]],color='Red'

      imsize = string(fix(initmaps.hst.subim_big*kpc_per_as),format='(I0)')
      cgtext,size_subim[1]*0.05,size_subim[2]*0.9,$
             textoidl(imsize+'\times'+imsize+' kpc'),$
             color='white'
      posbox1x[0] = truepos[0]+(truepos[2]-truepos[0])*$
                    hst_big_ifsfov[3,0]/size_subim[1]
      posbox1y[0] = truepos[1]+(truepos[3]-truepos[1])*$
                    hst_big_ifsfov[3,1]/size_subim[2]
      posbox2x[0] = truepos[0]+(truepos[2]-truepos[0])*$
                    hst_big_ifsfov[0,0]/size_subim[1]
      posbox2y[0] = truepos[1]+(truepos[3]-truepos[1])*$
                    hst_big_ifsfov[0,1]/size_subim[2]

;     HST continuum, IFS FOV
      size_subim = size(chst_fov)
      mapscl = chst_fov
      if size_subim[1] lt resampthresh OR size_subim[2] lt resampthresh then $
         mapscl = rebin(mapscl,size_subim[1]*samplefac,size_subim[2]*samplefac,/sample)
      cgloadct,65
      cgimage,mapscl,/keep,pos=pos_ifsfov[*,0],opos=truepos,$
             /noerase,missing_value=bad,missing_index=255,$
             missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
             /nodata,/noerase
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,$
                       center_nuclei_kpc_y,/toplab
      cgtext,xran_kpc[0]+(xran_kpc[1]-xran_kpc[0])*0.05,$
             yran_kpc[1]-(yran_kpc[1]-yran_kpc[0])*0.1,$
             'IFS FOV',/data,color='white'
      ifsf_plotcompass,xarr_kpc,yarr_kpc,carr=carr,/nolab,hsize=150d,hthick=2d
      posbox1x[1] = truepos[0]
      posbox1y[1] = truepos[3]
      posbox2x[1] = truepos[0]
      posbox2y[1] = truepos[1]

;     smoothed HST continuum, IFS FOV
      if dohstcolsm then begin
         map = cshst_fov_rb
         size_subim = size(map)
         zran = initmaps.hstcolsm.scllim
         dzran = zran[1]-zran[0]
         mapscl = cgimgscl(map,minval=zran[0],max=zran[1],$
                           stretch=initmaps.ct.stretch)
         if size_subim[1] lt resampthresh OR size_subim[2] lt resampthresh then $
            mapscl = rebin(mapscl,size_subim[1]*samplefac,size_subim[2]*samplefac,/sample)
         cgloadct,65
         cgimage,mapscl,/keep,pos=pos_ifsfov[*,1],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,$
                          center_nuclei_kpc_x,center_nuclei_kpc_y,/nolab
         cgtext,xran_kpc[0]+(xran_kpc[1]-xran_kpc[0])*0.05,$
                yran_kpc[1]-(yran_kpc[1]-yran_kpc[0])*0.1,$
                'IFS FOV, conv.',/data,color='white'
      endif

      cgplot,posbox1x,posbox1y,color='Red',$
         xsty=5,ysty=5,/noerase,xran=[0,1],yran=[0,1],pos=[0,0,1,1]
      cgoplot,posbox2x,posbox2y,color='Red'

   endif

   if tag_exist(initmaps,'ct') then begin

      cgarrow,xifsline[0],yifsline[0],xifsline[1],yifsline[1],$
         /norm,thick=8,hsize=0,color='Blue'
      cgtext,'IFS',xifsline_tpos,yifsline_tpos,chars=1d*charscale,$
         color='Blue',/norm,align=0.5

      ictlo = value_locate(datacube.wave,ctsumrange_tmp[0])
      icthi = value_locate(datacube.wave,ctsumrange_tmp[1])
      zran = initmaps.ct.scllim
      dzran = zran[1]-zran[0]
      if tag_exist(initmaps.ct,'domedian') then $
         ctmap = median(datacube.dat[*,*,ictlo:icthi],dim=3,/double)*$
                 double(icthi-ictlo+1) $
      else $
         ctmap = total(datacube.dat[*,*,ictlo:icthi],3)
      ctmap /= max(ctmap)
      if tag_exist(initmaps.ct,'beta') then beta=initmaps.ct.beta else beta=1d
      mapscl = cgimgscl(rebin(ctmap,dx*samplefac,dy*samplefac,/sample),$
                        minval=zran[0],max=zran[1],$
                        stretch=initmaps.ct.stretch,beta=beta)                        
      cgloadct,65
      cgimage,mapscl,/keep,pos=pos_ifsfov[*,fix(npanels_ifsfov) - 1],$
              opos=truepos,/noerase,missing_value=bad,missing_index=255,$
              missing_color='white'
      if tag_exist(initmaps.ct,'fitifspeak') AND $
         tag_exist(initmaps.ct,'fitifspeakwin_kpc') then begin
         nucfit_dwin_kpc = initmaps.ct.fitifspeakwin_kpc
         nucfit_halfdwin_pix = round(nucfit_dwin_kpc/kpc_per_pix/2d)
;           subsets of images for peak fitting, centered around (first) nucleus
         x_sub = round(center_nuclei[0,0]) + $
                 [-nucfit_halfdwin_pix,nucfit_halfdwin_pix]
         y_sub = round(center_nuclei[1,0]) + $
                 [-nucfit_halfdwin_pix,nucfit_halfdwin_pix]
         ctmap_center = ctmap[x_sub[0]:x_sub[1],$
                              y_sub[0]:y_sub[1]]
;           Circular moffat fit
         yfit = mpfit2dpeak(ctmap_center,a,/moffat,/circular)
;        Fitted peak coordinate in IFS pixels; single-offset coordinates,
;        [1,1] at a pixel center
         peakfit_pix = [a[4]+x_sub[0]+1,a[5]+y_sub[0]+1]
         peakfit_pix_ifs = peakfit_pix ; save for later
         peakfit_ifs_distance_from_nucleus_pix = peakfit_pix - $
                                                 center_nuclei[*,0]
         peakfit_ifs_distance_from_nucleus_kpc = $
            peakfit_ifs_distance_from_nucleus_pix * kpc_per_pix
         cgplot,[0],xsty=5,ysty=5,xran=[0.5,dx+0.5],$
            yran=[0.5,dy+0.5],position=truepos,$
            /nodata,/noerase
         cgoplot,peakfit_pix[0],peakfit_pix[1],psym=1,color='Red'
      endif else begin
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase,title=title_tmp
      endelse
      if npanels_ifsfov eq 1 then $
         cgtext,initdat.name,0.5,1d - yfrac_margin,$
                /norm,align=0.5,chars=1.25d*charscale
      if npanels_ifsfov gt 1 then begin
         nolab_tmp=1b
         toplab_tmp=0b
      endif else begin
         nolab_tmp=0b
         toplab_tmp=1b
      endelse
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,$
                       center_nuclei_kpc_x,center_nuclei_kpc_y,$
                       nolab=nolab_tmp,toplab=toplab_tmp
      cgtext,xran_kpc[0]+(xran_kpc[1]-xran_kpc[0])*0.05,$
             yran_kpc[1]-(yran_kpc[1]-yran_kpc[0])*0.1,$
             capifs,/data,color='white'
      if npanels_ifsfov eq 1 then $
         ifsf_plotcompass,xarr_kpc,yarr_kpc,carr=carr,/nolab,$
                          hsize=150d,hthick=2d

   endif

   cgps_close
   
   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Continuum radial profiles
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if tag_exist(initmaps,'ct') then begin
      
      npy = 2
      npx = 1
      if tag_exist(initdat,'decompose_qso_fit') then npx = 3
      if tag_exist(initdat,'remove_scattered') then npx = 4

;  Figure out correct image size in inches
      panel_in = 2d
      margin_in = 0.5d
      halfmargin_in = margin_in/2d
      xsize_in = panel_in * npx + margin_in
      aspectrat_fov=double(dx)/double(dy)
      ysize_in = margin_in*2d + panel_in * (1d + 1d/aspectrat_fov)
;  Sizes and positions of image windows in real and normalized coordinates
      pan_xfrac = panel_in/xsize_in
      pan_yfrac = panel_in/ysize_in
      ifs_yfrac = panel_in/aspectrat_fov/ysize_in
      mar_xfrac = margin_in/xsize_in
      mar_yfrac = margin_in/ysize_in
      hmar_xfrac = halfmargin_in/xsize_in
      hmar_yfrac = halfmargin_in/ysize_in
      pos_top = dblarr(4,npx)
      pos_bot = dblarr(4,npx)
      for i=0,npx-1 do begin
         pos_top[*,i] = [mar_xfrac+double(i)*pan_xfrac,$
                         1d - (hmar_yfrac+pan_yfrac),$
                         mar_xfrac+double(i+1)*pan_xfrac,$
                         1d - hmar_yfrac]
         pos_bot[*,i] = [mar_xfrac+double(i)*pan_xfrac,$
                         hmar_yfrac,$
                         mar_xfrac+double(i+1)*pan_xfrac,$
                         hmar_yfrac+ifs_yfrac]
      endfor

      zran = initmaps.ct.scllim_rad
;      dzran = zran[1]-zran[0]

      cgps_open,initdat.mapdir+initdat.label+'cont_rad.eps',charsize=1,/encap,$
                /inches,xs=xsize_in,ys=ysize_in,/qui,/nomatch
      
;     Total (continuum-only) model flux. Plot fits if decompose tags set, otherwise
;     plot total cube flux within specified range.
      if tag_exist(initdat,'decompose_qso_fit') then begin
;        Divide model flux by total # pixels for cases where total number of 
;        pixels varies by spaxel (due to contracted spectra at edges, e.g.)
         ctmap = total(contcube.qso_mod+contcube.host_mod,3) / contcube.npts
         ctsumrange_tmp = initdat.fitran
      endif else if tag_exist(initdat,'decompose_ppxf_fit') then begin
         tmpstel = contcube.stel_mod_tot
         ibdtmp = where(tmpstel eq bad,ctbdtmp)
         if ctbdtmp gt 0 then tmpstel[ibdtmp] = 0d
         tmppoly = contcube.poly_mod_tot
         ibdtmp = where(tmppoly eq bad,ctbdtmp)
         if ctbdtmp gt 0 then tmppoly[ibdtmp] = 0d
         ctmap = tmpstel+tmppoly
         ctsumrange_tmp = initdat.fitran
      endif else begin
         ictlo = value_locate(datacube.wave,initmaps.ct.sumrange[0])
         icthi = value_locate(datacube.wave,initmaps.ct.sumrange[1])
         ctmap = total(datacube.dat[*,*,ictlo:icthi],3)
         ctsumrange_tmp = initmaps.ct.sumrange
      endelse
      capran = string(ctsumrange_tmp[0],'-',ctsumrange_tmp[1],$
                      format='(I0,A0,I0)')
      if tag_exist(initmaps.ct,'sumrange_lab') then begin
         if initmaps.ct.sumrange_lab eq 'microns' then $
            capran = string(ctsumrange_tmp[0]/1d4,'-',ctsumrange_tmp[1]/1d4,$
                            format='(D0.2,A0,D0.2)')
      endif
      maxctmap = max(ctmap)
      ctmap /= maxctmap
      cgplot,map_rkpc_ifs,alog10(ctmap),yran=[-4,0],$
             xran=[0,max(map_rkpc_ifs)],/xsty,/ysty,psym=16,symsize=0.5d,$
             pos=pos_top[*,0],aspect=1d,title=initdat.name,$
             xtit = 'Radius (kpc)',ytit = 'log I/I$\downmax$'
      if tag_exist(initdat,'decompose_qso_fit') then begin
         cgoplot,psf1d_x,psf1d_y,color='Red'
      endif else if tag_exist(initmaps,'fit_empsf') then begin
         cgoplot,empsf1d_x,empsf1d_y,color='Red'
      endif else if tag_exist(initmaps,'ctradprof_psffwhm') then begin
         x = dindgen(101)/100d*max(map_rkpc_ifs)
         fwhm=initmaps.ctradprof_psffwhm * kpc_per_as
;        Gaussian
         y = alog10(gaussian(x,[1d,0d,fwhm/2.35]))
         cgoplot,x,y,color='Black'
;        Moffat, index = 1.5
         y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/1.5d)-1),1.5d]))
         cgoplot,x,y,color='Red',/linesty
;        Moffat, index = 2.5
         y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/2.5d)-1),2.5d]))
         cgoplot,x,y,color='Red'
;        Moffat, index = 5
         y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/5d)-1),5d]))
         cgoplot,x,y,color='Blue'
      endif     

;      mapscl = cgimgscl(rebin(ctmap,dx*samplefac,dy*samplefac,/sample),$
;                        minval=zran[0],max=zran[1],stretch=initmaps.ct.stretch)
      mapscl = cgimgscl(rebin(alog10(ctmap),dx*samplefac,dy*samplefac,/sample),$
                        minval=zran[0],max=zran[1],stretch=initmaps.ct.stretch)
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos_bot[*,0],opos=truepos,$
              /noerase,missing_value=bad,missing_index=255,$
              missing_color='white'
      cgplot,[0],xsty=5,ysty=5,xran=[0.5,dx+0.5],$
              yran=[0.5,dy+0.5],position=truepos,$
              /nodata,/noerase
      if tag_exist(initmaps.ct,'fitifspeak') AND $
         tag_exist(initmaps.ct,'fitifspeakwin_kpc') then $
          cgoplot,peakfit_pix_ifs[0],peakfit_pix_ifs[1],psym=1,color='Red'         
      cgtext,dx*0.05,dy*0.95,'Host Cont.+Quasar',color='white'
      cgtext,dx*0.05,dy*0.05,capran,/data,color='white'
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y

      if tag_exist(initdat,'decompose_qso_fit') then begin
         qso_map = total(contcube.qso_mod,3) / contcube.npts
;         qso_map /= max(qso_map)
         qso_map /= maxctmap
         cgplot,map_rkpc_ifs,alog10(qso_map),yran=[-4,0],$
                xran=[0,max(map_rkpc_ifs)],/xsty,/ysty,psym=16,symsize=0.5d,$
                pos=pos_top[*,1],/noerase,aspect=1d,$
                xtickformat='(A1)',ytickformat='(A1)'
;         if tag_exist(initdat,'decompose_qso_fit') then begin
         cgoplot,psf1d_x,psf1d_y,color='Red'
         cgtext,max(map_rkpc_ifs)*0.9d,-4d*0.1d,'FWHM='+$
                string(qso_fitpar[2]*initdat.platescale,format='(D0.2)')+$
                ' asec',align=1d
         cgtext,max(map_rkpc_ifs)*0.9d,-4d*0.2d,'FWHM='+$
                string(qso_fitpar[2]*initdat.platescale*kpc_per_as,format='(D0.2)')+$
                ' kpc',align=1d
         cgtext,max(map_rkpc_ifs)*0.9d,-4d*0.3d,'$\gamma$='+$
                string(qso_fitpar[7],format='(D0.1)'),align=1d
;         endif else if tag_exist(initmaps,'fit_empsf') then begin
;            cgoplot,empsf1d_x,empsf1d_y,color='Red'
;         endif else if tag_exist(initmaps,'ctradprof_psffwhm') then begin
;            x = dindgen(101)/100d*max(map_rkpc_ifs)
;            fwhm=initmaps.ctradprof_psffwhm * kpc_per_as
;;           Gaussian
;            y = alog10(gaussian(x,[1d,0d,fwhm/2.35]))
;            cgoplot,x,y,color='Black'
;;           Moffat, index = 1.5
;            y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/1.5d)-1),1.5d]))
;            cgoplot,x,y,color='Red',/linesty
;;           Moffat, index = 2.5
;            y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/2.5d)-1),2.5d]))
;            cgoplot,x,y,color='Red'
;;           Moffat, index = 5
;            y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/5d)-1),5d]))
;            cgoplot,x,y,color='Blue'
;         endif

;         mapscl = cgimgscl(rebin(qso_map,dx*samplefac,dy*samplefac,/sample),$
;                           minval=zran[0],max=zran[1],stretch=initmaps.ct.stretch)
         mapscl = cgimgscl(rebin(alog10(qso_map),dx*samplefac,dy*samplefac,/sample),$
                        minval=zran[0],max=zran[1],stretch=initmaps.ct.stretch)
         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos_bot[*,1],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,xran=[0.5,dx+0.5],$
                yran=[0.5,dy+0.5],position=truepos,$
                /nodata,/noerase
         if tag_exist(initmaps.ct,'fitifspeak') AND $
            tag_exist(initmaps.ct,'fitifspeakwin_kpc') then begin
            peakfit_pix = [qso_fitpar[4]+1,qso_fitpar[5]+1]
            peakfit_ifs_qso_distance_from_nucleus_pix = peakfit_pix - $
                                                        center_nuclei[*,0]
            peakfit_ifs_qso_distance_from_nucleus_kpc = $
                peakfit_ifs_qso_distance_from_nucleus_pix * kpc_per_pix
            cgoplot,peakfit_pix[0],peakfit_pix[1],psym=1,color='Red'         
         endif
         cgtext,dx*0.05,dy*0.95,'Quasar PSF',color='white'
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,$
                          center_nuclei_kpc_x,center_nuclei_kpc_y,/nolab

         host_map = total(contcube.host_mod,3) / contcube.npts
;         host_map /= max(host_map)
         host_map /= maxctmap
         cgplot,map_rkpc_ifs,alog10(host_map),yran=[-4,0],$
                xran=[0,max(map_rkpc_ifs)],/xsty,/ysty,psym=16,symsize=0.5d,$
                pos=pos_top[*,2],/noerase,aspect=1d,$
                xtickformat='(A1)',ytickformat='(A1)'
;         if tag_exist(initdat,'decompose_qso_fit') then begin
         cgoplot,psf1d_x,psf1d_y,color='Red'
;         endif else if tag_exist(initmaps,'fit_empsf') then begin
;            cgoplot,empsf1d_x,empsf1d_y,color='Red'
;         endif else if tag_exist(initmaps,'ctradprof_psffwhm') then begin
;            x = dindgen(101)/100d*max(map_rkpc_ifs)
;            fwhm=initmaps.ctradprof_psffwhm * kpc_per_as
;;           Gaussian
;            y = alog10(gaussian(x,[1d,0d,fwhm/2.35]))
;            cgoplot,x,y,color='Black'
;;           Moffat, index = 1.5
;            y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/1.5d)-1),1.5d]))
;            cgoplot,x,y,color='Red',/linesty
;;           Moffat, index = 2.5
;            y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/2.5d)-1),2.5d]))
;            cgoplot,x,y,color='Red'
;;           Moffat, index = 5
;            y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/5d)-1),5d]))
;            cgoplot,x,y,color='Blue'
;         endif
         
;         mapscl = cgimgscl(rebin(host_map,dx*samplefac,dy*samplefac,/sample),$
;                           minval=zran[0]),max=zran[1],stretch=initmaps.ct.stretch)
         mapscl = cgimgscl(rebin(alog10(host_map),dx*samplefac,dy*samplefac,/sample),$
                        minval=zran[0],max=zran[1],stretch=initmaps.ct.stretch)
         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos_bot[*,2],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,xran=[0.5,dx+0.5],$
                yran=[0.5,dy+0.5],position=truepos,$
                /nodata,/noerase
         cgtext,dx*0.05,dy*0.95,'Host Cont.',color='white'
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,$
                          center_nuclei_kpc_x,center_nuclei_kpc_y,/nolab

         if tag_exist(initdat,'remove_scattered') then begin
            
            scatt_map = total(contcube.poly_mod,3) / contcube.npts
;;           Use maximum flux for normalization unless it's much higher than 
;;           second highest flux
;            ifsort = reverse(sort(scatt_map))
;            if scatt_map[ifsort[0]]/scatt_map[ifsort[1]] gt 2 then $
;               scatt_map /= scatt_map[ifsort[1]] $
;            else scatt_map /= scatt_map[ifsort[0]]
            scatt_map /= maxctmap
            cgplot,map_rkpc_ifs,alog10(scatt_map),yran=[-4,0],$
                   xran=[0,max(map_rkpc_ifs)],/xsty,/ysty,psym=16,symsize=0.5d,$
                   pos=pos_top[*,3],/noerase,aspect=1d,$
                   xtickformat='(A1)',ytickformat='(A1)'
            cgoplot,psf1d_x,psf1d_y,color='Red'

;            mapscl = cgimgscl(rebin(scatt_map,dx*samplefac,dy*samplefac,/sample),$
;                              minval=zran[0],max=zran[1],$
;                              stretch=initmaps.ct.stretch)
            mapscl = cgimgscl(rebin(alog10(scatt_map),dx*samplefac,dy*samplefac,/sample),$
                        minval=zran[0],max=zran[1],stretch=initmaps.ct.stretch)
            cgloadct,65,/reverse
            cgimage,mapscl,/keep,pos=pos_bot[*,3],opos=truepos,$
                    /noerase,missing_value=bad,missing_index=255,$
                    missing_color='white'
            cgplot,[0],xsty=5,ysty=5,xran=[0.5,dx+0.5],$
                   yran=[0.5,dy+0.5],position=truepos,$
                   /nodata,/noerase
            cgtext,dx*0.05,dy*0.95,'Scattered Light',color='white'
            ifsf_plotaxesnuc,xran_kpc,yran_kpc,$
                             center_nuclei_kpc_x,center_nuclei_kpc_y,/nolab

         endif

      endif

      cgps_close

   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Continuum plots: Stellar velocity / sigma + E(B-V)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  Size of plot grid
   if tag_exist(initdat,'decompose_qso_fit') OR $
      tag_exist(initdat,'decompose_ppxf_fit') then npx = 2 else npx = 1
   if tag_exist(initdat,'ebv_star') then npx++
   if tag_exist(initdat,'mcerrors') then npy = 2 else npy = 1

;  Figure out correct image size in inches
   xpanel_in = 1.75d
   margin_in = 0.5d
   topmargin_in = 0.5d
   halfmargin_in = margin_in/2d
   xsize_in = xpanel_in*double(npx) + margin_in
   aspectrat_fov=double(dxwin)/double(dywin)
   ysize_in = xpanel_in/aspectrat_fov*double(npy) + $
              margin_in*(1.5d + (double(npy)-1)) + topmargin_in
;  Sizes and positions of image windows in real and normalized coordinates
   pan_xfrac = xpanel_in/xsize_in
   pan_yfrac = xpanel_in/aspectrat_fov/ysize_in
   sqpan_yfrac = xpanel_in/ysize_in
   mar_xfrac = margin_in/xsize_in
   mar_yfrac = margin_in/ysize_in
   topmar_yfrac = topmargin_in/ysize_in
   hmar_xfrac = halfmargin_in/xsize_in
   hmar_yfrac = halfmargin_in/ysize_in
   pos_top = dblarr(4,npx)
   pos_bot = dblarr(4,npx)
   pos_toptit = dblarr(2,npx)
   for i=0,npx-1 do begin
      if npy eq 1 then $
         pos_top[*,i] = [mar_xfrac+double(i)*pan_xfrac,$
                         mar_yfrac,$
                         mar_xfrac+double(i+1)*pan_xfrac,$
                         hmar_yfrac + pan_yfrac] $
      else begin
         pos_top[*,i] = [mar_xfrac+double(i)*pan_xfrac,$
                         1d - (mar_yfrac*0.5d + pan_yfrac + topmar_yfrac),$
                         mar_xfrac+double(i+1)*pan_xfrac,$
                         1d - (mar_yfrac*0.5d + topmar_yfrac)]
         pos_bot[*,i] = [mar_xfrac+double(i)*pan_xfrac,$
                         mar_yfrac,$
                         mar_xfrac+double(i+1)*pan_xfrac,$
                         mar_yfrac + pan_yfrac]
      endelse
   endfor
   pos_title = [mar_xfrac+double(npx)*pan_xfrac/2d,$
                1d - topmar_yfrac*0.75d]

   cgps_open,initdat.mapdir+initdat.label+'stel.eps',$
             charsize=1d*charscale,$
             /encap,/inches,xs=xsize_in,ys=ysize_in,/qui,/nomatch

   cbform = '(I0)' ; colorbar syntax
   stel_z = contcube.stel_z
;  For backwards compatibility ...
   if tag_exist(contcube,'stel_z_err') then begin
      stel_z_errlo = contcube.stel_z_err[*,*,0]
      stel_z_errhi = contcube.stel_z_err[*,*,1]
   endif else begin
      stel_z_errlo = stel_z * 0d
      stel_z_errhi = stel_z * 0d
   endelse
   igd = where(stel_z ne bad,ctgd)

   stel_vel = 0b
   if ctgd gt 0 then begin

;     Stellar velocity
      zdiff = dblarr(dx,dy)+bad
      stel_vel = dblarr(dx,dy)+bad
      stel_errvello = dblarr(dx,dy)+bad
      stel_errvelhi = dblarr(dx,dy)+bad
      zdiff[igd] = contcube.stel_z[igd] - initdat.zsys_gas
      stel_vel[igd] = c_kms * ((zdiff[igd]+1d)^2d - 1d) / $
                                 ((zdiff[igd]+1d)^2d + 1d)
      stel_errvello[igd] = stel_z_errlo[igd]*c_kms
      stel_errvelhi[igd] = stel_z_errhi[igd]*c_kms
      stel_errvel = [[[stel_errvello]],[[stel_errvelhi]]]
      
      map = stel_vel
      maperr = mean(stel_errvel,dim=3,/doub)

;     Set up range
      if hasrangefile then begin
         ithisline = where(rangeline eq 'stel' AND $
                           rangequant eq 'vel',ctthisline)
         if ctthisline eq 1 then auto=0b
      endif else auto=1b
      plotdat = ifsf_plotrange(auto=auto,$
                               mapgd=map[igd],divinit=100d,$
                               ncbdivmax=ncbdivmax,$
                               rline=rangeline,matline='stel',$
                               rcomp=rangecomp,$
                               rquant=rangequant,matquant='vel',$
                               rncbdiv=rangencbdiv,$
                               rlo=rangelo,rhi=rangehi)

      mapscl = bytscl(rebin(map[plotwin[0]-1:plotwin[2]-1,$
                                plotwin[1]-1:plotwin[3]-1],$
                            dxwin*samplefac,dywin*samplefac,/sample),$
                      min=plotdat[0],max=plotdat[1])
      cgloadct,74,/reverse
      cgimage,mapscl,/keep,pos=pos_top[*,0],opos=truepos,$
              missing_value=bad,missing_index=255,$
              missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
             /nodata,/noerase,xran=[0,dxwin],yran=[0,dywin]
;     Velocity contours
      if tag_exist(initmaps,'contourlevels') then begin
         key = 'stel_vel'
         if initmaps.contourlevels->haskey(key) then begin
            nlevels = n_elements(initmaps.contourlevels[key])
            cgcontour,map,dindgen(dxwin)+0.5,dindgen(dywin)+0.5,$
                      /overplot,color=0,c_linesty=2,c_thick=4,$
                      levels=initmaps.contourlevels[key],$
                      max=1000d
         endif
      endif
      ifsf_plotaxesnuc,xwinran_kpc,ywinran_kpc,center_nuclei_kpc_xwin,$
                       center_nuclei_kpc_ywin,/toplab
                       
;     Colorbar
      xoffset = pan_xfrac*0.05
      yoffset = mar_yfrac*0.2
      cbpos=[truepos[0]+xoffset,truepos[1]-yoffset,$
             truepos[2]-xoffset,truepos[1]]
      ticknames = string(dindgen(plotdat[3]+1)*$
                            plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
      ticknames[0]=' '
      ticknames[plotdat[3]]=' '
      cgcolorbar,position=cbpos,divisions=plotdat[3],$
                 ticknames=ticknames,charsize=0.6
;     Title
      title='v$\down50$'
      yoffset = mar_yfrac*0.85
      cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
             truepos[1]-yoffset,title,charsize=1.25,align=0.5,/norm

      if npy eq 2 then begin
         if hasrangefile then begin
            ithisline = where(rangeline eq 'stel' AND $
                              rangequant eq 'vel_err',ctthisline)
            if ctthisline eq 1 then auto=0b
         endif else auto=1b
         plotdat = ifsf_plotrange(auto=auto,$
                                  mapgd=maperr[igd],divinit=10d,$
                                  ncbdivmax=ncbdivmax,$
                                  rline=rangeline,matline='stel',$
                                  rcomp=rangecomp,$
                                  rquant=rangequant,matquant='vel_err',$
                                  rncbdiv=rangencbdiv,$
                                  rlo=rangelo,rhi=rangehi)
         mapscl = bytscl(rebin(maperr[plotwin[0]-1:plotwin[2]-1,$
                                      plotwin[1]-1:plotwin[3]-1],$
                               dxwin*samplefac,dywin*samplefac,/sample),$
                         min=plotdat[0],max=plotdat[1])
         cgloadct,74,/reverse
         cgimage,mapscl,/keep,pos=pos_bot[*,0],opos=truepos,$
                 missing_value=bad,missing_index=255,$
                 missing_color='white',/noerase
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase,xran=[0,dxwin],yran=[0,dywin]
         ifsf_plotaxesnuc,xwinran_kpc,ywinran_kpc,center_nuclei_kpc_xwin,$
                          center_nuclei_kpc_ywin,/noxlab
;        Colorbar
         xoffset = pan_xfrac*0.05
         yoffset = mar_yfrac*0.2
         cbpos=[truepos[0]+xoffset,truepos[1]-yoffset,$
                truepos[2]-xoffset,truepos[1]]
         ticknames = string(dindgen(plotdat[3]+1)*$
                            plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         ticknames[0]=' '
         ticknames[plotdat[3]]=' '
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,charsize=0.6
;        Title
         title='v$\down50$ err'
         yoffset = mar_yfrac*0.85
         cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
                truepos[1]-yoffset,title,charsize=1.25,align=0.5,/norm
      endif


   endif

   if tag_exist(initdat,'decompose_ppxf_fit') then begin

      map = contcube.stel_sigma
      maperr = contcube.stel_sigma_err
      igd = where(map ne bad,ctgd)
      maperr = mean(maperr,dim=3,/doub)

      if ctgd gt 0 then begin

;        Set up range
         if hasrangefile then begin
            ithisline = where(rangeline eq 'stel' AND $
                              rangequant eq 'vsig',ctthisline)
            if ctthisline eq 1 then auto=0b
         endif else auto=1b
         plotdat = ifsf_plotrange(auto=auto,$
                                  mapgd=map[igd],divinit=100d,$
                                  ncbdivmax=ncbdivmax,$
                                  rline=rangeline,matline='stel',$
                                  rcomp=rangecomp,$
                                  rquant=rangequant,matquant='vsig',$
                                  rncbdiv=rangencbdiv,$
                                  rlo=rangelo,rhi=rangehi)

         mapscl = bytscl(rebin(map[plotwin[0]-1:plotwin[2]-1,$
                                   plotwin[1]-1:plotwin[3]-1],$
                               dxwin*samplefac,dywin*samplefac,/sample),$
                         min=plotdat[0],max=plotdat[1])
         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos_top[*,1],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase,xran=[0,dxwin],yran=[0,dywin]
         ifsf_plotaxesnuc,xwinran_kpc,ywinran_kpc,center_nuclei_kpc_xwin,$
                          center_nuclei_kpc_ywin,/nolab
;        Colorbar
         xoffset = pan_xfrac*0.05
         yoffset = mar_yfrac*0.2
         cbpos=[truepos[0]+xoffset,truepos[1]-yoffset,$
                truepos[2]-xoffset,truepos[1]]
         ticknames = string(dindgen(plotdat[3]+1)*$
                            plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         ticknames[0]=' '
         ticknames[plotdat[3]]=' '
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,charsize=0.6
;        Title
         title='$\sigma$'
         yoffset = mar_yfrac*0.85
         cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
                truepos[1]-yoffset,title,charsize=1.25,align=0.5,/norm

         if npy eq 2 then begin
            if hasrangefile then begin
               ithisline = where(rangeline eq 'stel' AND $
                                 rangequant eq 'vsig_err',ctthisline)
            if ctthisline eq 1 then auto=0b
            endif else auto=1b
            plotdat = ifsf_plotrange(auto=auto,$
                                     mapgd=maperr[igd],divinit=10d,$
                                     ncbdivmax=ncbdivmax,$
                                     rline=rangeline,matline='stel',$
                                     rcomp=rangecomp,$
                                     rquant=rangequant,matquant='vsig_err',$
                                     rncbdiv=rangencbdiv,$
                                     rlo=rangelo,rhi=rangehi)
            mapscl = bytscl(rebin(maperr[plotwin[0]-1:plotwin[2]-1,$
                                         plotwin[1]-1:plotwin[3]-1],$
                                  dxwin*samplefac,dywin*samplefac,/sample),$
                            min=plotdat[0],max=plotdat[1])
            cgloadct,74,/reverse
            cgimage,mapscl,/keep,pos=pos_bot[*,1],opos=truepos,$
                    missing_value=bad,missing_index=255,$
                    missing_color='white',/noerase
            cgplot,[0],xsty=5,ysty=5,position=truepos,$
                   /nodata,/noerase,xran=[0,dxwin],yran=[0,dywin]
            ifsf_plotaxesnuc,xwinran_kpc,ywinran_kpc,center_nuclei_kpc_xwin,$
                             center_nuclei_kpc_ywin,/nolab
;           Colorbar
            xoffset = pan_xfrac*0.05
            yoffset = mar_yfrac*0.2
            cbpos=[truepos[0]+xoffset,truepos[1]-yoffset,$
                   truepos[2]-xoffset,truepos[1]]
            ticknames = string(dindgen(plotdat[3]+1)*$
                               plotdat[2]/double(plotdat[3]) - $
                               (plotdat[2] - plotdat[1]),format=cbform)
            ticknames[0]=' '
            ticknames[plotdat[3]]=' '
            cgcolorbar,position=cbpos,divisions=plotdat[3],$
                       ticknames=ticknames,charsize=0.6
;           Title
            title='$\sigma$ err'
            yoffset = mar_yfrac*0.85
            cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
                   truepos[1]-yoffset,title,charsize=1.25,align=0.5,/norm
         endif
      endif
   endif

   if tag_exist(initdat,'ebv_star') then begin

      map = contcube.stel_ebv
      igd = where(map ne bad,ctgd)
      stel_ebv = map
      stel_errebv = contcube.stel_ebv_err
      maperr = mean(stel_errebv,dim=3,/doub)
      cbform = '(D0.1)' ; colorbar syntax

      if ctgd gt 0 then begin

;        Set up range
         if hasrangefile then begin
            ithisline = where(rangeline eq 'stel' AND $
                              rangequant eq 'ebv',ctthisline)
            if ctthisline eq 1 then auto=0b
         endif else auto=1b
         plotdat = ifsf_plotrange(auto=auto,$
                                  mapgd=map[igd],divinit=0.2d,$
                                  ncbdivmax=ncbdivmax,$
                                  rline=rangeline,matline='stel',$
                                  rcomp=rangecomp,$
                                  rquant=rangequant,matquant='ebv',$
                                  rncbdiv=rangencbdiv,$
                                  rlo=rangelo,rhi=rangehi)

         plotdat_stel_ebv = plotdat

         mapscl = bytscl(rebin(map[plotwin[0]-1:plotwin[2]-1,$
                                   plotwin[1]-1:plotwin[3]-1],$
                               dxwin*samplefac,dywin*samplefac,/sample),$
                         min=plotdat[0],max=plotdat[1])
         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos_top[*,npx-1],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase,xran=[0,dxwin],yran=[0,dywin]
         ifsf_plotaxesnuc,xwinran_kpc,ywinran_kpc,center_nuclei_kpc_xwin,$
                          center_nuclei_kpc_ywin,/nolab
;        Colorbar
         xoffset = pan_xfrac*0.05
         yoffset = mar_yfrac*0.2
         cbpos=[truepos[0]+xoffset,truepos[1]-yoffset,$
                truepos[2]-xoffset,truepos[1]]
         ticknames = string(dindgen(plotdat[3]+1)*$
                            plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         ticknames[0]=' '
         ticknames[plotdat[3]]=' '
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,charsize=0.6
;        Title
         title='E(B-V)'
         yoffset = mar_yfrac*0.85
         cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
                truepos[1]-yoffset,title,charsize=1.25,align=0.5,/norm

         if npy eq 2 then begin
            if hasrangefile then begin
               ithisline = where(rangeline eq 'stel' AND $
                                 rangequant eq 'ebv_err',ctthisline)
            if ctthisline eq 1 then auto=0b
            endif else auto=1b
            plotdat = ifsf_plotrange(auto=auto,$
                                     mapgd=maperr[igd],divinit=0.05d,$
                                     ncbdivmax=ncbdivmax,$
                                     rline=rangeline,matline='stel',$
                                     rcomp=rangecomp,$
                                     rquant=rangequant,matquant='ebv_err',$
                                     rncbdiv=rangencbdiv,$
                                     rlo=rangelo,rhi=rangehi)
            mapscl = bytscl(rebin(maperr[plotwin[0]-1:plotwin[2]-1,$
                                         plotwin[1]-1:plotwin[3]-1],$
                                  dxwin*samplefac,dywin*samplefac,/sample),$
                            min=plotdat[0],max=plotdat[1])
            cgloadct,74,/reverse
            cgimage,mapscl,/keep,pos=pos_bot[*,npx-1],opos=truepos,$
                    missing_value=bad,missing_index=255,$
                    missing_color='white',/noerase
            cgplot,[0],xsty=5,ysty=5,position=truepos,$
                   /nodata,/noerase,xran=[0,dxwin],yran=[0,dywin]
            ifsf_plotaxesnuc,xwinran_kpc,ywinran_kpc,center_nuclei_kpc_xwin,$
                             center_nuclei_kpc_ywin,/nolab
;           Colorbar
            xoffset = pan_xfrac*0.05
            yoffset = mar_yfrac*0.2
            cbpos=[truepos[0]+xoffset,truepos[1]-yoffset,$
                   truepos[2]-xoffset,truepos[1]]
            ticknames = string(dindgen(plotdat[3]+1)*$
                               plotdat[2]/double(plotdat[3]) - $
                               (plotdat[2] - plotdat[1]),format='(D0.2)')
            ticknames[0]=' '
            ticknames[plotdat[3]]=' '
            cgcolorbar,position=cbpos,divisions=plotdat[3],$
                       ticknames=ticknames,charsize=0.6
;           Title
            title='E(B-V) err'
            yoffset = mar_yfrac*0.85
            cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
                   truepos[1]-yoffset,title,charsize=1.25,align=0.5,/norm
         endif
 
      endif
   endif else begin
      stel_ebv = 0d
      stel_errebv = 0d
   endelse

;  Title
   cgplot,[0],xsty=5,ysty=5,position=[0,0,1,1],/nodata,/noerase
   cgtext,pos_title[0],pos_title[1],initdat.name+': Stellar',$
          charsize=1.25d,align=0.5

   cgps_close


   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; EMISSION LINE PLOTS
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if ~ tag_exist(initdat,'noemlinfit') then begin
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plots of individual emission lines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if tag_exist(initmaps,'emlplot') then begin

;     Quantities to plot
      if tag_exist(initmaps.emlplot,'ftags') then ftags = initmaps.emlplot.ftags $
      else ftags = ['ftot','fpk']
      if tag_exist(initmaps.emlplot,'vtags') then vtags = initmaps.emlplot.vtags $
      else vtags = ['vsig','vpk']
      if tag_exist(initmaps.emlplot,'ftitles') then ftitles = initmaps.emlplot.ftitles $
      else ftitles = ['F$\downtot$','F$\downpk$']
      if tag_exist(initmaps.emlplot,'vtitles') then vtitles = initmaps.emlplot.vtitles $
      else vtitles = ['$\sigma$','v$\downpeak$']

;     Size of plot grid
      npx = max([n_elements(vtags),n_elements(ftags)])
      npy = 2

;     Figure out correct image size in inches
;      xpanel_in = 1.75d
      xpanel_in = 1.4d
      margin_in = 0.5d
      topmargin_in = 0.5d
      halfmargin_in = margin_in/1d
      xsize_in = xpanel_in*double(npx) + margin_in
      aspectrat_fov=double(dxwin)/double(dywin)
      ysize_in = (xpanel_in * 1d/aspectrat_fov + margin_in)*2d + topmargin_in
;     Sizes and positions of image windows in real and normalized coordinates
      pan_xfrac = xpanel_in/xsize_in
      pan_yfrac = xpanel_in/aspectrat_fov/ysize_in
      sqpan_yfrac = xpanel_in/ysize_in
      mar_xfrac = margin_in/xsize_in
      mar_yfrac = margin_in/ysize_in
      topmar_yfrac = topmargin_in/ysize_in
      hmar_xfrac = halfmargin_in/xsize_in
      hmar_yfrac = halfmargin_in/ysize_in
      pos_top = dblarr(4,npx)
      pos_toptit = dblarr(2,npx)
      pos_mid = dblarr(4,npx)
      pos_midtit = dblarr(2,npx)
      pos_bot = dblarr(4,npx)
      for i=0,npx-1 do begin
         pos_top[*,i] = [mar_xfrac+double(i)*pan_xfrac,$
                         1d - (mar_yfrac+pan_yfrac) - topmar_yfrac,$
                         mar_xfrac+double(i+1)*pan_xfrac,$
                         1d - mar_yfrac - topmar_yfrac]
         pos_bot[*,i] = [mar_xfrac+double(i)*pan_xfrac,$
                         mar_yfrac,$
                         mar_xfrac+double(i+1)*pan_xfrac,$
                         mar_yfrac+pan_yfrac]
      endfor
      pos_title = [mar_xfrac+double(npx)*pan_xfrac/2d,$
                   1d - topmar_yfrac*0.75d]
        
;     Loop through emission lines
      foreach line,lines_with_doublets do begin

         linelab = ifsf_linesyntax(line)
         cgps_open,initdat.mapdir+initdat.label+linelab+ '.eps',$
                   charsize=1,/encap,/inches,xs=xsize_in,ys=ysize_in,/qui,/nomatch

;        loop through plot types
         for j=0,npx-1 do begin

;           FLUXES

            cbform = '(D0.1)' ; colorbar syntax
            iplot = j ; plot index
            if j lt n_elements(ftags) then begin

;              Get map and scale
               map = emlflx[ftags[j],line]
               ibd = where(map eq bad OR map eq 0 OR ~ finite(map),ctbd)
               igd = where(map ne bad AND map ne 0 AND finite(map),ctgd)
               ctgdlin = ctgd ; saving for later ...
            
               if ctgd gt 0 then begin
               
                  if tag_exist(initmaps,'fluxfactor') then $
                     map[igd] *= initmaps.fluxfactor
                  if tag_exist(initmaps,'fluxunits') then $
                     map[igd] *= initmaps.fluxunits
                  if tag_exist(initmaps,'vornorm') then $
                     map[igd] /= initmaps.vornorm[igd]
                                 
                  zran=[0,1]
                  dzran = 1
                  ncbdiv = 5
                  zmax_flux = max(map[igd])

                  if hasrangefile then begin
                     ithisline = where(line eq rangeline AND $
                                       ftags[j] eq rangequant,ctthisline)
                     if ctthisline eq 1 then zmax_flux = rangehi[ithisline]
                  endif
               
                  map[igd] = map[igd]/zmax_flux[0]
                  if ctbd gt 0 then map[ibd] = bad
                  mapscl = bytscl(rebin(map[plotwin[0]-1:plotwin[2]-1,$
                                            plotwin[1]-1:plotwin[3]-1],$
                                  dxwin*samplefac,dywin*samplefac,/sample),$
                                  min=zran[0],max=zran[1])
   
;                 Plot image
                  cgloadct,65,/reverse
                  cgimage,mapscl,/keep,pos=pos_top[*,j],opos=truepos,$
                          noerase=iplot ne 0,missing_value=bad,missing_index=255,$
                          missing_color='white'
                  cgplot,[0],xsty=5,ysty=5,position=truepos,$
                         /nodata,/noerase,xran=[0,dxwin],yran=[0,dywin]
;                 Disk axes
                  linsty_da=[2,1]
                  if diskaxes_endpoints[0] ne 0b then $
                     for k=0,1 do cgoplot,diskaxes_endpoints[*,0,k]+0.5d,$
                                          diskaxes_endpoints[*,1,k]+0.5d,$
                                          thick=4,linesty=linsty_da[k]
                  if j eq 0 then $
                     ifsf_plotaxesnuc,xwinran_kpc,ywinran_kpc,center_nuclei_kpc_xwin,$
                                      center_nuclei_kpc_ywin,/toplab $
                  else $
                     ifsf_plotaxesnuc,xwinran_kpc,ywinran_kpc,center_nuclei_kpc_xwin,$
                                      center_nuclei_kpc_ywin,/nolab
;                 Colorbar
                  if j eq 1 then begin
                     xoffset = pan_xfrac*0.05
                     yoffset = mar_yfrac*0.2
                     cbpos=[truepos[0]+xoffset,truepos[3],truepos[2]-xoffset,$
                            truepos[3]+yoffset]
                     ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
                                        (dzran - zran[1]),format=cbform)
                     cgcolorbar,position=cbpos,divisions=ncbdiv,$
                                ticknames=ticknames,/top,charsize=0.6
                  endif
 
;                 Title
                  title=ftitles[j]
                  title += ' ('+string(zmax_flux,format='(E0.1)')+')'
                  yoffset = mar_yfrac*0.65
                  cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
                         truepos[3]+yoffset,title,charsize=0.85,align=0.5,/norm
 
               endif
            endif else ctgdlin=0

;           VELOCITIES
                          
            cbform = '(I0)' ; colorbar syntax

            if j lt n_elements(vtags) then begin

               iplot = npx+j ; plot index

;              Get map and scale
               if emlvel[vtags[j]].haskey(line) then begin
                  map = emlvel[vtags[j],line]
                  ibd = where(map eq bad OR ~ finite(map),ctbd)
                  igd = where(map ne bad AND map ne 0 AND finite(map),ctgd)
               endif else ctgd = 0

               if ctgd gt 0 then begin
 
                  hasrange = 0
                  if hasrangefile then begin
                     ithisline = where(line eq rangeline AND $
                                       vtags[j] eq rangequant,ctthisline)
                     if ctthisline eq 1 then begin
                        zran = [rangelo[ithisline],rangehi[ithisline]]
                        dzran = zran[1]-zran[0]
                        ncbdiv = rangencbdiv[ithisline]                   
                        ncbdiv = ncbdiv[0]
                        hasrange = 1
                     endif
                  endif
                  if ~hasrange then begin
                     zran = [min(map[igd]),max(map[igd])]
                     divarr = ifsf_cbdiv(zran,100d,ncbdivmax)
                     ncbdiv = divarr[0]
                     dzran = zran[1]-zran[0]
                  endif
 
                  if ctbd gt 0 then map[ibd] = bad
                  mapscl = bytscl(rebin(map[plotwin[0]-1:plotwin[2]-1,$
                                            plotwin[1]-1:plotwin[3]-1],$
                                  dxwin*samplefac,dywin*samplefac,/sample),$
                                  min=zran[0],max=zran[1])
   
;                 Plot image
                  if stregex(vtags[j],'sig',/bool) then cgloadct,65,/reverse $
                  else cgloadct,74,/reverse
                  cgimage,mapscl,/keep,pos=pos_bot[*,j],opos=truepos,$
                          noerase=iplot ne 0,missing_value=bad,missing_index=255,$
                          missing_color='white'
                  cgplot,[0],xsty=5,ysty=5,position=truepos,$
                         /nodata,/noerase,xran=[0,dxwin],yran=[0,dywin]
;                 Velocity contours
                  if tag_exist(initmaps,'contourlevels') then begin
                     key = line+'_'+vtags[j]
;                    Not sure why levels aren't being labeled
                     if initmaps.contourlevels->haskey(key) then begin
                        nlevels = n_elements(initmaps.contourlevels[key])
                        if tag_exist(initmaps,'argscontour') then $
                           cgcontour,map[plotwin[0]-1:plotwin[2]-1,$
                                         plotwin[1]-1:plotwin[3]-1],$
                                     dindgen(dxwin)+0.5,dindgen(dywin)+0.5,$
                                     /overplot,$
                                     levels=initmaps.contourlevels[key],$
                                     _extra = initmaps.argscontour $
                        else $
                           cgcontour,map[plotwin[0]-1:plotwin[2]-1,$
                                         plotwin[1]-1:plotwin[3]-1],$
                                     dindgen(dxwin)+0.5,dindgen(dywin)+0.5,$
                                     /overplot,color=0,c_linesty=2,c_thick=4,$
                                     levels=initmaps.contourlevels[key],$
                                     max=1000d
                     endif
                  endif
;;                 Cross section
;                  if xsec_endpoints[0] ne 0b then $
;                     for k=0,n_elements(initmaps.xsec.line)-1 do $
;                        if initmaps.xsec.line[k] eq line AND $
;                           initmaps.xsec.tag[k] eq vtags[j] then $
;                           cgoplot,xsec_endpoints[*,0,k]+0.5d,$
;                                   xsec_endpoints[*,1,k]+0.5d,$
;                                   thick=4,linesty=0
;;                 Disk axes
;                  linsty_da=[2,1]
;                  if diskaxes_endpoints[0] ne 0b then $
;                     for k=0,1 do cgoplot,diskaxes_endpoints[*,0,k]+0.5d,$
;                                          diskaxes_endpoints[*,1,k]+0.5d,$
;                                          thick=4,linesty=linsty_da[k]
                  ifsf_plotaxesnuc,xwinran_kpc,ywinran_kpc,center_nuclei_kpc_xwin,$
                                   center_nuclei_kpc_ywin,/nolab
;                 Colorbar
                  xoffset = pan_xfrac*0.05
                  yoffset = mar_yfrac*0.2
                  cbpos=[truepos[0]+xoffset,truepos[1]-yoffset,$
                         truepos[2]-xoffset,truepos[1]]
                  ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
                                     (dzran - zran[1]),format=cbform)
                  ticknames[0]=' '
                  ticknames[ncbdiv]=' '
                  cgcolorbar,position=cbpos,divisions=ncbdiv,$
                             ticknames=ticknames,charsize=0.6
;                 Title
                  title=vtitles[j]
                  yoffset = mar_yfrac*0.85
                  cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
                         truepos[1]-yoffset,title,charsize=1.25,align=0.5,/norm

               endif               
            endif
 
         endfor

;        Title
         cgplot,[0],xsty=5,ysty=5,position=[0,0,1,1],/nodata,/noerase
         cgtext,pos_title[0],pos_title[1],initdat.name+': '+linelabels[line],$
                charsize=1.25d,align=0.5,/normal ; 2020may06, DSNR, added /normal
 
         cgps_close

      endforeach

   endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Extinction plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if tag_exist(initmaps,'ebv') then begin
      if tag_exist(initmaps.ebv,'calc') then begin

         plotdat_ebv = hash()

;        Loop through types
         for i=0,n_elements(initmaps.ebv.calc)-1 do begin
            
            fluxtype = initmaps.ebv.calc[i]
            
            map = ebv[fluxtype]
            maperr = errebv[fluxtype]
            ibd = where(map eq bad OR ~ finite(map),ctbd)
            igd = where(map ne bad AND finite(map),ctgd)

            if ctgd gt 0 then begin

               npx = 1
;              Figure out correct image size in inches
               xpanel_in = 1.75d
               margin_in = 0.5d
               topmargin_in = 0.5d
               halfmargin_in = margin_in/1d
               xsize_in = xpanel_in*double(npx) + margin_in
               aspectrat_fov=double(dx)/double(dy)
               ysize_in = xpanel_in * 1d/aspectrat_fov + margin_in + topmargin_in
;              Sizes and positions of image windows in real and normalized coordinates
               pan_xfrac = xpanel_in/xsize_in
               pan_yfrac = xpanel_in/aspectrat_fov/ysize_in
               sqpan_yfrac = xpanel_in/ysize_in
               mar_xfrac = margin_in/xsize_in
               mar_yfrac = margin_in/ysize_in
               topmar_yfrac = topmargin_in/ysize_in
               hmar_xfrac = halfmargin_in/xsize_in
               hmar_yfrac = halfmargin_in/ysize_in
               pos = dblarr(4,npx)
               pos_tit = dblarr(2,npx)
               for j=0,npx-1 do $
                  pos[*,j] = [mar_xfrac+double(j)*pan_xfrac,$
                              1d - (mar_yfrac+pan_yfrac) - topmar_yfrac,$
                              mar_xfrac+double(j+1)*pan_xfrac,$
                              1d - mar_yfrac - topmar_yfrac]
               pos_title = [mar_xfrac+double(npx)*pan_xfrac/2d,$
                            1d - topmar_yfrac*0.75d]

;              Set up range
               if hasrangefile then begin
                  ithisline = where(rangeline eq 'ebv' AND $
                                    rangequant eq fluxtype,ctthisline)
                  if ctthisline eq 1 then auto=0b
               endif else auto=1b
               plotdat = $
                  ifsf_plotrange(auto=auto,$
                                 mapgd=map[igd],divinit=0.5d,$
                                 ncbdivmax=ncbdivmax,$
                                 rline=rangeline,matline='ebv',$
                                 rcomp=rangecomp,$
                                 rquant=rangequant,matquant=fluxtype,$
                                 rncbdiv=rangencbdiv,$
                                 rlo=rangelo,rhi=rangehi)

               plotdat_ebv[fluxtype] = plotdat

               cgps_open,initdat.mapdir+initdat.label+'ebv_'+fluxtype+'.eps',$
                         charsize=1,/encap,/inches,xs=xsize_in,ys=ysize_in,/qui,/nomatch
               cbform = '(D0.2)'

               mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
                               min=plotdat[0],max=plotdat[1])
               cgloadct,74,/reverse
               cgimage,mapscl,/keep,pos=pos,opos=truepos,$
                       /noerase,missing_value=bad,missing_index=255,$
                       missing_color='white'
               cgplot,[0],xsty=5,ysty=5,position=truepos,$
                      /nodata,/noerase
               ifsf_plotaxesnuc,xran_kpc,yran_kpc,$
                                center_nuclei_kpc_x,center_nuclei_kpc_y
               xoffset = pan_xfrac*0.05
               yoffset = mar_yfrac*0.2
               cbpos=[truepos[0]+xoffset,truepos[3],truepos[2]-xoffset,$
                      truepos[3]+yoffset]
               ticknames = string(dindgen(plotdat[3]+1)*$
                                  plotdat[2]/double(plotdat[3]) - $
                                  (plotdat[2] - plotdat[1]),format=cbform)
               cgcolorbar,position=cbpos,divisions=plotdat[3],$
                          ticknames=ticknames,/top,charsize=0.6
               yoffset = mar_yfrac*0.65
               cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
                      truepos[3]+yoffset,initmaps.ebv.titles[i],$
                      charsize=1.05,align=0.5,/norm


;              Title
               cgplot,[0],xsty=5,ysty=5,position=[0,0,1,1],/nodata,/noerase
               cgtext,pos_title[0],pos_title[1],$
                      initdat.name+': E(B-V)',$
                      charsize=1d,align=0.5

               cgps_close

            endif else begin
               plotdat_ebv[fluxtype] = 0b
            endelse
         endfor
      endif
   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Emission line vs. stellar extinction
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


   if tag_exist(initmaps,'ebv') AND $
      tag_exist(initdat,'ebv_star') then begin
      if tag_exist(initmaps.ebv,'calc') then begin
         for i=0,n_elements(initmaps.ebv.calc)-1 do begin
            
            fluxtype = initmaps.ebv.calc[i]
            map = ebv[fluxtype]
            maperr = errebv[fluxtype]
            ibd = where(map eq bad OR ~ finite(map),ctbd)
            igd = where(map ne bad AND finite(map),ctgd)

            if ctgd gt 0 then begin
               cgps_open,initdat.mapdir+initdat.label+'ebv_'+fluxtype+'_v_ebv_stel.eps',$
                         charsize=1,/encap,/inches,xs=7.5d,ys=7.5d,/qui,/nomatch
               xran = [min([plotdat_ebv[fluxtype,0],plotdat_stel_ebv[0]]),$
                       max([plotdat_ebv[fluxtype,1],plotdat_stel_ebv[1]])]
               yran=xran
               cgplot,[0],/xsty,/ysty,/nodata,xran=xran,yran=yran,$
                      xtit='gas E(B-V)',ytit='stellar E(B-V)'
               cgoplot,map,stel_ebv,psym=16,symsize=0.75d,$
                       err_xlow=maperr,err_xhi=maperr,$
                       err_ylow=stel_errebv[*,*,0],$
                       err_yhi=stel_errebv[*,*,1]
               cgoplot,xran,yran
; relationship between gas and stellar E(B-V) from Calzetti, 1997AIPC..408..403C
               cgoplot,xran,yran*0.44d,/linesty
               cgps_close
            endif

         endfor
      endif
   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Stellar extinction vs. HST color
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


   if dohstcolsm AND tag_exist(initdat,'ebv_star') then begin

      map = cshst_fov_rb
      igd = where(cshst_fov_rb ne bad AND stel_ebv ne bad)
      xran = [min(map[igd])*0.95d,max(map[igd])*1.05d]
      yran = plotdat_stel_ebv[0:1]
      if xran[0] lt 0 then xran[0]=0d

      cgps_open,initdat.mapdir+initdat.label+'color_v_ebv_stel.eps',$
                charsize=1,/encap,/inches,xs=7.5d,ys=7.5d,/qui,/nomatch
      cgplot,[0],/xsty,/ysty,/nodata,xran=xran,yran=yran,$
             xtit=initmaps.hstbl.label+'-'+initmaps.hstrd.label,$
             ytit='stellar E(B-V)',title=initdat.name
      cgoplot,map,stel_ebv,psym=16,symsize=0.75d,$
              err_ylow=stel_errebv[*,*,0],$
              err_yhi=stel_errebv[*,*,1],/err_clip
;     Calzetti magnitude difference between red and blue continuum filters,
;     assuming R_V = 4.05 (though other R_V give same answer) and E(B-V)=1
      calz_unred,[pivotbl,pivotrd],[1d,1d],1d,funred
;     y = mx+b, where y = E(B-V) and x = blue-red color
;     Slope is 1/magnitude difference for E(B-V)=1
      m = 1d/(2.5d*(alog10(funred[0])-alog10(funred[1])))
;     Solve for b using point y=0 for xintrinsic.
      if tag_exist(initdat,'intrinsic_color') then $
         xint = initdat.intrinsic_color $
      else xint = min(map[igd])
      b = -m*xint
      xmod = [min(map[igd]),max(map[igd])]
      ymod = m*xmod+b
      cgoplot,xmod,ymod
      cgps_close

   endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Gas extinction vs. HST color
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


   if dohstcolsm AND tag_exist(initmaps,'ebv') then begin
      if tag_exist(initmaps.ebv,'calc') then begin
         for i=0,n_elements(initmaps.ebv.calc)-1 do begin

            fluxtype = initmaps.ebv.calc[i]
            map = ebv[fluxtype]
            maperr = errebv[fluxtype]
            ibd = where(map eq bad OR ~ finite(map),ctbd)
            igd = where(map ne bad AND finite(map),ctgd)

            if ctgd gt 0 then begin
               cmap = cshst_fov_rb
               xran = [min(cmap[igd])*0.95d,max(cmap[igd])*1.05d]
               yran = plotdat_ebv[fluxtype,0:1]

               cgps_open,initdat.mapdir+initdat.label+'color_v_ebv_'+fluxtype+'.eps',$
                         charsize=1,/encap,/inches,xs=7.5d,ys=7.5d,/qui,/nomatch
               cgplot,[0],/xsty,/ysty,/nodata,xran=xran,yran=yran,$
                      xtit=initmaps.hstbl.label+'-'+initmaps.hstrd.label,$
                      ytit='gas E(B-V)'
               cgoplot,cmap[igd],map[igd],psym=16,symsize=0.75d,$
                       err_ylow=maperr[igd],err_yhi=maperr[igd]
;              Same model as above but with steeper slope
               xmod = [min(cmap[igd]),max(cmap[igd])]
               m /= 0.44d
               b = -m*xint
               ymod = m*xmod+b
               cgoplot,xmod,ymod
               cgps_close
            endif
         endfor
      endif
   endif
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Plots of line ratios
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if tag_exist(initmaps,'lr') then begin
      if tag_exist(initmaps.lr,'calc') then begin

;        Loop through types
         for i=0,n_elements(initmaps.lr.calc)-1 do begin

            fluxtype = initmaps.lr.calc[i]
            plotorder = ['o3hb','n2ha','s2ha','o1ha']
            lrtype = (lr[fluxtype].keys())->toarray()

;           Size of plot grid
            lrsort=!NULL
            npx=lr[fluxtype].count()
            is2 = where(lr[fluxtype].keys() eq 's2',cts2)
            if cts2 gt 0 then npx--
            npy=1
            for j=0,n_elements(plotorder)-1 do begin
               ilr = where(lrtype eq plotorder[j],ctlr)
               if ctlr eq 1 then begin
                  lrsort = [lrsort,lrtype[ilr]]
                  if lrtype[ilr] eq 'o3hb' then npy=2
               endif
            endfor

;           Figure out correct image size in inches
            xpanel_in = 1.75d
            margin_in = 0.5d
            topmargin_in = 0.5d
            halfmargin_in = margin_in/1d
            xsize_in = xpanel_in*double(npx) + margin_in
            aspectrat_fov=double(dx)/double(dy)
            ysize_in = (xpanel_in * 1d/aspectrat_fov + margin_in) + $
                       (xpanel_in + 2d*margin_in) + topmargin_in
;           Sizes and positions of image windows in real and normalized coordinates
            pan_xfrac = xpanel_in/xsize_in
            pan_yfrac = xpanel_in/aspectrat_fov/ysize_in
            sqpan_yfrac = xpanel_in/ysize_in
            mar_xfrac = margin_in/xsize_in
            mar_yfrac = margin_in/ysize_in
            topmar_yfrac = topmargin_in/ysize_in
            hmar_xfrac = halfmargin_in/xsize_in
            hmar_yfrac = halfmargin_in/ysize_in
            pos_top = dblarr(4,npx)
            pos_toptit = dblarr(2,npx)
            pos_mid = dblarr(4,npx)
            pos_midtit = dblarr(2,npx)
            pos_bot = dblarr(4,npx)
            for k=0,npx-1 do begin
               pos_top[*,k] = [mar_xfrac+double(k)*pan_xfrac,$
                         1d - (mar_yfrac+pan_yfrac) - topmar_yfrac,$
                         mar_xfrac+double(k+1)*pan_xfrac,$
                         1d - mar_yfrac - topmar_yfrac]
               pos_bot[*,k] = [mar_xfrac+double(k)*pan_xfrac,$
                         mar_yfrac,$
                         mar_xfrac+double(k+1)*pan_xfrac,$
                         mar_yfrac+sqpan_yfrac]
            endfor
            pos_title1 = [mar_xfrac+double(npx)*pan_xfrac/2d,$
                          1d - topmar_yfrac*0.45d]
            pos_title2 = [mar_xfrac+double(npx)*pan_xfrac/2d,$
                          1d - topmar_yfrac*0.85d]

            cgps_open,initdat.mapdir+initdat.label+'lr_'+fluxtype+'.eps',$
                      charsize=1,/encap,/inches,xs=xsize_in,ys=ysize_in,/qui,/nomatch
            cbform = '(D0.2)'

;           Loop through plot panels
            for j=0,npx-1 do begin

               map = lr[fluxtype,lrsort[j]]
               maperrlo = lrerrlo[fluxtype,lrsort[j]]
               maperrhi = lrerrhi[fluxtype,lrsort[j]]
               ibd = where(map eq bad OR ~ finite(map),ctbd)
               igd = where(map ne bad AND finite(map),ctgd)

;              Set up range
               if hasrangefile then begin
                  ithisline = where(rangeline eq lrsort[j] AND $
                                    rangequant eq fluxtype,ctthisline)
                  if ctthisline eq 1 then auto=0b
               endif else auto=1b
               plotdat = $
                  ifsf_plotrange(auto=auto,$
                                 mapgd=map[igd],divinit=0.5d,$
                                 ncbdivmax=ncbdivmax,$
                                 rline=rangeline,matline=lrsort[j],$
                                 rcomp=rangecomp,$
                                 rquant=rangequant,matquant=fluxtype,$
                                 rncbdiv=rangencbdiv,$
                                 rlo=rangelo,rhi=rangehi)

               mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
                               min=plotdat[0],max=plotdat[1])
               cgloadct,74,/reverse
               cgimage,mapscl,/keep,pos=pos_top[*,j],opos=truepos,$
                       /noerase,missing_value=bad,missing_index=255,$
                       missing_color='white'
               cgplot,[0],xsty=5,ysty=5,position=truepos,$
                      /nodata,/noerase
               if j eq 0 then nolab=0b else nolab=1b
               ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,$
                                   center_nuclei_kpc_y,nolab=nolab
               xoffset = pan_xfrac*0.05
               yoffset = mar_yfrac*0.2
               cbpos=[truepos[0]+xoffset,truepos[3],truepos[2]-xoffset,$
                      truepos[3]+yoffset]
               ticknames = string(dindgen(plotdat[3]+1)*$
                                  plotdat[2]/double(plotdat[3]) - $
                                  (plotdat[2] - plotdat[1]),format=cbform)
               cgcolorbar,position=cbpos,divisions=plotdat[3],$
                          ticknames=ticknames,/top,charsize=0.6
               yoffset = mar_yfrac*0.65
               if lrsort[j] eq 'o3hb' then panel_title='[OIII]/H$\beta$'
               if lrsort[j] eq 'n2ha' then panel_title='[NII]/H$\alpha$'
               if lrsort[j] eq 'o1ha' then panel_title='[OI]/H$\alpha$'
               if lrsort[j] eq 's2ha' then panel_title='[SII]/H$\alpha$'
               cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
                      truepos[3]+yoffset,panel_title,$
                      charsize=1.05,align=0.5,/norm

               if lrsort[j] eq 'o3hb' then begin
                  map_o3hb = map
                  maperrlo_o3hb = maperrlo
                  maperrhi_o3hb = maperrhi
                  igd_o3hb = igd
                  plotdat_o3hb = plotdat
               endif

               if npy eq 2 AND j ne 0 then begin

                  if lrsort[j] eq 'n2ha' then begin
                     xran = [-1.99d,0.99d]
                     yran = [-1.19d,1.49d]
                     xkew1 = 0.05d*indgen(110)-5d
                     ykew1 = 0.61d / (xkew1-0.47d)+1.19d
                     xkew2 = 0.05d*indgen(41)-2d
                     ykew2 = 0.61d / (xkew2-0.05d)+1.3d
                  endif
                  if lrsort[j] eq 'o1ha' then begin
                     xran = [-2.19d,0.79d]
                     yran = [-1.19d,1.49d]
                     xkew1 = 0.05*indgen(85)-5
                     ykew1 = 0.73d / (xkew1+0.59d)+1.33d
                     xkew2 = 0.5d*indgen(2)-1.1d
                     ykew2 = 1.18d*xkew2 + 1.30d
                  endif
                  if lrsort[j] eq 's2ha' then begin
                     xran = [-1.99d,0.99d]
                     yran = [-1.19d,1.49d]
                     xkew1 = 0.05*indgen(105)-5
                     ykew1 = 0.72d / (xkew1-0.32d)+1.30d
                     xkew2 = 0.5d*indgen(2)-0.4d
                     ykew2 = 1.89d*xkew2+0.76d
                  endif

                  igdvo = cgsetintersection(igd,igd_o3hb)

                  if j eq 1 then ytit = '[OIII]/H$\beta$' else ytit=''
                  if j gt 1 then ytickf='(A1)' else ytickf='(D0.1)'
                  cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,pos=pos_bot[*,j],$
                             /nodata,/noerase,xtit=panel_title,ytit=ytit,$
                             aspect=1d,ytickf=ytickf
                  cgoplot,map[igdvo],map_o3hb[igdvo],psym=16,symsize=0.01,$
                          err_xlow=maperrlo[igdvo],err_xhigh=maperrhi[igdvo],$
                          err_ylow=maperrlo_o3hb[igdvo],$
                          err_yhigh=maperrhi_o3hb[igdvo],$
                          err_color='Red',err_thick=2,/err_clip,$
                          err_width=0d
                  cgoplot,map[igdvo],map_o3hb[igdvo],psym=16,symsize=0.5d
                  cgoplot,xkew1,ykew1
                  cgoplot,xkew2,ykew2,linesty=1

               endif

            endfor

;           Title
            cgplot,[0],xsty=5,ysty=5,position=[0,0,1,1],/nodata,/noerase
            cgtext,pos_title1[0],pos_title1[1],$
                   initdat.name+':',$
                   charsize=1.25d,align=0.5
            cgtext,pos_title2[0],pos_title2[1],$
                   'Line Ratios',$
                   charsize=1.25d,align=0.5

            cgps_close

         endfor
      endif
   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Emission line radial profiles
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if tag_exist(initmaps,'emlplot') then begin
      if tag_exist(initmaps.emlplot,'radtags') AND $
         tag_exist(initmaps.emlplot,'radtitles') then begin

;        Quantities to plot
         radtags = initmaps.emlplot.radtags

;        Size of plot grid
         npx = n_elements(radtags)
         npy = 1

;        Figure out correct image size in inches
         xpanel_in = 1.75d
         margin_in = 0.5d
         topmargin_in = 0.25d
         halfmargin_in = margin_in/1d
         xsize_in = xpanel_in*double(npx) + margin_in
         ysize_in = xpanel_in + margin_in*2d + topmargin_in
      ;  Sizes and positions of image windows in real and normalized coordinates
         pan_xfrac = xpanel_in/xsize_in
         pan_yfrac = xpanel_in/ysize_in
         mar_xfrac = margin_in/xsize_in
         mar_yfrac = margin_in/ysize_in
         topmar_yfrac = topmargin_in/ysize_in
         hmar_xfrac = halfmargin_in/xsize_in
         hmar_yfrac = halfmargin_in/ysize_in
         pos_top = dblarr(4,npx)
         pos_toptit = dblarr(2,npx)
         pos_mid = dblarr(4,npx)
         pos_midtit = dblarr(2,npx)
         pos_bot = dblarr(4,npx)
         for i=0,npx-1 do $
            pos_bot[*,i] = [mar_xfrac+double(i)*pan_xfrac,$
                            mar_yfrac,$
                            mar_xfrac+double(i+1)*pan_xfrac,$
                            mar_yfrac+pan_yfrac]
         pos_title = [mar_xfrac+double(npx)*pan_xfrac/2d,$
                      1d - topmar_yfrac*1.5d]
         pos_xtit = [mar_xfrac+double(npx)*pan_xfrac/2d,$
                     mar_yfrac*0.25d]

         yran = [-4d,0d]
         cbform = '(D0.1)' ; colorbar syntax

;        Loop through emission lines
         foreach line,lines_with_doublets do begin

            linelab = ifsf_linesyntax(line)
            cgps_open,initdat.mapdir+initdat.label+linelab+'_rad.eps',$
                      charsize=1,/encap,/inches,xs=xsize_in,ys=ysize_in,/qui,/nomatch

            ichan=0
            for j=0,npx-1 do begin

;              Get map and scale
               map = emlflx[radtags[j],line]
               ibd = where(map eq bad OR map eq 0 OR ~ finite(map),ctbd)
               igd = where(map ne bad AND map ne 0 AND finite(map),ctgd)

               if ctgd gt 0 then begin

                  zmax_flux = max(map[igd])
                  maplin = map
                  maplin[igd] = alog10(maplin[igd]/zmax_flux[0])
                  if ctbd gt 0 then maplin[ibd] = bad

                  title = initmaps.emlplot.radtitles[j]
                  if j eq 0 then $
                     cgplot,map_rkpc_ifs,maplin,/xsty,/ysty,yran=yran,$
                            xran=[0,max(map_rkpc_ifs)],/noerase,$
                            title=title,pos=pos_bot[*,j],psym=16,symsize=1d,$
                            aspect=1.0d,ytit='log(F/F$\downmax$)' $
                  else $
                     cgplot,map_rkpc_ifs,maplin,/xsty,/ysty,yran=yran,$
                            xran=[0,max(map_rkpc_ifs)],/noerase,$
                            title=title,pos=pos_bot[*,j],psym=16,symsize=1d,$
                            aspect=1.0d,xtickn=replicate(' ',60),$
                            ytickn=replicate(' ',60)
                  if tag_exist(initdat,'decompose_qso_fit') then begin
                     cgoplot,psf1d_x,psf1d_y,color='Red'
                  endif else if tag_exist(initmaps,'fit_empsf') then begin
                     cgoplot,empsf1d_x,empsf1d_y,color='Red'
                  endif else if $
                     tag_exist(initmaps.emlplot,'rad_psffwhm') then begin
                     x = dindgen(101)/100d*max(map_rkpc_ifs)
                     fwhm=initmaps.emlinradprof_psffwhm * kpc_per_as
;                    Gaussian
                     y = alog10(gaussian(x,[1d,0d,fwhm/2.35]))
                     cgoplot,x,y,color='Black'
; Moffat index values chosen to match turbulence theory (5), IRAF default (2.5),
; and wingy profile (1.5). These are the same chosen by Trujillo et al. 2001.
;                    Moffat, index = 1.5
                     y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/1.5d)-1),1.5d]))
                     cgoplot,x,y,color='Red',/linesty
;                    Moffat, index = 2.5
                     y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/2.5d)-1),2.5d]))
                     cgoplot,x,y,color='Red'
;                    Moffat, index = 5
                     y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/5d)-1),5d]))
                     cgoplot,x,y,color='Blue'
                  endif

               endif
            endfor

;           Title and axis labels
            cgplot,[0],xsty=5,ysty=5,position=[0,0,1,1],/nodata,/noerase
            cgtext,pos_title[0],pos_title[1],initdat.name+': '+linelabels[line],$
                   charsize=1d,align=0.5
            cgtext,pos_xtit[0],pos_xtit[1],'Radius (kpc)',$
                   charsize=1.25d,align=0.5
            cgps_close

         endforeach
      endif
   endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; end emission line plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; NaD PLOTS
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if tag_exist(initdat,'donad') then begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; W_eq + velocities [absorption] [presently only works for 1 or 2 comp]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      npx_fov=3
      npx_sq=1
      npx = npx_fov + npx_sq
;      npy=1 + maxnadabsncomp_act

;     Figure out correct image size in inches
      xpanel_in = 1.75d
      margin_in = 0.5d
      topmargin_in = 0.5d
      halfmargin_in = margin_in/1d
      xsize_in = xpanel_in*double(npx_fov+npx_sq)+margin_in*(1d +double(npx_sq))
      aspectrat_fov=double(dx)/double(dy)
      ysize_in = (xpanel_in/aspectrat_fov + margin_in)*2d + topmargin_in
;      if npy eq 3 then ysize_in += xpanel_in/aspectrat_fov
;     Sizes and positions of image windows in real and normalized coordinates
      pan_xfrac = xpanel_in/xsize_in
      pan_yfrac = xpanel_in/aspectrat_fov/ysize_in
      sqpan_yfrac = xpanel_in/ysize_in
      mar_xfrac = margin_in/xsize_in
      mar_yfrac = margin_in/ysize_in
      topmar_yfrac = topmargin_in/ysize_in
      hmar_xfrac = halfmargin_in/xsize_in
      hmar_yfrac = halfmargin_in/ysize_in
      pos_top = dblarr(4,npx)
      pos_toptit = dblarr(2,npx)
;      pos_mid = dblarr(4,npx)
;      pos_midtit = dblarr(2,npx)
      pos_bot = dblarr(4,npx)
      for i=0,npx-1 do begin
         pos_top[*,i] = [mar_xfrac+double(i)*pan_xfrac,$
                         1d - (mar_yfrac+pan_yfrac) - topmar_yfrac,$
                         mar_xfrac+double(i+1)*pan_xfrac,$
                         1d - mar_yfrac - topmar_yfrac]
         pos_bot[*,i] = [mar_xfrac+double(i)*pan_xfrac,$
                         mar_yfrac,$
                         mar_xfrac+double(i+1)*pan_xfrac,$
                         mar_yfrac+pan_yfrac]
;         if npy eq 3 then $
;            pos_mid[*,i] = pos_bot[*,i] + [0d,pan_yfrac,0d,pan_yfrac]
      endfor
      if dx gt dy then pad_yfrac = xpanel_in*(1d - 1d/aspectrat_fov)/ysize_in $
      else pad_yfrac = 0d
      pos_top[*,3] = [mar_xfrac+double(npx_fov)*pan_xfrac+mar_xfrac,$
                      1d - (mar_yfrac/2d +pan_yfrac) - topmar_yfrac + pad_yfrac,$
                      mar_xfrac+double(npx_fov)*pan_xfrac+mar_xfrac+pan_xfrac,$
                      1d - mar_yfrac/2d - topmar_yfrac + pad_yfrac]
      pos_title = [mar_xfrac+double(npx)*pan_xfrac/2d,$
                   1d - topmar_yfrac*0.75d]

      cgps_open,initdat.mapdir+initdat.label+'NaDabs.eps',$
                charsize=1,/encap,/inches,xs=xsize_in,ys=ysize_in,/qui,/nomatch
      ranlin='NaDabs'
      cbform = '(I0)'
      panel_title_chars=1.05

;     Empirical W_eq
;-------------------------------------------------------------------------------
      ranlab = 'empweq'
      panel_title = textoidl('W_{eq}^{abs} (emp)')
      cbdivinit=1d
      map = nadcube.weq[*,*,0]
      maperr = nadcube.weq[*,*,1]
      igd = where(map ge initmaps.nadabsweq_snrthresh*maperr AND $
                  map gt 0d AND map ne bad)
      ibd = where(map lt initmaps.nadabsweq_snrthresh*maperr OR $
                  map eq 0d OR map eq bad)

      if hasrangefile then begin
         ithisline = where(rangeline eq ranlin AND $
                           rangequant eq ranlab,ctthisline)
         if ctthisline eq 1 then auto=0b
      endif else auto=1b
      plotdat = $
         ifsf_plotrange(auto=auto,$
                        mapgd=map[igd],divinit=cbdivinit,ncbdivmax=ncbdivmax,$
                        rline=rangeline,matline=ranlin,$
                        rquant=rangequant,matquant=ranlab,$
                        rncbdiv=rangencbdiv,rlo=rangelo,rhi=rangehi)

      map[ibd] = bad

;     Save some things for later
      igd_nadabs_empweq = igd
      ibd_nadabs_empweq = ibd
      plotdat_nadabs_empweq = plotdat
      map_nadabs_empweq = map
      map_nadabs_empweq_err = maperr

      mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
                      min=plotdat[0],max=plotdat[1])
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos_top[*,0],opos=truepos,$
              missing_value=bad,missing_index=255,$
              missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
             /nodata,/noerase
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,$
                       center_nuclei_kpc_x,center_nuclei_kpc_y,/noxlab
      xoffset = pan_xfrac*0.05
      yoffset = mar_yfrac*0.2
      cbpos=[truepos[0]+xoffset,truepos[3],truepos[2]-xoffset,$
             truepos[3]+yoffset]
      ticknames = string(dindgen(plotdat[3]+1)*$
                         plotdat[2]/double(plotdat[3]) - $
                         (plotdat[2] - plotdat[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=plotdat[3],$
                 ticknames=ticknames,/top,charsize=0.6
      yoffset = mar_yfrac*0.65
      cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
             truepos[3]+yoffset,panel_title,$
             charsize=panel_title_chars,align=0.5,/norm

;     Fitted W_eq
;-------------------------------------------------------------------------------
      ranlab = 'fitweq'
      panel_title = textoidl('W_{eq}^{abs} (fit)')
      cbdivinit=1d

      map = nadfit.weqabs[*,*,0]
      maperrlo = nadfit.weqabserr[*,*,0]
      maperrhi = nadfit.weqabserr[*,*,1]
      maperravg = (maperrlo + maperrhi) / 2d
      igd = where(map gt 0d AND $
                  map ne bad AND $
                  map ge initmaps.nadabsweq_snrthresh*maperravg,ctgdtmp)
      ibd = where(map eq 0d OR $
                  map eq bad OR $
                  map lt initmaps.nadabsweq_snrthresh*maperravg)
      ilowsnr = where(map ne bad AND $
                      map lt initmaps.nadabsweq_snrthresh*maperravg,ctgdtmp)

      if ctgdtmp gt 0 then begin

         if hasrangefile then begin
            ithisline = where(rangeline eq ranlin AND $
                              rangequant eq ranlab,ctthisline)
            if ctthisline eq 1 then auto=0b
         endif else auto=1b
         plotdat = $
            ifsf_plotrange(auto=auto,$
                           mapgd=map[igd],divinit=cbdivinit,ncbdivmax=ncbdivmax,$
                           rline=rangeline,matline=ranlin,$
                           rquant=rangequant,matquant=ranlab,$
                           rncbdiv=rangencbdiv,rlo=rangelo,rhi=rangehi)

         map[ibd] = bad

;        Save some things for later
         plotdat_nadabs_fitweq = plotdat

         mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
                         min=plotdat[0],max=plotdat[1])
         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos_top[*,1],opos=truepos,$
                 missing_value=bad,missing_index=255,$
                 missing_color='white',/noerase
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,$
                          center_nuclei_kpc_x,center_nuclei_kpc_y,/nolab
         xoffset = pan_xfrac*0.05
         yoffset = mar_yfrac*0.2
         cbpos=[truepos[0]+xoffset,truepos[3],truepos[2]-xoffset,$
                truepos[3]+yoffset]
         ticknames = string(dindgen(plotdat[3]+1)*$
                            plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,/top,charsize=0.6
         yoffset = mar_yfrac*0.65
         cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
                truepos[3]+yoffset,panel_title,$
                charsize=panel_title_chars,align=0.5,/norm
      endif

;     Save some things for later
      igd_nadabs_fitweq = igd
      ilowsnr_nadabs_fitweq = ilowsnr
      ibd_nadabs_fitweq = ibd
      map_nadabs_fitweq = map
      map_nadabs_fitweq_errlo = maperrlo
      map_nadabs_fitweq_errhi = maperrhi
      

;     S/N (fitted)
;-------------------------------------------------------------------------------
;      ranlab = 'empweqsnr'
;      panel_title = textoidl('W_{eq}^{abs}/\deltaW')
;      cbdivinit=5d
;
;      map = map_nadabs_empweq
;      maperr = map_nadabs_empweq_err
;      igd = igd_nadabs_empweq
;      map[igd] /= maperr[igd]

      ranlab = 'fitweqsnr'
      panel_title = textoidl('W_{eq}^{abs}/\deltaW')
      cbdivinit=5d

      map = map_nadabs_fitweq
      maperr = sqrt((map_nadabs_fitweq_errlo^2d + map_nadabs_fitweq_errhi^2d)/2d)
      igd = igd_nadabs_fitweq
      map[igd] /= maperr[igd]

      if hasrangefile then begin
         ithisline = where(rangeline eq ranlin AND $
                           rangequant eq ranlab,ctthisline)
         if ctthisline eq 1 then auto=0b
      endif else auto=1b
      plotdat = $
         ifsf_plotrange(auto=auto,$
                        mapgd=map[igd],divinit=cbdivinit,ncbdivmax=ncbdivmax,$
                        rline=rangeline,matline=ranlin,$
                        rquant=rangequant,matquant=ranlab,$
                        rncbdiv=rangencbdiv,rlo=rangelo,rhi=rangehi)


      mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
                      min=plotdat[0],max=plotdat[1])
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos_top[*,2],opos=truepos,$
              missing_value=bad,missing_index=255,$
              missing_color='white',/noerase
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
             /nodata,/noerase
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,$
                       center_nuclei_kpc_x,center_nuclei_kpc_y,/nolab
      xoffset = pan_xfrac*0.05
      yoffset = mar_yfrac*0.2
      cbpos=[truepos[0]+xoffset,truepos[3],truepos[2]-xoffset,$
             truepos[3]+yoffset]
      ticknames = string(dindgen(plotdat[3]+1)*$
                         plotdat[2]/double(plotdat[3]) - $
                         (plotdat[2] - plotdat[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=plotdat[3],$
                 ticknames=ticknames,/top,charsize=0.6
      yoffset = mar_yfrac*0.65
      cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
             truepos[3]+yoffset,panel_title,$
             charsize=panel_title_chars,align=0.5,/norm

;     Fitted W_eq vs. Empirical W_eq
;-------------------------------------------------------------------------------
      igderr = where(map_nadabs_fitweq ne 0d AND map_nadabs_fitweq ne bad AND $
                     map_nadabs_empweq ne 0d AND map_nadabs_empweq ne bad,ctgderr)

      if ctgderr gt 0 then begin
         ibderr = where(map_nadabs_fitweq eq 0d OR map_nadabs_fitweq eq bad OR $
                        map_nadabs_empweq eq 0d OR map_nadabs_empweq eq bad)
         map_nadabs_fitweq_errlo[ibderr] = 0d
         map_nadabs_fitweq_errhi[ibderr] = 0d
         map_nadabs_empweq_err[ibderr] = 0d

         xran = [0d,max([map_nadabs_fitweq[igd_nadabs_fitweq],$
                         map_nadabs_empweq[igd_nadabs_empweq]])*1.1d]
         yran = xran
         cgplot,map_nadabs_empweq,map_nadabs_fitweq,/xsty,/ysty,$
                xran=xran,yran=yran,psym=3,pos=pos_top[*,3],$
                xtit=textoidl('W_{eq}^{abs} (emp)'),$
                ytit=textoidl('W_{eq}^{abs} (fit)'),$
                /noerase,aspect=1d,$
                chars=0.8,err_width=0,err_color='Gray',/err_clip,$
                err_xlow=map_nadabs_empweq_err,err_xhigh=map_nadabs_empweq_err,$
                err_ylow=map_nadabs_fitweq_errlo,err_yhigh=map_nadabs_fitweq_errhi
         cgoplot,map_nadabs_empweq,map_nadabs_fitweq,psym=16,symsize=0.4
         cgoplot,xran,xran
         inz = where(map_nadabs_fitweq ne 0d AND map_nadabs_fitweq ne bad AND $
                     map_nadabs_empweq ne 0d AND map_nadabs_empweq ne bad)
         weq_rms = sqrt(mean((map_nadabs_fitweq[inz]-map_nadabs_empweq[inz])^2d))
         cgoplot,xran-weq_rms,xran,linesty=3
         cgoplot,xran+weq_rms,xran,linesty=3
      endif

;     Velocities
;-------------------------------------------------------------------------------
;      abslab = ['Abs (1 Comp)']
;      if maxnadabsncomp_act gt 1 then $
;         for i=1,nyabs-1 do $
;            abslab = [abslab,string('Abs (',i,' of ',nyabs-1,' Comp)',$
;                                    format=('(A0,I0,A0,I0,A0)'))]

;      for i=0,maxnadabsncomp_act-1 do begin
;        
;         if i eq 0 then pos_use = pos_bot
;         if i eq 1 then pos_use = pos_mid
        
;        sigma
;-------------------------------------------------------------------------------
         panel_title = '$\sigma$ (km/s)'
;         ranlab = 'fitc'+string(i+1,format='(I0)')+'vsig'
         ranlab = 'fitvsig'
         cbdivinit=100d
         
;         map = nadabssig[*,*,i]
         map = nadabscvdfvals['vsig']
         igd_thiscomp = where(map ne bad,ctgd_thiscomp)
         ibd_thiscomp = where(map eq bad,ctbd_thiscomp)
         if ctgd_thiscomp gt 0 then $
            igd = cgsetintersection(igd_nadabs_fitweq,igd_thiscomp) $
         else igd=igd_nadabs_fitweq
         if ctbd_thiscomp gt 0 then $
            ibd = cgsetunion(ibd_nadabs_fitweq,ibd_thiscomp) $
         else ibd = ibd_nadabs_fitweq

         map[ibd] = bad

         if hasrangefile then begin
            ithisline = where(rangeline eq ranlin AND $
                              rangequant eq ranlab,ctthisline)
            if ctthisline eq 1 then auto=0b
         endif else auto=1b
         plotdat = $
            ifsf_plotrange(auto=auto,$
                           mapgd=map[igd],divinit=cbdivinit,ncbdivmax=ncbdivmax,$
                           rline=rangeline,matline=ranlin,$
                           rquant=rangequant,matquant=ranlab,$
                           rncbdiv=rangencbdiv,rlo=rangelo,rhi=rangehi)

         mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
                         min=plotdat[0],max=plotdat[1])
         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos_bot[*,0],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,$
                          center_nuclei_kpc_y,/noxlab
         xoffset = pan_xfrac*0.05
         yoffset = mar_yfrac*0.2
         cbpos=[truepos[0]+xoffset,truepos[1]-yoffset,$
                truepos[2]-xoffset,truepos[1]]
         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         ticknames[0]=' '
         ticknames[plotdat[3]]=' '
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,charsize=0.6
         yoffset = mar_yfrac*0.85
         cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
                truepos[1]-yoffset,panel_title,$
                charsize=panel_title_chars,align=0.5,/norm
         
;        vpk
;-------------------------------------------------------------------------------
         panel_title = 'v$\downpk$ (km/s)'
;         ranlab = 'fitc'+string(i+1,format='(I0)')+'vpk'
         ranlab = 'fitvpk'
         cbdivinit=100d
         
;         map = nadabsvel[*,*,i]
         map = nadabscvdfvals['vpk']
         igd_thiscomp = where(map ne bad,ctgd_thiscomp)
         ibd_thiscomp = where(map eq bad,ctbd_thiscomp)
         if ctgd_thiscomp gt 0 then $
            igd = cgsetintersection(igd_nadabs_fitweq,igd_thiscomp) $
         else igd=igd_nadabs_fitweq
         if ctbd_thiscomp gt 0 then $
            ibd = cgsetunion(ibd_nadabs_fitweq,ibd_thiscomp) $
         else ibd = ibd_nadabs_fitweq

         map[ibd] = bad

         if hasrangefile then begin
            ithisline = where(rangeline eq ranlin AND $
                              rangequant eq ranlab,ctthisline)
            if ctthisline eq 1 then auto=0b
         endif else auto=1b
         plotdat = $
            ifsf_plotrange(auto=auto,$
                           mapgd=map[igd],divinit=cbdivinit,ncbdivmax=ncbdivmax,$
                           rline=rangeline,matline=ranlin,$
                           rquant=rangequant,matquant=ranlab,$
                           rncbdiv=rangencbdiv,rlo=rangelo,rhi=rangehi)

         mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
                         min=plotdat[0],max=plotdat[1])
         cgloadct,74,/reverse
         cgimage,mapscl,/keep,pos=pos_bot[*,1],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase,xran=[0,dx],yran=[0,dy]
         maptmp = map*0d
         iblue = where(map lt 0 AND map ne bad,ctblue)
         if ctblue gt 0 then maptmp[iblue] = 1d
         cgcontour,maptmp,dindgen(dx)+0.5,dindgen(dy)+0.5,$
                   /overplot,color=0,c_thick=4,$
                   levels=[1d]
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,$
                          center_nuclei_kpc_y,/nolab
         xoffset = pan_xfrac*0.05
         yoffset = mar_yfrac*0.2
         cbpos=[truepos[0]+xoffset,truepos[1]-yoffset,$
                truepos[2]-xoffset,truepos[1]]
         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         ticknames[0]=' '
         ticknames[plotdat[3]]=' '
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,charsize=0.6
         yoffset = mar_yfrac*0.85
         cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
                truepos[1]-yoffset,panel_title,$
                charsize=panel_title_chars,align=0.5,/norm

         
;        v50
;-------------------------------------------------------------------------------
         panel_title = 'v$\down50$ (km/s)'
;         ranlab = 'fitc'+string(i+1,format='(I0)')+'v%50'
         ranlab = 'fitv%50'
         cbdivinit=100d
         
;         map = nadabsvel[*,*,i]
         map = nadabscvdfvals['v%50']
         igd_thiscomp = where(map ne bad,ctgd_thiscomp)
         ibd_thiscomp = where(map eq bad,ctbd_thiscomp)
         if ctgd_thiscomp gt 0 then $
            igd = cgsetintersection(igd_nadabs_fitweq,igd_thiscomp) $
         else igd=igd_nadabs_fitweq
         if ctbd_thiscomp gt 0 then $
            ibd = cgsetunion(ibd_nadabs_fitweq,ibd_thiscomp) $
         else ibd = ibd_nadabs_fitweq

         map[ibd] = bad

         if hasrangefile then begin
            ithisline = where(rangeline eq ranlin AND $
                              rangequant eq ranlab,ctthisline)
            if ctthisline eq 1 then auto=0b
         endif else auto=1b
         plotdat = $
            ifsf_plotrange(auto=auto,$
                           mapgd=map[igd],divinit=cbdivinit,ncbdivmax=ncbdivmax,$
                           rline=rangeline,matline=ranlin,$
                           rquant=rangequant,matquant=ranlab,$
                           rncbdiv=rangencbdiv,rlo=rangelo,rhi=rangehi)

         mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
                         min=plotdat[0],max=plotdat[1])
         cgloadct,74,/reverse
         cgimage,mapscl,/keep,pos=pos_bot[*,2],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase,xran=[0,dx],yran=[0,dy]
         maptmp = map*0d
         iblue = where(map lt 0 AND map ne bad,ctblue)
         if ctblue gt 0 then maptmp[iblue] = 1d
         cgcontour,maptmp,dindgen(dx)+0.5,dindgen(dy)+0.5,$
                   /overplot,color=0,c_thick=4,$
                   levels=[1d]
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,$
                          center_nuclei_kpc_y,/nolab
         xoffset = pan_xfrac*0.05
         yoffset = mar_yfrac*0.2
         cbpos=[truepos[0]+xoffset,truepos[1]-yoffset,$
                truepos[2]-xoffset,truepos[1]]
         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         ticknames[0]=' '
         ticknames[plotdat[3]]=' '
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,charsize=0.6
         yoffset = mar_yfrac*0.85
         cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
                truepos[1]-yoffset,panel_title,$
                charsize=panel_title_chars,align=0.5,/norm
         
;        v98
;-------------------------------------------------------------------------------
         panel_title = 'v$\down98$ (km/s)'
;         ranlab = 'fitc'+string(i+1,format='(I0)')+'v%98'
         ranlab = 'fitv%98'
         cbdivinit=200d
         
;         map = nadabsv98[*,*,i]
         map = nadabscvdfvals['v%98']
         igd_thiscomp = where(map ne bad,ctgd_thiscomp)
         ibd_thiscomp = where(map eq bad,ctbd_thiscomp)
         if ctgd_thiscomp gt 0 then $
            igd = cgsetintersection(igd_nadabs_fitweq,igd_thiscomp) $
         else igd=igd_nadabs_fitweq
         if ctbd_thiscomp gt 0 then $
            ibd = cgsetunion(ibd_nadabs_fitweq,ibd_thiscomp) $
         else ibd = ibd_nadabs_fitweq

         map[ibd] = bad

         if hasrangefile then begin
            ithisline = where(rangeline eq ranlin AND $
                              rangequant eq ranlab,ctthisline)
            if ctthisline eq 1 then auto=0b
         endif else auto=1b
         plotdat = $
            ifsf_plotrange(auto=auto,$
                           mapgd=map[igd],divinit=cbdivinit,ncbdivmax=ncbdivmax,$
                           rline=rangeline,matline=ranlin,$
                           rquant=rangequant,matquant=ranlab,$
                           rncbdiv=rangencbdiv,rlo=rangelo,rhi=rangehi)

         mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
                         min=plotdat[0],max=plotdat[1])
         cgloadct,74,/reverse
         cgimage,mapscl,/keep,pos=pos_bot[*,3],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,$
                          center_nuclei_kpc_y,/nolab
         xoffset = pan_xfrac*0.05
         yoffset = mar_yfrac*0.2
         cbpos=[truepos[0]+xoffset,truepos[1]-yoffset,$
                truepos[2]-xoffset,truepos[1]]
         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         ticknames[0]=' '
         ticknames[plotdat[3]]=' '
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,charsize=0.6
         yoffset = mar_yfrac*0.85
         cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
                truepos[1]-yoffset,panel_title,$
                charsize=panel_title_chars,align=0.5,/norm
;
;      endfor




;     Title
;-------------------------------------------------------------------------------
      cgplot,[0],xsty=5,ysty=5,position=[0,0,1,1],/nodata,/noerase
      cgtext,pos_title[0],pos_title[1],$
             initdat.name+': NaI D absorption',$
             charsize=1.25d,align=0.5

      cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; W_eq + velocities [emission] [presently only works for 1 or 2 comp]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      npx_fov=3
      npx_sq=1
      npx = npx_fov + npx_sq
;      npy=1 + maxnadabsncomp_act

;     Figure out correct image size in inches
      xpanel_in = 1.75d
      margin_in = 0.5d
      topmargin_in = 0.5d
      halfmargin_in = margin_in/1d
      xsize_in = xpanel_in*double(npx_fov+npx_sq)+margin_in*(1d +double(npx_sq))
      aspectrat_fov=double(dx)/double(dy)
      ysize_in = (xpanel_in/aspectrat_fov + margin_in)*3d + topmargin_in
;     Sizes and positions of image windows in real and normalized coordinates
      pan_xfrac = xpanel_in/xsize_in
      pan_yfrac = xpanel_in/aspectrat_fov/ysize_in
      sqpan_yfrac = xpanel_in/ysize_in
      mar_xfrac = margin_in/xsize_in
      mar_yfrac = margin_in/ysize_in
      topmar_yfrac = topmargin_in/ysize_in
      hmar_xfrac = halfmargin_in/xsize_in
      hmar_yfrac = halfmargin_in/ysize_in
      pos_top = dblarr(4,npx)
      pos_toptit = dblarr(2,npx)
      pos_mid = dblarr(4,npx)
      pos_midtit = dblarr(2,npx)
      pos_bot = dblarr(4,npx)
      for i=0,npx-1 do begin
         pos_top[*,i] = [mar_xfrac+double(i)*pan_xfrac,$
                         1d - (mar_yfrac+pan_yfrac) - topmar_yfrac,$
                         mar_xfrac+double(i+1)*pan_xfrac,$
                         1d - mar_yfrac - topmar_yfrac]
         pos_mid[*,i] = [mar_xfrac+double(i)*pan_xfrac,$
                         mar_yfrac*2d + pan_yfrac,$
                         mar_xfrac+double(i+1)*pan_xfrac,$
                         (mar_yfrac+pan_yfrac)*2d]
         pos_bot[*,i] = [mar_xfrac+double(i)*pan_xfrac,$
                         mar_yfrac,$
                         mar_xfrac+double(i+1)*pan_xfrac,$
                         mar_yfrac+pan_yfrac]
      endfor
      if dx gt dy then pad_yfrac = xpanel_in*(1d - 1d/aspectrat_fov)/ysize_in $
      else pad_yfrac = 0d
      pos_top[*,3] = [mar_xfrac+double(npx_fov)*pan_xfrac+mar_xfrac,$
                      1d - (mar_yfrac/2d +pan_yfrac) - topmar_yfrac + pad_yfrac,$
                      mar_xfrac+double(npx_fov)*pan_xfrac+mar_xfrac+pan_xfrac,$
                      1d - mar_yfrac/2d - topmar_yfrac + pad_yfrac]
      pos_mid[*,3] = [mar_xfrac+double(npx_fov)*pan_xfrac+mar_xfrac,$
                      1d - (mar_yfrac/2d +2d*pan_yfrac) - topmar_yfrac + pad_yfrac,$
                      mar_xfrac+double(npx_fov)*pan_xfrac+mar_xfrac+pan_xfrac,$
                      1d - mar_yfrac/2d - topmar_yfrac + pad_yfrac - pan_yfrac]
      pos_title = [mar_xfrac+double(npx)*pan_xfrac/2d,$
                   1d - topmar_yfrac*0.75d]

      cgps_open,initdat.mapdir+initdat.label+'NaDem.eps',$
                charsize=1,/encap,/inches,xs=xsize_in,ys=ysize_in,/qui,/nomatch
      ranlin='NaDem'
      cbform = '(I0)'
      panel_title_chars=1.05

;     Empirical W_eq
;-------------------------------------------------------------------------------
      ranlab = 'empweq'
      panel_title = textoidl('W_{eq}^{em} (emp)')
      cbdivinit=1d
      map = nadcube.weq[*,*,2]
      maperr = nadcube.weq[*,*,3]
      igd = where(abs(map) ge abs(initmaps.nademweq_snrthresh*maperr) AND $
                  map lt 0d,ctgdtmp_emp)
      ibd = where(abs(map) lt abs(initmaps.nademweq_snrthresh*maperr) OR $
                  map ge 0d,ctbdtmp)

      if ctgdtmp_emp gt 0 then begin

      if hasrangefile then begin
         ithisline = where(rangeline eq ranlin AND $
                           rangequant eq ranlab,ctthisline)
         if ctthisline eq 1 then auto=0b
      endif else auto=1b
      plotdat = $
         ifsf_plotrange(auto=auto,$
                        mapgd=map[igd],divinit=cbdivinit,ncbdivmax=ncbdivmax,$
                        rline=rangeline,matline=ranlin,$
                        rquant=rangequant,matquant=ranlab,$
                        rncbdiv=rangencbdiv,rlo=rangelo,rhi=rangehi)

      map[ibd] = bad

;     Save some things for later
      plotdat_nadem_empweq = plotdat

      mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
                      min=plotdat[0],max=plotdat[1])
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos_top[*,0],opos=truepos,$
              missing_value=bad,missing_index=255,$
              missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
             /nodata,/noerase
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,$
                       center_nuclei_kpc_x,center_nuclei_kpc_y,/noxlab
      xoffset = pan_xfrac*0.05
      yoffset = mar_yfrac*0.2
      cbpos=[truepos[0]+xoffset,truepos[3],truepos[2]-xoffset,$
             truepos[3]+yoffset]
      ticknames = string(dindgen(plotdat[3]+1)*$
                         plotdat[2]/double(plotdat[3]) - $
                         (plotdat[2] - plotdat[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=plotdat[3],$
                 ticknames=ticknames,/top,charsize=0.6
      yoffset = mar_yfrac*0.65
      cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
             truepos[3]+yoffset,panel_title,$
             charsize=panel_title_chars,align=0.5,/norm

      endif 
      
;     Save some things for later
      igd_nadem_empweq = igd
      ibd_nadem_empweq = ibd
      map_nadem_empweq = map
      map_nadem_empweq_err = maperr
      

;     Fitted W_eq
;-------------------------------------------------------------------------------
      ranlab = 'fitweq'
      panel_title = textoidl('W_{eq}^{em} (fit)')
      cbdivinit=1d

      map = nadfit.weqem[*,*,0]
      maperrlo = nadfit.weqemerr[*,*,0]
      maperrhi = nadfit.weqemerr[*,*,1]
      maperravg = (maperrlo + maperrhi) / 2d
      igd = where(map lt 0d AND $
                  maperrhi gt 0d AND $
                  abs(map) ge abs(initmaps.nademweq_snrthresh*maperravg),ctgdtmp_fit)
      ilowsnr = where(map le 0d AND $
                      abs(map) lt abs(initmaps.nademweq_snrthresh*maperravg),ctgdtmp_fit)
      ibd = where(map ge 0d OR $
                  abs(map) lt abs(initmaps.nademweq_snrthresh*maperravg))

      if ctgdtmp_fit gt 0 then begin


         if hasrangefile then begin
            ithisline = where(rangeline eq ranlin AND $
                              rangequant eq ranlab,ctthisline)
            if ctthisline eq 1 then auto=0b
         endif else auto=1b
         plotdat = $
            ifsf_plotrange(auto=auto,$
                           mapgd=map[igd],divinit=cbdivinit,ncbdivmax=ncbdivmax,$
                           rline=rangeline,matline=ranlin,$
                           rquant=rangequant,matquant=ranlab,$
                           rncbdiv=rangencbdiv,rlo=rangelo,rhi=rangehi)

         map[ibd] = bad
         maperrlo[ibd] = 0d
         maperrhi[ibd] = 0d

;        Save some things for later
         plotdat_nadem_fitweq = plotdat

         mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
                         min=plotdat[0],max=plotdat[1])
         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos_top[*,1],opos=truepos,$
                 missing_value=bad,missing_index=255,$
                 missing_color='white',/noerase
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,$
                          center_nuclei_kpc_x,center_nuclei_kpc_y,/nolab
         xoffset = pan_xfrac*0.05
         yoffset = mar_yfrac*0.2
         cbpos=[truepos[0]+xoffset,truepos[3],truepos[2]-xoffset,$
                truepos[3]+yoffset]
         ticknames = string(dindgen(plotdat[3]+1)*$
                            plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,/top,charsize=0.6
         yoffset = mar_yfrac*0.65
         cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
                truepos[3]+yoffset,panel_title,$
                charsize=panel_title_chars,align=0.5,/norm

;     S/N (fitted)
;-------------------------------------------------------------------------------

      ranlab = 'fitweqsnr'
      panel_title = textoidl('W_{eq}^{em}/\deltaW')
      cbdivinit=5d

      mapsnr = map
      maperr = sqrt((maperrlo^2d + maperrhi^2d)/2d)
      mapsnr[igd] /= maperr[igd]
      mapsnr[igd] = abs(mapsnr[igd])

      if hasrangefile then begin
         ithisline = where(rangeline eq ranlin AND $
                           rangequant eq ranlab,ctthisline)
         if ctthisline eq 1 then auto=0b
      endif else auto=1b
      plotdat = $
         ifsf_plotrange(auto=auto,$
                        mapgd=mapsnr[igd],divinit=cbdivinit,ncbdivmax=ncbdivmax,$
                        rline=rangeline,matline=ranlin,$
                        rquant=rangequant,matquant=ranlab,$
                        rncbdiv=rangencbdiv,rlo=rangelo,rhi=rangehi)


      mapscl = bytscl(rebin(mapsnr,dx*samplefac,dy*samplefac,/sample),$
                      min=plotdat[0],max=plotdat[1])
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos_top[*,2],opos=truepos,$
              missing_value=bad,missing_index=255,$
              missing_color='white',/noerase
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
             /nodata,/noerase
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,$
                       center_nuclei_kpc_x,center_nuclei_kpc_y,/nolab
      xoffset = pan_xfrac*0.05
      yoffset = mar_yfrac*0.2
      cbpos=[truepos[0]+xoffset,truepos[3],truepos[2]-xoffset,$
             truepos[3]+yoffset]
      ticknames = string(dindgen(plotdat[3]+1)*$
                         plotdat[2]/double(plotdat[3]) - $
                         (plotdat[2] - plotdat[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=plotdat[3],$
                 ticknames=ticknames,/top,charsize=0.6
      yoffset = mar_yfrac*0.65
      cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
             truepos[3]+yoffset,panel_title,$
             charsize=panel_title_chars,align=0.5,/norm

;     Fitted W_eq vs. Empirical W_eq
;-------------------------------------------------------------------------------
      igderr = where(map ne 0d AND map ne bad AND $
                     map_nadem_empweq ne 0d AND map_nadem_empweq ne bad,ctgderr)

      if ctgderr gt 0 then begin
         ibderr = where(map eq 0d OR map eq bad OR $
                        map_nadem_empweq eq 0d OR map_nadem_empweq eq bad)
         maperrlo[ibderr] = 0d
         maperrhi[ibderr] = 0d
         map_nadem_empweq_err[ibderr] = 0d

         xran = [min([map[igd],$
                      map_nadem_empweq[igd_nadem_empweq]])*1.1d,0d]
         yran = xran
         cgplot,map_nadem_empweq,map,/xsty,/ysty,$
                xran=xran,yran=yran,psym=3,pos=pos_top[*,3],$
                xtit=textoidl('W_{eq}^{em} (emp)'),$
                ytit=textoidl('W_{eq}^{em} (fit)'),$
                /noerase,aspect=1d,$
                chars=0.8,err_width=0,err_color='Gray',/err_clip,$
                err_xlow=map_nadem_empweq_err,err_xhigh=map_nadem_empweq_err,$
                err_ylow=maperrlo,err_yhigh=maperrhi
         cgoplot,map_nadem_empweq,map,psym=16,symsize=0.4
         cgoplot,xran,xran
         inz = where(map ne 0d AND map ne bad AND $
                     map_nadem_empweq ne 0d AND map_nadem_empweq ne bad)
         weq_rms = sqrt(mean((map[inz]-map_nadem_empweq[inz])^2d))
         cgoplot,xran-weq_rms,xran,linesty=3
         cgoplot,xran+weq_rms,xran,linesty=3
      endif

   endif

;  Save some things for later
   igd_nadem_fitweq = igd
   ilowsnr_nadem_fitweq = ilowsnr
   ibd_nadem_fitweq = ibd
   map_nadem_fitweq = map
   map_nadem_fitweq_errlo = maperrlo
   map_nadem_fitweq_errhi = maperrhi
   

;     Empirical Flux
;-------------------------------------------------------------------------------
      ranlab = 'empflux'
      panel_title = textoidl('I^{em} (emp)')
      cbdivinit=1d

      if ctgdtmp_emp gt 0 then begin

      map = nadcube.emflux[*,*,0]
      maperr = nadcube.emflux[*,*,1]
      igd = igd_nadem_empweq
      ibd = ibd_nadem_empweq
      if tag_exist(initmaps,'fluxunits') then begin
         map[igd] *= initmaps.fluxunits
         maperr[igd] *= initmaps.fluxunits
      endif
      if tag_exist(initmaps,'vornorm') then begin
         map[igd] /= initmaps.vornorm[igd]
         maperr[igd] /= initmaps.vornorm[igd]
      endif


      if hasrangefile then begin
         ithisline = where(rangeline eq ranlin AND $
                           rangequant eq ranlab,ctthisline)
         if ctthisline eq 1 then auto=0b
      endif else auto=1b
      plotdat = $
         ifsf_plotrange(auto=auto,$
                        mapgd=map[igd],divinit=cbdivinit,ncbdivmax=ncbdivmax,$
                        rline=rangeline,matline=ranlin,$
                        rquant=rangequant,matquant=ranlab,$
                        rncbdiv=rangencbdiv,rlo=rangelo,rhi=rangehi)

      map[ibd] = bad

;     Save some things for later
      plotdat_nadem_empflux = plotdat

      mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
                      min=plotdat[0],max=plotdat[1])
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos_mid[*,0],opos=truepos,$
              missing_value=bad,missing_index=255,$
              missing_color='white',/noerase
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
             /nodata,/noerase
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,$
                       center_nuclei_kpc_x,center_nuclei_kpc_y,/noxlab
      xoffset = pan_xfrac*0.05
      yoffset = mar_yfrac*0.2
      cbpos=[truepos[0]+xoffset,truepos[1]-yoffset,$
             truepos[2]-xoffset,truepos[1]]
      ticknames = string(dindgen(plotdat[3]+1)*$
                         plotdat[2]/double(plotdat[3]) - $
                         (plotdat[2] - plotdat[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=plotdat[3],$
                 ticknames=ticknames,charsize=0.6
      yoffset = mar_yfrac*0.85
      cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
             truepos[1]-yoffset,panel_title,$
             charsize=panel_title_chars,align=0.5,/norm

      endif

;     Save some things for later
      igd_nadem_empflux = igd
      ibd_nadem_empflux = ibd
      map_nadem_empflux = map
      map_nadem_empflux_err = maperr

;     Fitted Flux
;-------------------------------------------------------------------------------
      ranlab = 'fitflux'
      panel_title = textoidl('I^{em} (fit)')
      cbdivinit=1d

      map = nadfit.totfluxem[*,*,0]
      maperrlo = nadfit.totfluxemerr[*,*,0]
      maperrhi = nadfit.totfluxemerr[*,*,1]
      maperravg = (maperrlo + maperrhi) / 2d

      if ctgdtmp_fit gt 0 then begin

         igd = igd_nadem_fitweq
         ibd = ibd_nadem_fitweq
         if tag_exist(initmaps,'fluxunits') then begin
            map[igd] *= initmaps.fluxunits
            maperrlo[igd] *= initmaps.fluxunits
            maperrhi[igd] *= initmaps.fluxunits
            maperravg[igd] *= initmaps.fluxunits
         endif
         if tag_exist(initmaps,'vornorm') then begin
            map[igd] /= initmaps.vornorm[igd]
            maperrlo[igd] /= initmaps.vornorm[igd]
            maperrhi[igd] /= initmaps.vornorm[igd]
            maperravg[igd] /= initmaps.vornorm[igd]
         endif

         if hasrangefile then begin
            ithisline = where(rangeline eq ranlin AND $
                              rangequant eq ranlab,ctthisline)
            if ctthisline eq 1 then auto=0b
         endif else auto=1b
         plotdat = $
            ifsf_plotrange(auto=auto,$
                           mapgd=map[igd],divinit=cbdivinit,ncbdivmax=ncbdivmax,$
                           rline=rangeline,matline=ranlin,$
                           rquant=rangequant,matquant=ranlab,$
                           rncbdiv=rangencbdiv,rlo=rangelo,rhi=rangehi)

         map[ibd] = bad

;        Save some things for later
         igd_nadem_fitflux = igd
         ibd_nadem_fitflux = ibd
         plotdat_nadem_fitflux = plotdat
         map_nadem_fitflux = map
         map_nadem_fitflux_errlo = maperrlo
         map_nadem_fitflux_errhi = maperrhi

         mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
                         min=plotdat[0],max=plotdat[1])
         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos_mid[*,1],opos=truepos,$
                 missing_value=bad,missing_index=255,$
                 missing_color='white',/noerase
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,$
                          center_nuclei_kpc_x,center_nuclei_kpc_y,/nolab
         xoffset = pan_xfrac*0.05
         yoffset = mar_yfrac*0.2
         cbpos=[truepos[0]+xoffset,truepos[1]-yoffset,$
                truepos[2]-xoffset,truepos[1]]
         ticknames = string(dindgen(plotdat[3]+1)*$
                            plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,charsize=0.6
         yoffset = mar_yfrac*0.85
         cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
                truepos[1]-yoffset,panel_title,$
                charsize=panel_title_chars,align=0.5,/norm

;     S/N (fitted)
;-------------------------------------------------------------------------------

      ranlab = 'fitfluxsnr'
      panel_title = textoidl('I^{em}/\deltaI^{em}')
      cbdivinit=5d

      map = map_nadem_fitflux
      maperr = sqrt((map_nadem_fitflux_errlo^2d + map_nadem_fitflux_errhi^2d)/2d)
      igd = igd_nadem_fitflux
      map[igd] /= maperr[igd]

      if hasrangefile then begin
         ithisline = where(rangeline eq ranlin AND $
                           rangequant eq ranlab,ctthisline)
         if ctthisline eq 1 then auto=0b
      endif else auto=1b
      plotdat = $
         ifsf_plotrange(auto=auto,$
                        mapgd=map[igd],divinit=cbdivinit,ncbdivmax=ncbdivmax,$
                        rline=rangeline,matline=ranlin,$
                        rquant=rangequant,matquant=ranlab,$
                        rncbdiv=rangencbdiv,rlo=rangelo,rhi=rangehi)


      mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
                      min=plotdat[0],max=plotdat[1])
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos_mid[*,2],opos=truepos,$
              missing_value=bad,missing_index=255,$
              missing_color='white',/noerase
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
             /nodata,/noerase
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,$
                       center_nuclei_kpc_x,center_nuclei_kpc_y,/nolab
      xoffset = pan_xfrac*0.05
      yoffset = mar_yfrac*0.2
      cbpos=[truepos[0]+xoffset,truepos[1]-yoffset,$
             truepos[2]-xoffset,truepos[1]]
      ticknames = string(dindgen(plotdat[3]+1)*$
                         plotdat[2]/double(plotdat[3]) - $
                         (plotdat[2] - plotdat[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=plotdat[3],$
                 ticknames=ticknames,charsize=0.6
      yoffset = mar_yfrac*0.85
      cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
             truepos[1]-yoffset,panel_title,$
             charsize=panel_title_chars,align=0.5,/norm

;     Fitted W_eq vs. Empirical W_eq
;-------------------------------------------------------------------------------
      igderr = where(map_nadem_fitflux ne 0d AND map_nadem_fitflux ne bad AND $
                     map_nadem_empflux ne 0d AND map_nadem_empflux ne bad,ctgderr)

      if ctgderr gt 0 then begin
         ibderr = where(map_nadem_fitflux eq 0d OR map_nadem_fitflux eq bad OR $
                        map_nadem_empflux eq 0d OR map_nadem_empflux eq bad)
         map_nadem_fitflux_errlo[ibderr] = 0d
         map_nadem_fitflux_errhi[ibderr] = 0d
         map_nadem_empflux_err[ibderr] = 0d

         xran = [0d,max([map_nadem_fitflux[igd_nadem_fitflux],$
                      map_nadem_empflux[igd_nadem_empflux]])*1.1d]
         yran = xran
         cgplot,map_nadem_empflux,map_nadem_fitflux,/xsty,/ysty,$
                xran=xran,yran=yran,psym=3,pos=pos_mid[*,3],$
                xtit=textoidl('I^{em} (emp)'),$
                ytit=textoidl('I^{em} (fit)'),$
                /noerase,aspect=1d,$
                chars=0.8,err_width=0,err_color='Gray',/err_clip,$
                err_xlow=map_nadem_empflux_err,err_xhigh=map_nadem_empflux_err,$
                err_ylow=map_nadem_fitflux_errlo,err_yhigh=map_nadem_fitflux_errhi
         cgoplot,map_nadem_empflux,map_nadem_fitflux,psym=16,symsize=0.4
         cgoplot,xran,xran
         inz = where(map_nadem_fitflux ne 0d AND map_nadem_fitflux ne bad AND $
                     map_nadem_empflux ne 0d AND map_nadem_empflux ne bad)
         flux_rms = sqrt(mean((map_nadem_fitflux[inz]-map_nadem_empflux[inz])^2d))
         cgoplot,xran-flux_rms,xran,linesty=3
         cgoplot,xran+flux_rms,xran,linesty=3
      endif

;     Velocities
;-------------------------------------------------------------------------------
        
;        sigma
;-------------------------------------------------------------------------------
         panel_title = '$\sigma$ (km/s)'
         ranlab = 'fitvsig'
         cbdivinit=100d
         
         map = nademcvdfvals['vsig']
         igd_thiscomp = where(map ne bad,ctgd_thiscomp)
         ibd_thiscomp = where(map eq bad,ctbd_thiscomp)
         if ctgd_thiscomp gt 0 then $
            igd = cgsetintersection(igd_nadem_fitweq,igd_thiscomp) $
         else igd=igd_nadem_fitweq
         if ctbd_thiscomp gt 0 then $
            ibd = cgsetunion(ibd_nadem_fitweq,ibd_thiscomp) $
         else ibd = ibd_nadem_fitweq

         map[ibd] = bad

         if hasrangefile then begin
            ithisline = where(rangeline eq ranlin AND $
                              rangequant eq ranlab,ctthisline)
            if ctthisline eq 1 then auto=0b
         endif else auto=1b
         plotdat = $
            ifsf_plotrange(auto=auto,$
                           mapgd=map[igd],divinit=cbdivinit,ncbdivmax=ncbdivmax,$
                           rline=rangeline,matline=ranlin,$
                           rquant=rangequant,matquant=ranlab,$
                           rncbdiv=rangencbdiv,rlo=rangelo,rhi=rangehi)

         mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
                         min=plotdat[0],max=plotdat[1])
         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos_bot[*,0],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,$
                          center_nuclei_kpc_y,/noxlab
         xoffset = pan_xfrac*0.05
         yoffset = mar_yfrac*0.2
         cbpos=[truepos[0]+xoffset,truepos[1]-yoffset,$
                truepos[2]-xoffset,truepos[1]]
         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         ticknames[0]=' '
         ticknames[plotdat[3]]=' '
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,charsize=0.6
         yoffset = mar_yfrac*0.85
         cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
                truepos[1]-yoffset,panel_title,$
                charsize=panel_title_chars,align=0.5,/norm
         
;        vpk
;-------------------------------------------------------------------------------
         panel_title = 'v$\downpk$ (km/s)'
;         ranlab = 'fitc'+string(i+1,format='(I0)')+'vpk'
         ranlab = 'fitvpk'
         cbdivinit=100d
         
;         map = nademvel[*,*,i]
         map = nademcvdfvals['vpk']
         igd_thiscomp = where(map ne bad,ctgd_thiscomp)
         ibd_thiscomp = where(map eq bad,ctbd_thiscomp)
         if ctgd_thiscomp gt 0 then $
            igd = cgsetintersection(igd_nadem_fitweq,igd_thiscomp) $
         else igd=igd_nadem_fitweq
         if ctbd_thiscomp gt 0 then $
            ibd = cgsetunion(ibd_nadem_fitweq,ibd_thiscomp) $
         else ibd = ibd_nadem_fitweq

         map[ibd] = bad

         if hasrangefile then begin
            ithisline = where(rangeline eq ranlin AND $
                              rangequant eq ranlab,ctthisline)
            if ctthisline eq 1 then auto=0b
         endif else auto=1b
         plotdat = $
            ifsf_plotrange(auto=auto,$
                           mapgd=map[igd],divinit=cbdivinit,ncbdivmax=ncbdivmax,$
                           rline=rangeline,matline=ranlin,$
                           rquant=rangequant,matquant=ranlab,$
                           rncbdiv=rangencbdiv,rlo=rangelo,rhi=rangehi)

         mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
                         min=plotdat[0],max=plotdat[1])
         cgloadct,74,/reverse
         cgimage,mapscl,/keep,pos=pos_bot[*,1],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase,xran=[0,dx],yran=[0,dy]
         maptmp = map*0d
         iblue = where(map lt 0 AND map ne bad,ctblue)
         if ctblue gt 0 then maptmp[iblue] = 1d
         cgcontour,maptmp,dindgen(dx)+0.5,dindgen(dy)+0.5,$
                   /overplot,color=0,c_thick=4,$
                   levels=[1d]
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,$
                          center_nuclei_kpc_y,/nolab
         xoffset = pan_xfrac*0.05
         yoffset = mar_yfrac*0.2
         cbpos=[truepos[0]+xoffset,truepos[1]-yoffset,$
                truepos[2]-xoffset,truepos[1]]
         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         ticknames[0]=' '
         ticknames[plotdat[3]]=' '
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,charsize=0.6
         yoffset = mar_yfrac*0.85
         cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
                truepos[1]-yoffset,panel_title,$
                charsize=panel_title_chars,align=0.5,/norm

         
;        v50
;-------------------------------------------------------------------------------
         panel_title = 'v$\down50$ (km/s)'
;         ranlab = 'fitc'+string(i+1,format='(I0)')+'v%50'
         ranlab = 'fitv%50'
         cbdivinit=100d
         
;         map = nademvel[*,*,i]
         map = nademcvdfvals['v%50']
         igd_thiscomp = where(map ne bad,ctgd_thiscomp)
         ibd_thiscomp = where(map eq bad,ctbd_thiscomp)
         if ctgd_thiscomp gt 0 then $
            igd = cgsetintersection(igd_nadem_fitweq,igd_thiscomp) $
         else igd=igd_nadem_fitweq
         if ctbd_thiscomp gt 0 then $
            ibd = cgsetunion(ibd_nadem_fitweq,ibd_thiscomp) $
         else ibd = ibd_nadem_fitweq

         map[ibd] = bad

         if hasrangefile then begin
            ithisline = where(rangeline eq ranlin AND $
                              rangequant eq ranlab,ctthisline)
            if ctthisline eq 1 then auto=0b
         endif else auto=1b
         plotdat = $
            ifsf_plotrange(auto=auto,$
                           mapgd=map[igd],divinit=cbdivinit,ncbdivmax=ncbdivmax,$
                           rline=rangeline,matline=ranlin,$
                           rquant=rangequant,matquant=ranlab,$
                           rncbdiv=rangencbdiv,rlo=rangelo,rhi=rangehi)

         mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
                         min=plotdat[0],max=plotdat[1])
         cgloadct,74,/reverse
         cgimage,mapscl,/keep,pos=pos_bot[*,2],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase,xran=[0,dx],yran=[0,dy]
         maptmp = map*0d
         iblue = where(map lt 0 AND map ne bad,ctblue)
         if ctblue gt 0 then maptmp[iblue] = 1d
         cgcontour,maptmp,dindgen(dx)+0.5,dindgen(dy)+0.5,$
                   /overplot,color=0,c_thick=4,$
                   levels=[1d]
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,$
                          center_nuclei_kpc_y,/nolab
         xoffset = pan_xfrac*0.05
         yoffset = mar_yfrac*0.2
         cbpos=[truepos[0]+xoffset,truepos[1]-yoffset,$
                truepos[2]-xoffset,truepos[1]]
         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         ticknames[0]=' '
         ticknames[plotdat[3]]=' '
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,charsize=0.6
         yoffset = mar_yfrac*0.85
         cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
                truepos[1]-yoffset,panel_title,$
                charsize=panel_title_chars,align=0.5,/norm
         
;        v98
;-------------------------------------------------------------------------------
         panel_title = 'v$\down98$ (km/s)'
;         ranlab = 'fitc'+string(i+1,format='(I0)')+'v%98'
         ranlab = 'fitv%98'
         cbdivinit=200d
         
;         map = nademv98[*,*,i]
         map = nademcvdfvals['v%98']
         igd_thiscomp = where(map ne bad,ctgd_thiscomp)
         ibd_thiscomp = where(map eq bad,ctbd_thiscomp)
         if ctgd_thiscomp gt 0 then $
            igd = cgsetintersection(igd_nadem_fitweq,igd_thiscomp) $
         else igd=igd_nadem_fitweq
         if ctbd_thiscomp gt 0 then $
            ibd = cgsetunion(ibd_nadem_fitweq,ibd_thiscomp) $
         else ibd = ibd_nadem_fitweq

         map[ibd] = bad

         if hasrangefile then begin
            ithisline = where(rangeline eq ranlin AND $
                              rangequant eq ranlab,ctthisline)
            if ctthisline eq 1 then auto=0b
         endif else auto=1b
         plotdat = $
            ifsf_plotrange(auto=auto,$
                           mapgd=map[igd],divinit=cbdivinit,ncbdivmax=ncbdivmax,$
                           rline=rangeline,matline=ranlin,$
                           rquant=rangequant,matquant=ranlab,$
                           rncbdiv=rangencbdiv,rlo=rangelo,rhi=rangehi)

         mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
                         min=plotdat[0],max=plotdat[1])
         cgloadct,74,/reverse
         cgimage,mapscl,/keep,pos=pos_bot[*,3],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,$
                          center_nuclei_kpc_y,/nolab
         xoffset = pan_xfrac*0.05
         yoffset = mar_yfrac*0.2
         cbpos=[truepos[0]+xoffset,truepos[1]-yoffset,$
                truepos[2]-xoffset,truepos[1]]
         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         ticknames[0]=' '
         ticknames[plotdat[3]]=' '
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,charsize=0.6
         yoffset = mar_yfrac*0.85
         cgtext,truepos[0]+(truepos[2]-truepos[0])/2d,$
                truepos[1]-yoffset,panel_title,$
                charsize=panel_title_chars,align=0.5,/norm
;
;      endfor
       endif




;     Title
;-------------------------------------------------------------------------------
      cgplot,[0],xsty=5,ysty=5,position=[0,0,1,1],/nodata,/noerase
      cgtext,pos_title[0],pos_title[1],$
             initdat.name+': NaI D emission',$
             charsize=1.25d,align=0.5

      cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; W_eq (empirical)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;      if maxnademncomp_act gt 0 then ny=3 else ny=1
;
;      cgps_open,initdat.mapdir+initdat.label+'NaDempweq.eps',charsize=1,/encap,$
;         /inches,xs=plotquantum*2,ys=plotquantum*ny*aspectrat,/qui
;
;      pos = cglayout([2,ny],ixmar=[2,2],iymar=[2,2],oxmar=[0,0],oymar=[0,0],$
;         xgap=0,ygap=0,unit=!D.X_PX_CM/3.0)
;      cbform = '(I0)'
;
;;
;;     ABSORPTION
;;
;      map = nadcube.weq[*,*,0]
;      igd = where(map ge initmaps.nadabsweq_snrthresh*nadcube.weq[*,*,1] AND $
;                  map gt 0d AND map ne bad)
;      ibd = where(map lt initmaps.nadabsweq_snrthresh*nadcube.weq[*,*,1] OR $
;                  map eq 0d OR map eq bad)
;
;;     Set up range
;;     Check for manual range first ...
;      hasrange = 0
;      if hasrangefile then begin
;         ithisline = where(rangeline eq 'NaDabs' AND $
;                           rangequant eq 'empweq',ctthisline)
;         if ctthisline eq 1 then begin
;            zran = [rangelo[ithisline],rangehi[ithisline]]
;            dzran = zran[1]-zran[0]
;            ncbdiv = rangencbdiv[ithisline]
;            ncbdiv = ncbdiv[0]
;            hasrange = 1
;         endif
;      endif
;;     otherwise set it automagically.
;      if ~hasrange then begin
;         zran = [min(map[igd]),max(map[igd])]
;         divarr = ifsf_cbdiv(zran,2d,ncbdivmax)
;         ncbdiv = divarr[0]
;         dzran = zran[1]-zran[0]
;      endif
;
;;     replace bad points
;      map[ibd] = bad
;
;;     Save some things for later
;      zran_empweqabs = zran
;      map_empweqabs = map
;
;      mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;                      min=zran[0],max=zran[1])
;      cgloadct,65,/reverse
;      cgimage,mapscl,/keep,pos=pos[*,0],opos=truepos,$
;         missing_value=bad,missing_index=255,$
;         missing_color='white'
;      cgplot,[0],xsty=5,ysty=5,position=truepos,$
;;         /nodata,/noerase,title='W$\down eq$(NaD abs, $\angstrom$)'
;         /nodata,/noerase,title=textoidl('W_{eq}^{abs}')+' ($\angstrom$)'
;      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
;         (dzran - zran[1]),format=cbform)
;      cgcolorbar,position=cbpos,divisions=ncbdiv,$
;         ticknames=ticknames,/ver,/right,charsize=0.6
;
;;     Error
;      err = nadcube.weq[*,*,1]
;      map[igd] /= err[igd]
;
;;     Set up range
;;     Check for manual range first ...
;      hasrange = 0
;      if hasrangefile then begin
;         ithisline = where(rangeline eq 'NaDabs' AND $
;            rangequant eq 'weqsnr',ctthisline)
;         if ctthisline eq 1 then begin
;            zran = [rangelo[ithisline],rangehi[ithisline]]
;            dzran = zran[1]-zran[0]
;            ncbdiv = rangencbdiv[ithisline]
;            ncbdiv = ncbdiv[0]
;            hasrange = 1
;         endif
;      endif
;;     otherwise set it automagically.
;      if ~hasrange then begin
;         zran = [min(map[igd]),max(map[igd])]
;         divarr = ifsf_cbdiv(zran,10d,ncbdivmax)
;         ncbdiv = divarr[0]
;         dzran = zran[1]-zran[0]
;      endif
;
;      mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;         min=zran[0],max=zran[1])
;      cgloadct,65,/reverse
;      cgimage,mapscl,/keep,pos=pos[*,1],opos=truepos,$
;         /noerase,missing_value=bad,missing_index=255,$
;         missing_color='white'
;      cgplot,[0],xsty=5,ysty=5,position=truepos,$
;         /nodata,/noerase,title=textoidl('W_{eq}^{abs}/\deltaW')
;      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
;         (dzran - zran[1]),format=cbform)
;      cgcolorbar,position=cbpos,divisions=ncbdiv,$
;         ticknames=ticknames,/ver,/right,charsize=0.6
;
;;
;;     EMISSION
;;
;
;      if ny gt 1 then begin
;
;      map = nadcube.weq[*,*,2]
;      igd = where(abs(map) ge initmaps.nademweq_snrthresh*nadcube.weq[*,*,3] AND $
;                  abs(map) gt 0d AND map ne bad)
;      ibd = where(abs(map) lt initmaps.nademweq_snrthresh*nadcube.weq[*,*,3] OR $
;                  map eq 0d OR map eq bad)
;
;;     Set up range
;;     Check for manual range first ...
;      hasrange = 0
;      if hasrangefile then begin
;         ithisline = where(rangeline eq 'NaDem' AND $
;            rangequant eq 'weq',ctthisline)
;         if ctthisline eq 1 then begin
;            zran = [rangelo[ithisline],rangehi[ithisline]]
;            dzran = zran[1]-zran[0]
;            ncbdiv = rangencbdiv[ithisline]
;            ncbdiv = ncbdiv[0]
;            hasrange = 1
;         endif
;      endif
;;     otherwise set it automagically.
;      if ~hasrange then begin
;         zran = [min(map[igd]),max(map[igd])]
;         divarr = ifsf_cbdiv(zran,2d,ncbdivmax)
;         ncbdiv = divarr[0]
;         dzran = zran[1]-zran[0]
;      endif
;
;;     replace bad points
;      map[ibd] = bad
;
;;     Save some things for later
;      zran_empweqem = zran
;      map_empweqem = map
;
;      mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;         min=zran[0],max=zran[1])
;      cgloadct,65
;      cgimage,mapscl,/keep,pos=pos[*,2],opos=truepos,$
;         /noerase,missing_value=bad,missing_index=255,$
;         missing_color='white'
;      cgplot,[0],xsty=5,ysty=5,position=truepos,$
;         /nodata,/noerase,title=textoidl('W_{eq}^{em}')+' ($\angstrom$)'
;      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
;         (dzran - zran[1]),format=cbform)
;      cgcolorbar,position=cbpos,divisions=ncbdiv,$
;         ticknames=ticknames,/ver,/right,charsize=0.6
;
;;     Error
;      err = nadcube.weq[*,*,3]
;      map[igd] /= -err[igd]
;
;;     Set up range
;;     Check for manual range first ...
;      hasrange = 0
;      if hasrangefile then begin
;         ithisline = where(rangeline eq 'NaDem' AND $
;            rangequant eq 'weqsnr',ctthisline)
;         if ctthisline eq 1 then begin
;            zran = [rangelo[ithisline],rangehi[ithisline]]
;            dzran = zran[1]-zran[0]
;            ncbdiv = rangencbdiv[ithisline]
;            ncbdiv = ncbdiv[0]
;            hasrange = 1
;         endif
;      endif
;;     otherwise set it automagically.
;      if ~hasrange then begin
;         zran = [min(map[igd]),max(map[igd])]
;         divarr = ifsf_cbdiv(zran,2d,ncbdivmax)
;         ncbdiv = divarr[0]
;         dzran = zran[1]-zran[0]
;      endif
;
;      mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;         min=zran[0],max=zran[1])
;      cgloadct,65,/reverse
;      cgimage,mapscl,/keep,pos=pos[*,3],opos=truepos,$
;         /noerase,missing_value=bad,missing_index=255,$
;         missing_color='white'
;      cgplot,[0],xsty=5,ysty=5,position=truepos,$
;         /nodata,/noerase,title=textoidl('W_{eq}^{em}/\deltaW')
;      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
;         (dzran - zran[1]),format=cbform)
;      cgcolorbar,position=cbpos,divisions=ncbdiv,$
;         ticknames=ticknames,/ver,/right,charsize=0.6
;
;
;
;;     Flux
;      cbform = '(D0.2)'
;      map = nadcube.emflux[*,*,0]
;      igd = where(abs(nadcube.weq[*,*,2]) ge $
;                  initmaps.nademweq_snrthresh*nadcube.weq[*,*,3] AND $
;                  abs(map) gt 0d AND map ne bad)
;      ibd = where(abs(nadcube.weq[*,*,2]) lt $
;                  initmaps.nademweq_snrthresh*nadcube.weq[*,*,3] OR $
;                  map eq 0d OR map eq bad)
;
;      if tag_exist(initmaps,'fluxfactor') then $
;         map[igd] *= initmaps.fluxfactor
;         
;;     Set up range
;;     Check for manual range first ...
;      hasrange = 0
;      if hasrangefile then begin
;         ithisline = where(rangeline eq 'NaDem' AND $
;            rangequant eq 'flux',ctthisline)
;         if ctthisline eq 1 then begin
;            zran = [rangelo[ithisline],rangehi[ithisline]]
;            dzran = zran[1]-zran[0]
;            ncbdiv = rangencbdiv[ithisline]
;            ncbdiv = ncbdiv[0]
;            hasrange = 1
;         endif
;      endif
;;     otherwise set it automagically.
;      if ~hasrange then begin
;         zran = [min(map[igd]),max(map[igd])]
;         divarr = ifsf_cbdiv(zran,2d,ncbdivmax)
;         ncbdiv = divarr[0]
;         dzran = zran[1]-zran[0]
;      endif
;      
;;     replace bad points
;      map[ibd] = bad
;
;;     Save some things for later
;      zran_empfluxem = zran
;      map_empfluxem = map
;
;      mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;         min=zran[0],max=zran[1])
;      cgloadct,65,/reverse
;      cgimage,mapscl,/keep,pos=pos[*,4],opos=truepos,$
;         /noerase,missing_value=bad,missing_index=255,$
;         missing_color='white'
;      cgplot,[0],xsty=5,ysty=5,position=truepos,$
;         /nodata,/noerase,$
;         title=textoidl('I^{em} / 10^{-16}')
;      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
;         (dzran - zran[1]),format=cbform)
;      cgcolorbar,position=cbpos,divisions=ncbdiv,$
;         ticknames=ticknames,/ver,/right,charsize=0.6
;
;      endif
;
;      cgps_close
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; W_eq (fits)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;      if maxnademncomp_act gt 0 then ny=3 else ny=1
;
;      cgps_open,initdat.mapdir+initdat.label+'NaDfitweq.eps',charsize=1,/encap,$
;         /inches,xs=plotquantum*2,ys=plotquantum*ny*aspectrat,/qui
;
;      pos = cglayout([2,ny],ixmar=[3,3],iymar=[3,3],oxmar=[0,0],oymar=[0,0],$
;         xgap=0,ygap=0,unit=!D.X_PX_CM/3.0)
;      cbform = '(I0)'
;
;;
;;     ABSORPTION
;;
;      map = nadfit.weqabs[*,*,0]
;      maperrlo = nadfit.weqabserr[*,*,0]
;      maperrhi = nadfit.weqabserr[*,*,1]
;      maperravg = (maperrlo + maperrhi) / 2d
;      maperrlo_emp = nadcube.weq[*,*,1]
;      maperrhi_emp = nadcube.weq[*,*,1]
;      igd = where(map gt 0d AND $
;                  map ne bad AND $
;                  map ge initmaps.nadabsweq_snrthresh*maperravg)
;      ibd = where(map eq 0d OR $
;                  map eq bad OR $
;                  map lt initmaps.nadabsweq_snrthresh*maperravg)
;      igd_nadfitabsweq = igd
;      ibd_nadfitabsweq = ibd
;
;;     Set up range
;;     Check for manual range first ...
;      hasrange = 0
;      if hasrangefile then begin
;         ithisline = where(rangeline eq 'NaDabs' AND $
;                           rangequant eq 'fitweq',ctthisline)
;         if ctthisline eq 1 then begin
;            zran = [rangelo[ithisline],rangehi[ithisline]]
;            dzran = zran[1]-zran[0]
;            ncbdiv = rangencbdiv[ithisline]
;            ncbdiv = ncbdiv[0]
;            hasrange = 1
;         endif
;      endif
;;     otherwise set it automagically.
;      if ~hasrange then begin
;         zran = [min(map[igd]),max(map[igd])]
;         divarr = ifsf_cbdiv(zran,2d,ncbdivmax)
;         ncbdiv = divarr[0]
;         dzran = zran[1]-zran[0]
;      endif
;      zran_fitweqabs = zran
;
;;     replace bad points
;      map[ibd] = bad
;
;      mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;                      min=zran[0],max=zran[1])
;      cgloadct,65,/reverse
;      cgimage,mapscl,/keep,pos=pos[*,0],opos=truepos,$
;         missing_value=bad,missing_index=255,missing_color='white'
;      cgplot,[0],xsty=5,ysty=5,position=truepos,$
;         /nodata,/noerase,title=textoidl('W_{eq}^{abs}')+' ($\angstrom$)'
;      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;      ifsf_plotcompass,xarr_kpc,yarr_kpc
;      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
;         (dzran - zran[1]),format=cbform)
;      cgcolorbar,position=cbpos,divisions=ncbdiv,$
;         ticknames=ticknames,/ver,/right,charsize=0.6
;
;      ibderr = where(map eq 0d OR map eq bad OR $
;                     map_empweqabs eq 0d OR map_empweqabs eq bad)
;      maperrlo[ibderr] = 0d
;      maperrhi[ibderr] = 0d
;      maperrlo_emp[ibderr] = 0d
;      maperrhi_emp[ibderr] = 0d
;
;      igd_empweqabs = where(map_empweqabs ne bad)
;      xran = [0d,max([map[igd],map_empweqabs[igd_empweqabs]])*1.1d]
;      yran = xran
;      cgplot,map_empweqabs,map,/xsty,/ysty,$
;         xran=xran,yran=yran,psym=3,pos=pos[*,1],$
;         xtit=textoidl('W_{eq}^{abs} (empirical)'),$
;         ytit=textoidl('W_{eq}^{abs} (fit)'),$
;         /noerase,aspect=1d,$
;         chars=0.8,err_width=0,err_color='Gray',/err_clip,$
;         err_xlow=maperrlo_emp,err_xhigh=maperrhi_emp,$
;         err_ylow=maperrlo,err_yhigh=maperrhi
;      cgoplot,map_empweqabs,map,psym=16,symsize=0.4
;      cgoplot,xran,xran
;      inz = where(map ne 0d AND map ne bad AND $
;                  map_empweqabs ne 0d AND map_empweqabs ne bad)
;      weq_rms = sqrt(mean((map[inz]-map_empweqabs[inz])^2d))
;      cgoplot,xran-weq_rms,xran,linesty=3
;      cgoplot,xran+weq_rms,xran,linesty=3
;;
;;     EMISSION
;;
;
;      if ny gt 1 then begin
;
;      map = nadfit.weqem[*,*,0]
;      maperrlo = nadfit.weqemerr[*,*,0]
;      maperrhi = nadfit.weqemerr[*,*,1]
;      maperravg = (maperrlo + maperrhi)/2d
;      maperrlo_emp = nadcube.weq[*,*,3]
;      maperrhi_emp = nadcube.weq[*,*,3]
;;      igd = where(abs(map) gt 0d AND map ne bad)
;;      ibd = where(map eq 0d OR map eq bad)
;      igd = where(abs(map) gt 0d AND $
;                  map ne bad)
;      ibd = where(map eq 0d OR $
;                  map eq bad)
;      igd_nadfitemweq = igd
;      ibd_nadfitemweq = ibd
;
;;     Set up range
;;     Check for manual range first ...
;      hasrange = 0
;      if hasrangefile then begin
;         ithisline = where(rangeline eq 'NaDem' AND $
;            rangequant eq 'weq',ctthisline)
;         if ctthisline eq 1 then begin
;            zran = [rangelo[ithisline],rangehi[ithisline]]
;            dzran = zran[1]-zran[0]
;            ncbdiv = rangencbdiv[ithisline]
;            ncbdiv = ncbdiv[0]
;            hasrange = 1
;         endif
;      endif
;;     otherwise set it automagically.
;      if ~hasrange then begin
;         zran = [min(map[igd]),max(map[igd])]
;         divarr = ifsf_cbdiv(zran,2d,ncbdivmax)
;         ncbdiv = divarr[0]
;         dzran = zran[1]-zran[0]
;      endif
;
;;     replace bad points
;      map[ibd] = 0d
;
;      mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;         min=zran[0],max=zran[1])
;      cgloadct,65
;      cgimage,mapscl,/keep,pos=pos[*,2],opos=truepos,$
;         /noerase,missing_value=bad,missing_index=255,$
;         missing_color='white'
;      cgplot,[0],xsty=5,ysty=5,position=truepos,$
;         /nodata,/noerase,title=textoidl('W_{eq}^{em}')+' ($\angstrom$)'
;      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
;         (dzran - zran[1]),format=cbform)
;      cgcolorbar,position=cbpos,divisions=ncbdiv,$
;         ticknames=ticknames,/ver,/right,charsize=0.6
;
;      ibderr = where(map eq 0d OR map eq bad OR $
;                     map_empweqem eq 0d OR map_empweqem eq bad)
;      maperrlo[ibderr] = 0d
;      maperrhi[ibderr] = 0d
;      maperrlo_emp[ibderr] = 0d
;      maperrhi_emp[ibderr] = 0d
;
;      igd_empweqem = where(map_empweqem ne bad)
;      xran = [min([map[igd],map_empweqem[igd_empweqem]])*1.1d,0d]
;      yran = xran
;      cgplot,map_empweqem,map,/xsty,/ysty,xran=xran,yran=yran,psym=3,$
;         pos=pos[*,3],/noerase,aspect=1d,$
;         xtit=textoidl('W_{eq}^{em} (empirical)'),$
;         ytit=textoidl('W_{eq}^{em} (fit)'),$    
;         chars=0.8d,err_color='Gray',err_width=0,/err_clip,$
;         err_xlow=maperrlo_emp,err_xhigh=maperrhi_emp,$
;         err_ylow=maperrlo,err_yhigh=maperrhi
;      cgoplot,map_empweqem,map,psym=16,symsize=0.4d
;      cgoplot,xran,xran
;      inz = where(map ne 0d AND map ne bad AND $
;                  map_empweqem ne 0d AND map_empweqem ne bad)
;      weq_rms = sqrt(mean((map[inz]-map_empweqem[inz])^2d))
;      cgoplot,xran-weq_rms,xran,linesty=3
;      cgoplot,xran+weq_rms,xran,linesty=3
;
;
;;     Flux
;      cbform = '(D0.2)'
;      map = nadfit.totfluxem[*,*,0]
;      maperrlo = nadfit.totfluxemerr[*,*,0]
;      maperrhi = nadfit.totfluxemerr[*,*,1]
;      maperrlo_emp = nadcube.emflux[*,*,1]
;      maperrhi_emp = nadcube.emflux[*,*,1]
;;      igd = where(abs(map) gt 0d AND map ne bad)
;;      ibd = where(map eq 0d OR map eq bad)
;      igd = igd_nadfitemweq
;      ibd = ibd_nadfitemweq
;
;      if tag_exist(initmaps,'fluxfactor') then $
;         map[igd] *= initmaps.fluxfactor
;
;;     Set up range
;;     Check for manual range first ...
;      hasrange = 0
;      if hasrangefile then begin
;         ithisline = where(rangeline eq 'NaDem' AND $
;            rangequant eq 'flux',ctthisline)
;         if ctthisline eq 1 then begin
;            zran = [rangelo[ithisline],rangehi[ithisline]]
;            dzran = zran[1]-zran[0]
;            ncbdiv = rangencbdiv[ithisline]
;            ncbdiv = ncbdiv[0]
;            hasrange = 1
;         endif
;      endif
;;     otherwise set it automagically.
;      if ~hasrange then begin
;         zran = [0,max(map[igd])]
;         divarr = ifsf_cbdiv(zran,initmaps.nademflux_cbint,ncbdivmax)
;         ncbdiv = divarr[0]
;         dzran = zran[1]-zran[0]
;      endif
;      
;;     replace bad points
;      map[ibd] = bad
;
;      mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;         min=zran[0],max=zran[1])
;      cgloadct,65,/reverse
;      cgimage,mapscl,/keep,pos=pos[*,4],opos=truepos,$
;         /noerase,missing_value=bad,missing_index=255,$
;         missing_color='white'
;      cgplot,[0],xsty=5,ysty=5,position=truepos,$
;         /nodata,/noerase,title='I$\upem$ / 10$\up-16$'
;      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
;         (dzran - zran[1]),format=cbform)
;      cgcolorbar,position=cbpos,divisions=ncbdiv,$
;         ticknames=ticknames,/ver,/right,charsize=0.6
;
;      igderr = where(map ne 0d AND map ne bad AND $
;                     map_empfluxem ne 0d AND map_empfluxem ne bad)
;      ibderr = where(map eq 0d OR map eq bad OR $
;                     map_empfluxem eq 0d OR map_empfluxem eq bad)
;      maperrlo[ibderr] = 0d
;      maperrhi[ibderr] = 0d
;      maperrlo_emp[ibderr] = 0d
;      maperrhi_emp[ibderr] = 0d
;      if tag_exist(initmaps,'fluxfactor') then begin
;         maperrlo[igderr] *= initmaps.fluxfactor
;         maperrhi[igderr] *= initmaps.fluxfactor
;         maperrlo_emp[igderr] *= initmaps.fluxfactor
;         maperrhi_emp[igderr] *= initmaps.fluxfactor
;      endif
;
;      igd_empfluxem = where(map_empfluxem ne bad)
;      xyran = [0d,max([map[igd],map_empfluxem[igd_empfluxem]])*1.1d]
;      cgplot,map_empfluxem,map,/xsty,/ysty,xran=xyran,yran=xyran,psym=3,$
;         pos=pos[*,5],xtit='I$\upem$ / 10$\up-16$ (empirical)',$
;         ytit='I$\upem$ / 10$\up-16$ (fit)',/noerase,aspect=1d,$
;         chars=0.8d,err_color='Gray',err_width=0,/err_clip,$
;         err_xlow=maperrlo_emp,err_xhigh=maperrhi_emp,$
;         err_ylow=maperrlo,err_yhigh=maperrhi
;      cgoplot,map_empfluxem,map,psym=16,symsize=0.4d
;      cgoplot,xyran,xyran
;      inz = where(map ne 0d AND map ne bad AND $
;                  map_empfluxem ne 0d AND map_empfluxem ne bad)
;      weq_rms = sqrt(mean((map[inz]-map_empfluxem[inz])^2d))
;      cgoplot,xyran-weq_rms,xyran,linesty=3
;      cgoplot,xyran+weq_rms,xyran,linesty=3
;
;      endif
;
;      cgps_close
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; components (fits)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      cgps_open,initdat.mapdir+initdat.label+'NaDfitncomp.eps',charsize=1.5,$
         /encap,$
         /inches,xs=plotquantum*2,ys=plotquantum*2*aspectrat,/qui,/nomatch
      pos = cglayout([1,1],ixmar=[0,0],iymar=[0,0],oxmar=[3,5],oymar=[6,1],$
                     xgap=0,ygap=0,aspect=1,unit=!D.X_PX_CM/3.0)
      map = nadabsncomp
;      ibd = where(map eq bad)
      ibd = ibd_nadabs_fitweq
      map[ibd] = 0l
      
      zran = [0,2]
      ncbdiv = 2
      dzran = zran[1]-zran[0]
      cbform = '(I0)'

      mapscl = cgimgscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
                        minval=zran[0],maxval=zran[1],ncolors=3)
      cgloadct,8,/brewer,ncolors=3
      cgimage,mapscl,pos=pos[*,0]
      cgplot,[0],xsty=5,ysty=5,pos=pos[*,0],$
         /nodata,/noerase ;,title=textoidl('N_{comp}^{abs}')
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
      cbpos=[pos[2,0],pos[1,0],pos[2,0]+0.02,pos[3,0]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgDCBar,ncolors=3,position=cbpos,labels=ticknames,/ver,/right,$
         charsize=1.5,spacing=0.5
      cgtext,0.96,0.57,textoidl('# of Absorption Components'),orient=270,/normal,$
             align=0.5
         
      cgps_close
      
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; velocity (empirical)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;      if maxnademncomp_act gt 0 then ny=3 else ny=1
;
;      cgps_open,initdat.mapdir+initdat.label+'NaDempvel.eps',charsize=1,/encap,$
;         /inches,xs=plotquantum*3,ys=plotquantum*ny*aspectrat,/qui
;
;      pos = cglayout([3,ny],ixmar=[3,3],iymar=[3,3],oxmar=[0,0],oymar=[0,0],$
;         xgap=0,ygap=0,unit=!D.X_PX_CM/3.0)
;      cbform = '(I0)'
;
;;
;;     ABSORPTION
;;
;      map = nadcube.vel[*,*,0]
;      igd = where(nadcube.weq[*,*,0] ge initmaps.nadabsweq_snrthresh*nadcube.weq[*,*,1] $
;                  AND map gt 0d AND map ne bad)
;      ibd = where(nadcube.weq[*,*,0] lt initmaps.nadabsweq_snrthresh*nadcube.weq[*,*,1] $
;                  OR map eq 0d OR map eq bad)
;
;;     Set up range
;;     Check for manual range first ...
;      hasrange = 0
;      if hasrangefile then begin
;         ithisline = where(rangeline eq 'NaDabs' AND $
;            rangequant eq 'empvwidth',ctthisline)
;         if ctthisline eq 1 then begin
;            zran = [rangelo[ithisline],rangehi[ithisline]]
;            dzran = zran[1]-zran[0]
;            ncbdiv = rangencbdiv[ithisline]
;            ncbdiv = ncbdiv[0]
;            hasrange = 1
;         endif
;      endif
;;     otherwise set it automagically.
;      if ~hasrange then begin
;         zran = [min(map[igd]),max(map[igd])]
;         divarr = ifsf_cbdiv(zran,200d,ncbdivmax)
;         ncbdiv = divarr[0]
;         dzran = zran[1]-zran[0]
;      endif
;
;;     replace bad points with "bad"
;      map[ibd] = bad
;
;      mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;         min=zran[0],max=zran[1])
;      cgloadct,65,/reverse
;      cgimage,mapscl,/keep,pos=pos[*,0],opos=truepos,$
;         noerase=i ne 0,missing_value=bad,missing_index=255,$
;         missing_color='white'
;      cgplot,[0],xsty=5,ysty=5,position=truepos,$
;         /nodata,/noerase,title='$\Delta$v(NaD abs, km/s)'
;      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
;         (dzran - zran[1]),format=cbform)
;      cgcolorbar,position=cbpos,divisions=ncbdiv,$
;         ticknames=ticknames,/ver,/right,charsize=0.6
;
;      map = nadcube.vel[*,*,1]
;      igd = where(nadcube.weq[*,*,0] ge initmaps.nadabsweq_snrthresh*nadcube.weq[*,*,1] $
;                  AND map ne bad)
;      ibd = where(nadcube.weq[*,*,0] lt initmaps.nadabsweq_snrthresh*nadcube.weq[*,*,1] $
;                  OR map eq bad)
;
;;     Set up range
;;     Check for manual range first ...
;      hasrange = 0
;      if hasrangefile then begin
;         ithisline = where(rangeline eq 'NaDabs' AND $
;            rangequant eq 'empvavg',ctthisline)
;         if ctthisline eq 1 then begin
;            zran = [rangelo[ithisline],rangehi[ithisline]]
;            dzran = zran[1]-zran[0]
;            ncbdiv = rangencbdiv[ithisline]
;            ncbdiv = ncbdiv[0]
;            hasrange = 1
;         endif
;      endif
;;     otherwise set it automagically.
;      if ~hasrange then begin
;         zran = [min(map[igd]),max(map[igd])]
;         divarr = ifsf_cbdiv(zran,100d,ncbdivmax)
;         ncbdiv = divarr[0]
;         dzran = zran[1]-zran[0]
;      endif
;
;;     replace bad points with "bad"
;      map[ibd] = bad
;
;      mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;         min=zran[0],max=zran[1])
;      cgloadct,74,/reverse
;      cgimage,mapscl,/keep,pos=pos[*,1],opos=truepos,$
;         noerase=i ne 0,missing_value=bad,missing_index=255,$
;         missing_color='white'
;      cgplot,[0],xsty=5,ysty=5,position=truepos,$
;         /nodata,/noerase,title='<v>(NaD abs, km/s)'
;      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
;         (dzran - zran[1]),format=cbform)
;      cgcolorbar,position=cbpos,divisions=ncbdiv,$
;         ticknames=ticknames,/ver,/right,charsize=0.6
;
;
;      map = nadcube.vel[*,*,2]
;      igd = where(nadcube.weq[*,*,0] ge initmaps.nadabsweq_snrthresh*nadcube.weq[*,*,1] $
;                  AND map ne bad)
;      ibd = where(nadcube.weq[*,*,0] lt initmaps.nadabsweq_snrthresh*nadcube.weq[*,*,1] $
;                  OR map eq bad)
;
;;     Set up range
;;     Check for manual range first ...
;      hasrange = 0
;      if hasrangefile then begin
;         ithisline = where(rangeline eq 'NaDabs' AND $
;            rangequant eq 'empvmax',ctthisline)
;         if ctthisline eq 1 then begin
;            zran = [rangelo[ithisline],rangehi[ithisline]]
;            dzran = zran[1]-zran[0]
;            ncbdiv = rangencbdiv[ithisline]
;            ncbdiv = ncbdiv[0]
;            hasrange = 1
;         endif
;      endif
;;     otherwise set it automagically.
;      if ~hasrange then begin
;         zran = [min(map[igd]),max(map[igd])]
;         divarr = ifsf_cbdiv(zran,100d,ncbdivmax)
;         ncbdiv = divarr[0]
;         dzran = zran[1]-zran[0]
;      endif
;
;;     replace bad points with "bad"
;      map[ibd] = bad
;
;      mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;         min=zran[0],max=zran[1])
;      cgloadct,74,/reverse
;      cgimage,mapscl,/keep,pos=pos[*,2],opos=truepos,$
;         noerase=i ne 0,missing_value=bad,missing_index=255,$
;         missing_color='white'
;      cgplot,[0],xsty=5,ysty=5,position=truepos,$
;         /nodata,/noerase,title='v$\downmax$(NaD abs, km/s)'
;      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
;         (dzran - zran[1]),format=cbform)
;      cgcolorbar,position=cbpos,divisions=ncbdiv,$
;         ticknames=ticknames,/ver,/right,charsize=0.6
;
;;
;;     EMISSION
;;
;
;      if ny gt 1 then begin
;
;      map = nadcube.vel[*,*,3]
;      igd = where(abs(nadcube.weq[*,*,2]) ge $
;                  initmaps.nademweq_snrthresh*nadcube.weq[*,*,3] $
;                  AND map gt 0d AND map ne bad)
;      ibd = where(abs(nadcube.weq[*,*,2]) lt $
;                  initmaps.nademweq_snrthresh*nadcube.weq[*,*,3] $
;                  OR map eq 0d OR map eq bad)
;
;;     Set up range
;;     Check for manual range first ...
;      hasrange = 0
;      if hasrangefile then begin
;         ithisline = where(rangeline eq 'NaDem' AND $
;            rangequant eq 'empvelwid',ctthisline)
;         if ctthisline eq 1 then begin
;            zran = [rangelo[ithisline],rangehi[ithisline]]
;            dzran = zran[1]-zran[0]
;            ncbdiv = rangencbdiv[ithisline]
;            ncbdiv = ncbdiv[0]
;            hasrange = 1
;         endif
;      endif
;;     otherwise set it automagically.
;      if ~hasrange then begin
;         zran = [min(map[igd]),max(map[igd])]
;         divarr = ifsf_cbdiv(zran,200d,ncbdivmax)
;         ncbdiv = divarr[0]
;         dzran = zran[1]-zran[0]
;      endif
;
;;     replace bad points with "bad"
;      map[ibd] = bad
;
;      mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;         min=zran[0],max=zran[1])
;      cgloadct,65,/reverse
;      cgimage,mapscl,/keep,pos=pos[*,3],opos=truepos,$
;         noerase=i ne 0,missing_value=bad,missing_index=255,$
;         missing_color='white'
;      cgplot,[0],xsty=5,ysty=5,position=truepos,$
;         /nodata,/noerase,title='$\Delta$v(NaD em, km/s)'
;      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
;         (dzran - zran[1]),format=cbform)
;      cgcolorbar,position=cbpos,divisions=ncbdiv,$
;         ticknames=ticknames,/ver,/right,charsize=0.6
;
;      map = nadcube.vel[*,*,4]
;      igd = where(abs(nadcube.weq[*,*,2]) ge $
;                  initmaps.nademweq_snrthresh*nadcube.weq[*,*,3] $
;                  AND map ne bad)
;      ibd = where(abs(nadcube.weq[*,*,2]) lt $
;                  initmaps.nademweq_snrthresh*nadcube.weq[*,*,3] $
;                  OR map eq bad)
;
;;     Set up range
;;     Check for manual range first ...
;      hasrange = 0
;      if hasrangefile then begin
;         ithisline = where(rangeline eq 'NaDem' AND $
;            rangequant eq 'empvelavg',ctthisline)
;         if ctthisline eq 1 then begin
;            zran = [rangelo[ithisline],rangehi[ithisline]]
;            dzran = zran[1]-zran[0]
;            ncbdiv = rangencbdiv[ithisline]
;            ncbdiv = ncbdiv[0]
;            hasrange = 1
;         endif
;      endif
;;     otherwise set it automagically.
;      if ~hasrange then begin
;         zran = [min(map[igd]),max(map[igd])]
;         divarr = ifsf_cbdiv(zran,100d,ncbdivmax)
;         ncbdiv = divarr[0]
;         dzran = zran[1]-zran[0]
;      endif
;
;;     replace bad points with "bad"
;      map[ibd] = bad
;
;      mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;         min=zran[0],max=zran[1])
;      cgloadct,74,/reverse
;      cgimage,mapscl,/keep,pos=pos[*,4],opos=truepos,$
;         noerase=i ne 0,missing_value=bad,missing_index=255,$
;         missing_color='white'
;      cgplot,[0],xsty=5,ysty=5,position=truepos,$
;         /nodata,/noerase,title='<v>(NaD em, km/s)'
;      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
;         (dzran - zran[1]),format=cbform)
;      cgcolorbar,position=cbpos,divisions=ncbdiv,$
;         ticknames=ticknames,/ver,/right,charsize=0.6
;
;      map = nadcube.vel[*,*,5]
;      igd = where(abs(nadcube.weq[*,*,2]) ge $
;                  initmaps.nademweq_snrthresh*nadcube.weq[*,*,3] $
;                  AND map ne bad)
;      ibd = where(abs(nadcube.weq[*,*,2]) lt $
;                  initmaps.nademweq_snrthresh*nadcube.weq[*,*,3] $
;                  OR map eq bad)
;
;;     Set up range
;;     Check for manual range first ...
;      hasrange = 0
;      if hasrangefile then begin
;         ithisline = where(rangeline eq 'NaDem' AND $
;            rangequant eq 'empvelmax',ctthisline)
;         if ctthisline eq 1 then begin
;            zran = [rangelo[ithisline],rangehi[ithisline]]
;            dzran = zran[1]-zran[0]
;            ncbdiv = rangencbdiv[ithisline]
;            ncbdiv = ncbdiv[0]
;            hasrange = 1
;         endif
;      endif
;;     otherwise set it automagically.
;      if ~hasrange then begin
;         zran = [min(map[igd]),max(map[igd])]
;         divarr = ifsf_cbdiv(zran,100d,ncbdivmax)
;         ncbdiv = divarr[0]
;         dzran = zran[1]-zran[0]
;      endif
;
;;     replace bad points with "bad"
;      map[ibd] = bad
;
;      mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;         min=zran[0],max=zran[1])
;      cgloadct,74,/reverse
;      cgimage,mapscl,/keep,pos=pos[*,5],opos=truepos,$
;         noerase=i ne 0,missing_value=bad,missing_index=255,$
;         missing_color='white'
;      cgplot,[0],xsty=5,ysty=5,position=truepos,$
;         /nodata,/noerase,title='v$\downmax$(NaD em, km/s)'
;      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
;         (dzran - zran[1]),format=cbform)
;      cgcolorbar,position=cbpos,divisions=ncbdiv,$
;         ticknames=ticknames,/ver,/right,charsize=0.6
;
;      endif
;
;      cgps_close
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; velocity (fits)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;      nx = 3
;      nyabs=0
;      if maxnadabsncomp_act gt 0 then begin
;         if donadabsonecomp then nyabs=1
;         if donadabsmulticomp then nyabs+=maxnadabsncomp_act
;      endif
;      nyem=0
;      if maxnademncomp_act gt 0 then begin
;         if donademonecomp then nyem=1
;         if donademmulticomp then nyem+=maxnademncomp_act
;      endif
;      ny = nyabs+nyem
;
;      cgps_open,initdat.mapdir+initdat.label+'NaDfitvel.eps',charsize=1,/encap,$
;                /inches,xs=plotquantum*nx,ys=plotquantum*ny*aspectrat,/qui
;
;      pos = cglayout([nx,ny],ixmar=[2,3],iymar=[2,2],oxmar=[2,0],oymar=[0,2],$
;                     xgap=0,ygap=0,unit=!D.X_PX_CM/3.0)
;      cbform = '(I0)'
;      abslab = ['Abs (1 Comp)']
;      if nyabs gt 1 then $
;         for i=1,nyabs-1 do $
;            abslab = [abslab,string('Abs (',i,' of ',nyabs-1,' Comp)',$
;                                    format=('(A0,I0,A0,I0,A0)'))]
;      emlab = ['Em (1 Comp)']
;      if nyem gt 1 then $
;         for i=1,nyem-1 do $
;            emlab = [emlab,string('Em (',i,' of ',nyem-1,' Comp)',$
;                                  format=('(A0,I0,A0,I0,A0)'))]
;
;;
;;     ABSORPTION
;;
;      for i=0,nyabs-1 do begin
;;        sigma
;         map = nadabssig[*,*,i]
;;         igd = where(map ne bad)
;;         ibd = where(map eq bad)
;         igd_thiscomp = where(map ne bad,ctgd_thiscomp)
;         ibd_thiscomp = where(map eq bad,ctbd_thiscomp)
;         if ctgd_thiscomp gt 0 then $
;            igd = cgsetintersection(igd_nadfitabsweq,igd_thiscomp) $
;         else igd=igd_nadfitabsweq
;         if ctbd_thiscomp gt 0 then $
;            ibd = cgsetunion(ibd_nadfitabsweq,ibd_thiscomp) $
;         else ibd = ibd_nadfitabsweq
;         map[ibd] = bad
;;        Set up range
;         ranlab = 'fitc'+string(i,format='(I0)')+'vsig'
;         if hasrangefile then begin
;            ithisline = where(rangeline eq 'NaDabs' AND $
;                              rangequant eq ranlab,ctthisline)
;            if ctthisline eq 1 then auto=0b
;         endif else auto=1b
;         plotdat = $
;            ifsf_plotrange(auto=auto,$
;                           mapgd=map[igd],divinit=100d,ncbdivmax=ncbdivmax,$
;                           rline=rangeline,matline='NaDabs',$
;                           rquant=rangequant,matquant=ranlab,$
;                           rncbdiv=rangencbdiv,rlo=rangelo,rhi=rangehi)
;         mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;                         min=plotdat[0],max=plotdat[1])
;         cgloadct,65,/reverse
;         cgimage,mapscl,/keep,pos=pos[*,0+i*nx],opos=truepos,$
;                 noerase=i ne 0,missing_value=bad,missing_index=255,$
;                 missing_color='white'
;         cgplot,[0],xsty=5,ysty=5,position=truepos,$
;                /nodata,/noerase
;         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;         if i eq 0 then begin
;            cgtext,(xran_kpc[0]+xran_kpc[1])/2d,$
;                   yran_kpc[1]+0.1*(yran_kpc[1]-yran_kpc[0]),$
;                  '$\sigma$ (km/s)',chars=1.25,align=0.5
;            ifsf_plotcompass,xarr_kpc,yarr_kpc
;         endif
;         cgtext,xran_kpc[0]-0.17*(xran_kpc[1]-xran_kpc[0]),$
;                (yran_kpc[0]+yran_kpc[1])/2d,$
;                abslab[i],align=0.5,orient=90d,chars=1.25
;         cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
;                            (plotdat[2] - plotdat[1]),format=cbform)
;         cgcolorbar,position=cbpos,divisions=plotdat[3],$
;                    ticknames=ticknames,/ver,/right,charsize=0.6
;         
;;        v50
;         map = nadabsvel[*,*,i]
;;         igd = where(map ne bad)
;;         ibd = where(map eq bad)
;         igd_thiscomp = where(map ne bad,ctgd_thiscomp)
;         ibd_thiscomp = where(map eq bad,ctbd_thiscomp)
;         if ctgd_thiscomp gt 0 then $
;            igd = cgsetintersection(igd_nadfitabsweq,igd_thiscomp) $
;         else igd=igd_nadfitabsweq
;         if ctbd_thiscomp gt 0 then $
;            ibd = cgsetunion(ibd_nadfitabsweq,ibd_thiscomp) $
;         else ibd = ibd_nadfitabsweq
;         map[ibd] = bad
;;        Set up range
;         ranlab = 'fitc'+string(i,format='(I0)')+'v%50'
;         if hasrangefile then begin
;            ithisline = where(rangeline eq 'NaDabs' AND $
;                              rangequant eq ranlab,ctthisline)
;            if ctthisline eq 1 then auto=0b
;         endif else auto=1b
;         plotdat = $
;            ifsf_plotrange(auto=auto,$
;                           mapgd=map[igd],divinit=200d,ncbdivmax=ncbdivmax,$
;                           rline=rangeline,matline='NaDabs',$
;                           rcomp=rangecomp,matquant=ranlab,$
;                           rncbdiv=rangencbdiv,rlo=rangelo,rhi=rangehi)
;         mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;                         min=plotdat[0],max=plotdat[1])
;;         cgloadct,22,/brewer,/reverse
;         cgloadct,74,/reverse
;         cgimage,mapscl,/keep,pos=pos[*,1+i*nx],opos=truepos,$
;                 /noerase,missing_value=bad,missing_index=255,$
;                 missing_color='white'
;         cgplot,[0],xsty=5,ysty=5,position=truepos,$
;                /nodata,/noerase
;         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;         if i eq 0 then $
;            cgtext,(xran_kpc[0]+xran_kpc[1])/2d,$
;                   yran_kpc[1]+0.1*(yran_kpc[1]-yran_kpc[0]),$
;                  'v$\down50$ (km/s)',chars=1.25,align=0.5
;         cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
;                            (plotdat[2] - plotdat[1]),format=cbform)
;         cgcolorbar,position=cbpos,divisions=plotdat[3],$
;                    ticknames=ticknames,/ver,/right,charsize=0.6
;;        v98
;         map = nadabsv98[*,*,i]
;;         igd = where(map ne bad)
;;         ibd = where(map eq bad)
;         igd_thiscomp = where(map ne bad,ctgd_thiscomp)
;         ibd_thiscomp = where(map eq bad,ctbd_thiscomp)
;         if ctgd_thiscomp gt 0 then $
;            igd = cgsetintersection(igd_nadfitabsweq,igd_thiscomp) $
;         else igd=igd_nadfitabsweq
;         if ctbd_thiscomp gt 0 then $
;            ibd = cgsetunion(ibd_nadfitabsweq,ibd_thiscomp) $
;         else ibd = ibd_nadfitabsweq
;         map[ibd] = bad
;;        Set up range
;         ranlab = 'fitc'+string(i,format='(I0)')+'v%98'
;         if hasrangefile then begin
;            ithisline = where(rangeline eq 'NaDabs' AND $
;                              rangequant eq ranlab,ctthisline)
;            if ctthisline eq 1 then auto=0b
;         endif else auto=1b
;         plotdat = $
;            ifsf_plotrange(auto=auto,$
;                           mapgd=map[igd],divinit=200d,ncbdivmax=ncbdivmax,$
;                           rline=rangeline,matline='NaDabs',$
;                           rcomp=rangecomp,matquant=ranlab,$
;                           rncbdiv=rangencbdiv,rlo=rangelo,rhi=rangehi)
;         mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;            min=plotdat[0],max=plotdat[1])
;;         cgloadct,22,/brewer,/reverse
;         cgloadct,74,/reverse
;         cgimage,mapscl,/keep,pos=pos[*,2+i*nx],opos=truepos,$
;                 /noerase,missing_value=bad,missing_index=255,$
;                 missing_color='white'
;         cgplot,[0],xsty=5,ysty=5,position=truepos,$
;                /nodata,/noerase
;         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;         if i eq 0 then $
;            cgtext,(xran_kpc[0]+xran_kpc[1])/2d,$
;                   yran_kpc[1]+0.1*(yran_kpc[1]-yran_kpc[0]),$
;                  'v$\down98$ (km/s)',chars=1.25,align=0.5
;         cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
;                            (plotdat[2] - plotdat[1]),format=cbform)
;         cgcolorbar,position=cbpos,divisions=plotdat[3],$
;                    ticknames=ticknames,/ver,/right,charsize=0.6
;;
;      endfor
;;
;;     EMISSION
;;
;      for i=0,nyem-1 do begin
;;        sigma
;         map = nademsig[*,*,i]
;;         igd = where(map ne bad)
;;         ibd = where(map eq bad)
;         igd_thiscomp = where(map ne bad)
;         ibd_thiscomp = where(map eq bad)
;         igd = cgsetintersection(igd_nadfitemweq,igd_thiscomp)
;         ibd = cgsetunion(ibd_nadfitemweq,ibd_thiscomp)
;         map[ibd] = bad
;         plotdat = ifsf_plotrange(/auto,mapgd=map[igd],divinit=200d,$
;                                  ncbdivmax=ncbdivmax)
;         mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;            min=plotdat[0],max=plotdat[1])
;         cgloadct,65,/reverse
;         cgimage,mapscl,/keep,pos=pos[*,nyabs*nx+i*nx],opos=truepos,$
;                 /noerase,missing_value=bad,missing_index=255,$
;                 missing_color='white'
;         cgplot,[0],xsty=5,ysty=5,position=truepos,$
;                /nodata,/noerase
;         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;         cgtext,xran_kpc[0]-0.17*(xran_kpc[1]-xran_kpc[0]),$
;                (yran_kpc[0]+yran_kpc[1])/2d,$
;                emlab[i],align=0.5,orient=90d,chars=1.25
;         cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
;                            (plotdat[2] - plotdat[1]),format=cbform)
;         cgcolorbar,position=cbpos,divisions=plotdat[3],$
;                    ticknames=ticknames,/ver,/right,charsize=0.6
;;        v50
;         map = nademvel[*,*,i]
;;         igd = where(map ne bad)
;;         ibd = where(map eq bad)
;         igd_thiscomp = where(map ne bad)
;         ibd_thiscomp = where(map eq bad)
;         igd = cgsetintersection(igd_nadfitemweq,igd_thiscomp)
;         ibd = cgsetunion(ibd_nadfitemweq,ibd_thiscomp)
;         map[ibd] = bad
;         plotdat = ifsf_plotrange(/auto,mapgd=map[igd],divinit=200d,$
;                                  ncbdivmax=ncbdivmax)
;         mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;            min=plotdat[0],max=plotdat[1])
;         cgloadct,22,/brewer,/reverse
;;         cgloadct,65,/reverse
;         cgimage,mapscl,/keep,pos=pos[*,1+nyabs*nx+i*nx],opos=truepos,$
;                 /noerase,missing_value=bad,missing_index=255,$
;                 missing_color='white'
;         cgplot,[0],xsty=5,ysty=5,position=truepos,$
;                /nodata,/noerase
;         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;         cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
;                            (plotdat[2] - plotdat[1]),format=cbform)
;         cgcolorbar,position=cbpos,divisions=plotdat[3],$
;                    ticknames=ticknames,/ver,/right,charsize=0.6
;;        v98
;         map = nademv98[*,*,i]
;;         igd = where(map ne bad)
;;         ibd = where(map eq bad)
;         igd_thiscomp = where(map ne bad)
;         ibd_thiscomp = where(map eq bad)
;         igd = cgsetintersection(igd_nadfitemweq,igd_thiscomp)
;         ibd = cgsetunion(ibd_nadfitemweq,ibd_thiscomp)
;         map[ibd] = bad
;         plotdat = ifsf_plotrange(/auto,mapgd=map[igd],divinit=200d,$
;                                  ncbdivmax=ncbdivmax)
;         mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;            min=plotdat[0],max=plotdat[1])
;         cgloadct,22,/brewer,/reverse
;;         cgloadct,65,/reverse
;         cgimage,mapscl,/keep,pos=pos[*,2+nyabs*nx+i*nx],opos=truepos,$
;                 /noerase,missing_value=bad,missing_index=255,$
;                 missing_color='white'
;         cgplot,[0],xsty=5,ysty=5,position=truepos,$
;                /nodata,/noerase
;         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;         cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
;                            (plotdat[2] - plotdat[1]),format=cbform)
;         cgcolorbar,position=cbpos,divisions=plotdat[3],$
;                    ticknames=ticknames,/ver,/right,charsize=0.6
;;
;      endfor
;
;      cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; tau and C_f
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;      nx = 2
;      ny = 1
;;      if donadabsonecomp then ny=1
;      if donadabsmulticomp then ny+=maxnadabsncomp_act
;
;      cgps_open,initdat.mapdir+initdat.label+'NaDfit_taucf.eps',charsize=1,/encap,$
;                /inches,xs=plotquantum*nx,ys=plotquantum*ny*aspectrat,/qui
;
;      pos = cglayout([nx,ny],ixmar=[2,3],iymar=[2,2],oxmar=[2,0],oymar=[0,2],$
;                     xgap=0,ygap=0,unit=!D.X_PX_CM/3.0)
;      abslab = ['Abs (1 Comp)']
;      if nyabs gt 1 then $
;         for i=1,nyabs do $
;            abslab = [abslab,string('Abs (',i,' of ',nyabs,' Comp)',$
;                                    format=('(A0,I0,A0,I0,A0)'))]
;
;      for i=0,ny-1 do begin
;;        tau
;         map = nadabstau[*,*,i]
;;         igd = where(map ne bad)
;;         ibd = where(map eq bad)
;         igd_thiscomp = where(map ne bad)
;         ibd_thiscomp = where(map eq bad)
;         igd = cgsetintersection(igd_nadfitabsweq,igd_thiscomp)
;         ibd = cgsetunion(ibd_nadfitabsweq,ibd_thiscomp)
;         map[ibd]=bad
;         plotdat = ifsf_plotrange(/auto,mapgd=map[igd],divinit=0.5d,$
;                                  ncbdivmax=ncbdivmax)
;         mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;                         min=plotdat[0],max=plotdat[1])
;         cgloadct,65,/reverse
;         cgimage,mapscl,/keep,pos=pos[*,0+i*nx],opos=truepos,$
;                 noerase=i ne 0,missing_value=bad,missing_index=255,$
;                 missing_color='white'
;         cgplot,[0],xsty=5,ysty=5,position=truepos,$
;                /nodata,/noerase
;         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;         if i eq 0 then begin
;            cgtext,(xran_kpc[0]+xran_kpc[1])/2d,$
;                   yran_kpc[1]+0.1*(yran_kpc[1]-yran_kpc[0]),$
;                  '$\tau$',chars=1.25,align=0.5
;            ifsf_plotcompass,xarr_kpc,yarr_kpc
;         endif
;         cgtext,xran_kpc[0]-0.17*(xran_kpc[1]-xran_kpc[0]),$
;                (yran_kpc[0]+yran_kpc[1])/2d,$
;                abslab[i],align=0.5,orient=90d,chars=1.25
;         cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
;                            (plotdat[2] - plotdat[1]),format=cbform)
;         cbform = '(I0)'
;         cgcolorbar,position=cbpos,divisions=plotdat[3],$
;                    ticknames=ticknames,/ver,/right,charsize=0.6
;         
;;        C_f
;         map = nadabscf[*,*,i]
;;         igd = where(map ne bad)
;;         ibd = where(map eq bad)
;         igd_thiscomp = where(map ne bad)
;         ibd_thiscomp = where(map eq bad)
;         igd = cgsetintersection(igd_nadfitabsweq,igd_thiscomp)
;         ibd = cgsetunion(ibd_nadfitabsweq,ibd_thiscomp)
;         map[ibd]=bad
;         plotdat = ifsf_plotrange(/auto,mapgd=map[igd],divinit=0.2d,$
;                                  ncbdivmax=ncbdivmax)
;         mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;                         min=plotdat[0],max=plotdat[1])
;         cgloadct,65,/reverse
;         cgimage,mapscl,/keep,pos=pos[*,1+i*nx],opos=truepos,$
;                 /noerase,missing_value=bad,missing_index=255,$
;                 missing_color='white'
;         cgplot,[0],xsty=5,ysty=5,position=truepos,$
;                /nodata,/noerase
;         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;         if i eq 0 then $
;            cgtext,(xran_kpc[0]+xran_kpc[1])/2d,$
;                   yran_kpc[1]+0.1*(yran_kpc[1]-yran_kpc[0]),$
;                  'C$\downf$',chars=1.25,align=0.5
;         cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
;                            (plotdat[2] - plotdat[1]),format=cbform)
;         cbform = '(D0.1)'
;         cgcolorbar,position=cbpos,divisions=plotdat[3],$
;                    ticknames=ticknames,/ver,/right,charsize=0.6
;      endfor
;
;      cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; end NaD plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   endif



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Outflow maps
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;   if tag_exist(initmaps,'compof') AND $
;    ~ tag_exist(initmaps,'noemlinfit') then begin
;
;;     Quantities to plot
;      if tag_exist(initmaps,'cvdf_oftags') then tags=initmaps.cvdf_oftags $
;      else tags=['f','sig','v50','v84','v98']
;      if tag_exist(initmaps,'cvdf_oftagtypes') then $
;         tagtypes=initmaps.cvdf_oftagtypes $
;      else tagtypes=['flux','vel','vel','vel','vel']
;      if tag_exist(initmaps,'cvdf_oftitles') then $
;         titles=initmaps.cvdf_oftitles $
;      else titles = ['F$\downtot$','$\sigma$','v$\down50$',$
;                     'v$\down84$','v$\down98$']
;
;;     Size of plot grid
;      if tag_exist(initmaps,'cvdf_ofnpx') then npx=initmaps.cvdf_ofnpx $
;      else npx=3
;      if tag_exist(initmaps,'cvdf_ofnpy') then npy=initmaps.cvdf_ofnpy $
;      else npy=2
;      
;;     Loop through emission lines
;      foreach line,lines_with_doublets do begin
;
;         linelab = ifsf_linesyntax(line)
;         cgps_open,initdat.mapdir+initdat.label+linelab+'_of_c.eps',charsize=1,/encap,$
;                   /inches,xs=plotquantum*npx,ys=plotquantum*npy*aspectrat,/qui
;         pos = cglayout([npx,npy],ixmar=[3,3],iymar=[3,3],oxmar=[0,0],oymar=[0,0],$
;                        xgap=0,ygap=0,unit=!D.X_PX_CM/3.0)
;
;;        loop through plot types
;         for j=0,n_elements(tags)-1 do begin
;
;            iplot = j ; plot index
;;           Get map and scale
;            itag = where(strcmp(ofpars_tags,tags[j],/fold_case) eq 1,cttag)
;            if cttag ne 0 then begin
;
;            map = ofpars[line].(itag)
;            ibd = where(map eq bad AND ~ finite(map),ctbd)
;            inan = where(~finite(map),ctnan)
;            igd = where(map ne bad AND map ne 0 AND finite(map),ctgd)
;
;            title=titles[j]
;            
;            if ctgd gt 0 then begin
;               
;               if tagtypes[j] eq 'flux' then begin
; 
;                  cbform = '(D0.1)' ; colorbar syntax
;               
;                  if tag_exist(initmaps,'fluxfactor') then $
;                     map[igd] *= initmaps.fluxfactor
;                                 
;                  zran=[0,1]
;                  dzran = 1
;                  ncbdiv = 5
;                  zmax_flux = max(map[igd])
;
;                  if hasrangefile then begin
;                     ithisline = where(line eq rangeline AND $
;                                       'of_'+tags[j] eq rangequant,ctthisline)
;                     if ctthisline eq 1 then zmax_flux = rangelo[ithisline]
;                  endif
;                
;                  map[igd] = map[igd]/zmax_flux[0]
;                  title += ' ('+string(zmax_flux,format='(E0.2)')+')'
;                  cgloadct,65,/reverse
;
;               endif else begin
;       
;                  cbform = '(I0)' ; colorbar syntax
;                  hasrange = 0
;                  if hasrangefile then begin
;                     ithisline = where(line eq rangeline AND $
;                                       'of_'+tags[j] eq rangequant,ctthisline)
;                     if ctthisline eq 1 then begin
;                        zran = [rangelo[ithisline],rangehi[ithisline]]
;                        dzran = zran[1]-zran[0]
;                        ncbdiv = rangencbdiv[ithisline]                   
;                        ncbdiv = ncbdiv[0]
;                        hasrange = 1
;                     endif
;                  endif
;                  if ~hasrange then begin
;                     zran = [min(map[igd]),max(map[igd])]
;                     divarr = ifsf_cbdiv(zran,100d,ncbdivmax)
;                     ncbdiv = divarr[0]
;                     dzran = zran[1]-zran[0]
;                  endif
;
;                  cgloadct,74,/reverse
; 
;               endelse
; 
;               if ctnan gt 0 then map[inan] = bad
;               mapscl = bytscl(rebin(map,dx*samplefac,dy*samplefac,/sample),$
;                               min=zran[0],max=zran[1])
;
;;              Plot image
;               cgimage,mapscl,/keep,pos=pos[*,iplot],opos=truepos,$
;                       noerase=iplot ne 0,missing_value=bad,missing_index=255,$
;                       missing_color='white'
;               cgplot,[0],xsty=5,ysty=5,position=truepos,$
;                         /nodata,/noerase,title=title
;               ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;;              Colorbar
;               cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
;               ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
;                                  (dzran - zran[1]),format=cbform)
;               cgcolorbar,position=cbpos,divisions=ncbdiv,$
;                          ticknames=ticknames,/ver,/right,charsize=0.6
;
;            endif
;            endif
;            
;         endfor
; 
;         cgps_close
;
;      endforeach
;
;   endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Compute outflow properties
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   m_nad = 0d
   dmdt_nad = 0d
   p_nad = 0d
   dpdt_nad = 0d
   e_nad = 0d
   dedt_nad = 0d
   Rmax_kpc_nad = 0d

;   m_errlo_arr_nad = !NULL
;   dmdt_errlo_arr_nad = !NULL
;   p_errlo_arr_nad = !NULL
;   dpdt_errlo_arr_nad = !NULL
;   e_errlo_arr_nad = !NULL
;   dedt_errlo_arr_nad = !NULL
;   m_errhi_arr_nad = !NULL
;   dmdt_errhi_arr_nad = !NULL
;   p_errhi_arr_nad = !NULL
;   dpdt_errhi_arr_nad = !NULL
;   e_errhi_arr_nad = !NULL
;   dedt_errhi_arr_nad = !NULL
   m_err_nad = dblarr(2)
   dmdt_err_nad = dblarr(2)
   p_err_nad = dblarr(2)
   dpdt_err_nad = dblarr(2)
   e_err_nad = dblarr(2)
   dedt_err_nad = dblarr(2)

   m_ha = 0d
   dmdt_ha = 0d
   p_ha = 0d
   dpdt_ha = 0d
   e_ha = 0d
   dedt_ha = 0d
   Rmax_kpc_ha = 0d

   m_err_ha = 0d
   dmdt_err_ha = 0d
   p_err_ha = 0d
   dpdt_err_ha = 0d
   e_err_ha = 0d
   dedt_err_ha = 0d


   ;;;;;;;;;;;
   ; NEUTRAL ;
   ;;;;;;;;;;;

;  Sum over outflowing components to get outflow properties

   if tag_exist(initdat,'donad') AND tag_exist(initmaps,'outflow_abs') then begin

      ;;;;;;;;;;;;;;
      ; CRFW model ;
      ;;;;;;;;;;;;;;
  
      Rcm_nad = 3.08567802d18 * initmaps.outflow_abs.R * 1d3
;     radius of wind in cm
      Rpix_nad = (initmaps.outflow_abs.R / kpc_per_as)/initdat.platescale
;     radius of wind in spaxels
      map_dphi_sintheta_nad = 1/Rpix_nad
      map_dtheta_nad = asin((map_r+0.5)/Rpix_nad) - asin((map_r-0.5)/Rpix_nad)
      map_domega_nad = map_dphi_sintheta_nad*map_dtheta_nad
      map_costheta_nad = cos(asin(map_r/Rpix_nad))

;     Cycle through components
      for i=0,maxnadabsncomp_act do begin
         cnhcf =  nadabscnhcf[*,*,i]
         vel = nadabsvel[*,*,i]*1d5
         sig = nadabssig[*,*,i]*1d5
         errcnhcf =  reform(errnadabscnhcf[*,*,i,*],dx,dy,2)
         errvel = rebin(errnadabsvel[*,*,i,*],dx,dy,2)*1d5
         errsig = rebin(errnadabssig[*,*,i,*],dx,dy,2)*1d5
         igdof = where(cnhcf ne bad AND cnhcf gt 0 $
                       AND cnhcf gt errcnhcf[*,*,0] AND vel le 0,ctgd)
         if ctgd gt 0 then begin

;           Make sure velocity error on low end isn't greater than velocity
            ibdvel = where(errvel[*,*,0] gt abs(vel),ctbdvel)
            if ctbdvel gt 0 then errvel[ibdvel] = abs(vel[ibdvel])

;           Calculate max radius
            Rmax_kpc_nad_c = max(map_rkpc_ifs[igdof])
            Rmax_kpc_nad = Rmax_kpc_nad_c > Rmax_kpc_nad ? Rmax_kpc_nad_c : Rmax_kpc_nad

;           Now discard above specified radius
            igdof = where(cnhcf ne bad AND $
                          cnhcf gt 0 AND $
                          cnhcf gt errcnhcf[*,*,0] AND $
                          vel le 0 AND $
                          map_rkpc_ifs lt initmaps.outflow_abs.R,ctgd)

            map_vrad_comp_nad = dblarr(dx,dy)
            map_sig_comp_nad = dblarr(dx,dy)
            map_m_comp_nad = dblarr(dx,dy)
            map_dmdt_comp_nad = dblarr(dx,dy)
            map_p_comp_nad = dblarr(dx,dy)
            map_dpdt_comp_nad = dblarr(dx,dy)
            map_e_comp_nad = dblarr(dx,dy)
            map_dedt_comp_nad = dblarr(dx,dy)
            map_vrad_comp_err_nad = dblarr(dx,dy,2)
            map_sig_comp_err_nad = dblarr(dx,dy,2)
            map_m_comp_err_nad = dblarr(dx,dy,2)
            map_dmdt_comp_err_nad = dblarr(dx,dy,2)
            map_p_comp_err_nad = dblarr(dx,dy,2)
            map_dpdt_comp_err_nad = dblarr(dx,dy,2)
            map_e_comp_err_nad = dblarr(dx,dy,2)
            map_dedt_comp_err_nad = dblarr(dx,dy,2)

            map_vrad_comp_nad[igdof] = $
               abs(vel[igdof]) / map_costheta_nad[igdof]
            map_sig_comp_nad[igdof] = sig[igdof]
            map_m_comp_nad[igdof] = $
               mumpsm*Rcm_nad^2d*cnhcf[igdof]*map_domega_nad[igdof]
            map_dmdt_comp_nad[igdof] = $
               mumpsm*Rcm_nad*cnhcf[igdof]*map_vrad_comp_nad[igdof]*speryr*$
               map_domega_nad[igdof]
            map_p_comp_nad[igdof] = $
               map_m_comp_nad[igdof] * msun * map_vrad_comp_nad[igdof]
                  ; in dyne s
            map_dpdt_comp_nad[igdof] = $
               map_dmdt_comp_nad[igdof] * msun / speryr * $
               map_vrad_comp_nad[igdof] * c_cms / lsun
                  ; dp/dt*c, in Lsun
            map_e_comp_nad[igdof] = $
               map_m_comp_nad[igdof] * msun * $
               (0.5d*map_vrad_comp_nad[igdof]^2d + $
               1.5d*map_sig_comp_nad[igdof]^2d)
                  ; in erg
            map_dedt_comp_nad[igdof] = $
               map_dmdt_comp_nad[igdof] * msun / speryr * $
               (0.5d*map_vrad_comp_nad[igdof]^2d + $
               1.5d*map_sig_comp_nad[igdof]^2d)
                  ; in erg/s

            for j=0,1 do begin
               xyind = array_indices(cnhcf,igdof)
               xigdof = xyind[0,*]
               yigdof = xyind[1,*]
               for k=0,ctgd-1 do begin
                  map_vrad_comp_err_nad[xigdof[k],yigdof[k],j] = $
                    abs(errvel[xigdof[k],yigdof[k],j]) / map_costheta_nad[igdof[k]]
                  map_sig_comp_err_nad[xigdof[k],yigdof[k],j] = errsig[xigdof[k],yigdof[k],j]
                  map_m_comp_err_nad[xigdof[k],yigdof[k],j] = $
                     map_m_comp_nad[xigdof[k],yigdof[k]] * $
                     errcnhcf[xigdof[k],yigdof[k],j]/cnhcf[xigdof[k],yigdof[k]]
                  map_dmdt_comp_err_nad[xigdof[k],yigdof[k],j] = $
                     map_dmdt_comp_nad[igdof[k]] * $
                     sqrt((errcnhcf[xigdof[k],yigdof[k],j]/cnhcf[igdof[k]])^2d + $
                          (map_vrad_comp_err_nad[xigdof[k],yigdof[k],j]/$
                           map_vrad_comp_nad[igdof[k]])^2d)
                  map_p_comp_err_nad[xigdof[k],yigdof[k],j] = $
                     map_p_comp_nad[igdof[k]] * $
                     sqrt((errcnhcf[xigdof[k],yigdof[k],j]/cnhcf[igdof[k]])^2d + $
                          (map_vrad_comp_err_nad[xigdof[k],yigdof[k],j]/$
                           map_vrad_comp_nad[igdof[k]])^2d)
                  map_dpdt_comp_err_nad[xigdof[k],yigdof[k],j] = $
                     map_dpdt_comp_nad[igdof[k]] * $
                     sqrt((errcnhcf[xigdof[k],yigdof[k],j]/cnhcf[igdof[k]])^2d + $
                          2d*(map_vrad_comp_err_nad[xigdof[k],yigdof[k],j]/$
                           map_vrad_comp_nad[igdof[k]])^2d)
                  map_e_comp_err_nad[xigdof[k],yigdof[k],j] = $
                     sqrt($
                        (map_e_comp_nad[xigdof[k],yigdof[k]]*$
                         map_m_comp_err_nad[xigdof[k],yigdof[k],j]/$
                         map_m_comp_nad[xigdof[k],yigdof[k]])^2d + $
                        (msun*map_m_comp_nad[xigdof[k],yigdof[k]]*$
                         map_vrad_comp_nad[xigdof[k],yigdof[k]]*$
                         map_vrad_comp_err_nad[xigdof[k],yigdof[k],j])^2d + $
                        (3d*msun*map_m_comp_nad[xigdof[k],yigdof[k]]*$
                         map_sig_comp_nad[xigdof[k],yigdof[k]]*$
                         map_sig_comp_err_nad[xigdof[k],yigdof[k],j])^2d)
                  map_dedt_comp_err_nad[xigdof[k],yigdof[k],j] = $
                     sqrt($
                        (map_dedt_comp_nad[xigdof[k],yigdof[k]]*$
                         map_dmdt_comp_err_nad[xigdof[k],yigdof[k],j]/$
                         map_dmdt_comp_nad[xigdof[k],yigdof[k]])^2d + $
                        (msun/speryr*map_dmdt_comp_nad[xigdof[k],yigdof[k]]*$
                         map_vrad_comp_nad[xigdof[k],yigdof[k]]*$
                         map_vrad_comp_err_nad[xigdof[k],yigdof[k],j])^2d + $
                        (3d*msun/speryr*map_dmdt_comp_nad[xigdof[k],yigdof[k]]*$
                         map_sig_comp_nad[xigdof[k],yigdof[k]]*$
                         map_sig_comp_err_nad[xigdof[k],yigdof[k],j])^2d)
               endfor
            endfor

            m_nad += 2d*total(map_m_comp_nad)
            dmdt_nad += 2d*total(map_dmdt_comp_nad)
            p_nad += 2d*total(map_p_comp_nad)
            dpdt_nad += 2d*total(map_dpdt_comp_nad)
            e_nad += 2d*total(map_e_comp_nad)
            dedt_nad += 2d*total(map_dedt_comp_nad)

;            Errors using median relative error across spaxels
;            tmperrlo = map_m_comp_err_nad[*,*,0]
;            tmperrhi = map_m_comp_err_nad[*,*,1]            
;            m_errlo_arr_nad = [m_errlo_arr_nad,tmperrlo[igdof]/map_m_comp_nad[igdof]]
;            m_errhi_arr_nad = [m_errhi_arr_nad,tmperrhi[igdof]/map_m_comp_nad[igdof]]
;
;            tmperrlo = map_dmdt_comp_err_nad[*,*,0]
;            tmperrhi = map_dmdt_comp_err_nad[*,*,1]
;            dmdt_errlo_arr_nad = [dmdt_errlo_arr_nad,tmperrlo[igdof]/map_dmdt_comp_nad[igdof]]
;            dmdt_errhi_arr_nad = [dmdt_errhi_arr_nad,tmperrhi[igdof]/map_dmdt_comp_nad[igdof]]
;
;            tmperrlo = map_p_comp_err_nad[*,*,0]
;            tmperrhi = map_p_comp_err_nad[*,*,1]
;            p_errlo_arr_nad = [p_errlo_arr_nad,tmperrlo[igdof]/map_p_comp_nad[igdof]]
;            p_errhi_arr_nad = [p_errhi_arr_nad,tmperrhi[igdof]/map_p_comp_nad[igdof]]
;
;            tmperrlo = map_dpdt_comp_err_nad[*,*,0]
;            tmperrhi = map_dpdt_comp_err_nad[*,*,1]
;            dpdt_errlo_arr_nad = [dpdt_errlo_arr_nad,tmperrlo[igdof]/map_dpdt_comp_nad[igdof]]
;            dpdt_errhi_arr_nad = [dpdt_errhi_arr_nad,tmperrhi[igdof]/map_dpdt_comp_nad[igdof]]

;           fix low errors that are too big! probably should do this earlier, but OK
;           for now
            map_m_comp_err_nad_lo = map_m_comp_err_nad[*,*,0]
            ibdlo = where(map_m_comp_err_nad_lo gt map_m_comp_nad,ctbdlo)
            if ctbdlo gt 0 then begin
               map_m_comp_err_nad_lo[ibdlo] = map_m_comp_nad[ibdlo]
               map_m_comp_err_nad[*,*,0] = map_m_comp_err_nad_lo
            endif
            map_dmdt_comp_err_nad_lo = map_dmdt_comp_err_nad[*,*,0]
            ibdlo = where(map_dmdt_comp_err_nad_lo gt map_dmdt_comp_nad,ctbdlo)
            if ctbdlo gt 0 then begin
               map_dmdt_comp_err_nad_lo[ibdlo] = map_dmdt_comp_nad[ibdlo]
               map_dmdt_comp_err_nad[*,*,0] = map_dmdt_comp_err_nad_lo
            endif
            map_p_comp_err_nad_lo = map_p_comp_err_nad[*,*,0]
            ibdlo = where(map_p_comp_err_nad_lo gt map_p_comp_nad,ctbdlo)
            if ctbdlo gt 0 then begin
               map_p_comp_err_nad_lo[ibdlo] = map_p_comp_nad[ibdlo]
               map_p_comp_err_nad[*,*,0] = map_p_comp_err_nad_lo
            endif
            map_dpdt_comp_err_nad_lo = map_dpdt_comp_err_nad[*,*,0]
            ibdlo = where(map_dpdt_comp_err_nad_lo gt map_dpdt_comp_nad,ctbdlo)
            if ctbdlo gt 0 then begin
               map_dpdt_comp_err_nad_lo[ibdlo] = map_dpdt_comp_nad[ibdlo]
               map_dpdt_comp_err_nad[*,*,0] = map_dpdt_comp_err_nad_lo
            endif
            map_e_comp_err_nad_lo = map_e_comp_err_nad[*,*,0]
            ibdlo = where(map_e_comp_err_nad_lo gt map_e_comp_nad,ctbdlo)
            if ctbdlo gt 0 then begin
               map_e_comp_err_nad_lo[ibdlo] = map_e_comp_nad[ibdlo]
               map_e_comp_err_nad[*,*,0] = map_e_comp_err_nad_lo
            endif
            map_dedt_comp_err_nad_lo = map_dedt_comp_err_nad[*,*,0]
            ibdlo = where(map_dedt_comp_err_nad_lo gt map_dedt_comp_nad,ctbdlo)
            if ctbdlo gt 0 then begin
               map_dedt_comp_err_nad_lo[ibdlo] = map_dedt_comp_nad[ibdlo]
               map_dedt_comp_err_nad[*,*,0] = map_dedt_comp_err_nad_lo
            endif

;           Errors using usual error estimation
            m_err_nad += total(total(map_m_comp_err_nad^2d,1),1)
            dmdt_err_nad += total(total(map_dmdt_comp_err_nad^2d,1),1)
            p_err_nad += total(total(map_p_comp_err_nad^2d,1),1)
            dpdt_err_nad += total(total(map_dpdt_comp_err_nad^2d,1),1)
            e_err_nad += total(total(map_e_comp_err_nad^2d,1),1)
            dedt_err_nad += total(total(map_dedt_comp_err_nad^2d,1),1)


         endif
 
      endfor

      m_err_nad = 2d*sqrt(m_err_nad)      
      dmdt_err_nad = 2d*sqrt(dmdt_err_nad)
      p_err_nad = 2d*sqrt(p_err_nad)
      dpdt_err_nad = 2d*sqrt(dpdt_err_nad)
      e_err_nad = 2d*sqrt(e_err_nad)
      dedt_err_nad = 2d*sqrt(dedt_err_nad)
;      m_err_nad = [median(m_errlo_arr_nad),median(m_errhi_arr_nad)]*m_nad
;      dmdt_err_nad = [median(dmdt_errlo_arr_nad),median(dmdt_errhi_arr_nad)]*dmdt_nad
;      p_err_nad = [median(p_errlo_arr_nad),median(p_errhi_arr_nad)]*p_nad
;      dpdt_err_nad = [median(dpdt_errlo_arr_nad),median(dpdt_errhi_arr_nad)]*dpdt_nad
      
   endif  

   if ~ tag_exist(initdat,'noemlinfit') then begin
    
      if linelist.haskey('Halpha') then ofcompline = 'Halpha' $
      else if linelist.haskey('Hbeta') then ofcompline = 'Hbeta' $
      else begin
         linelistkeys = linelist.keys()
         ofcompline = linelistkeys[0]
         message,'Halpha and Hbeta absent; setting OFCOMPLINE to random line.',/cont
      endelse
   endif else begin
      ofcompline = ''
   endelse
   if tag_exist(initmaps,'outflow_eml') then begin

      ;;;;;;;;;;;
      ; IONIZED ;
      ;;;;;;;;;;;

;;     Use CVDF to get outflow properties. Use Halpha if available, otherwise
;;     Hbeta.
;;     Blueshifted
;      map_fha = ofpars[ofcompline].f
;      igd = where(map_fha ne bad,ctgd)
;;     Convert from flux/arcsec^2 back to straight flux if requested
;      if ctgd gt 0 then map_fha[igd] *= initdat.platescale
;      map_vha = ofpars[ofcompline].v50*1d5
;      map_sigha = ofpars[ofcompline].sig*1d5
;                                ; velocities in cm/s
;
;;     Redshifted
;      mapr_fha = ofpars[ofcompline].fr
;      igd = where(mapr_fha ne bad,ctgd)
;;     Convert from flux/arcsec^2 back to straight flux if requested
;      if ctgd gt 0 then mapr_fha[igd] *= initdat.platescale
;      mapr_vha = ofpars[ofcompline].v50r*1d5
;      mapr_sigha = ofpars[ofcompline].sigr*1d5
;                                ; velocities in cm/s

;     Use second component to get outflow properties. Use Halpha if available, otherwise
;     Hbeta.

;     Blueshifted
      map_fha = emlflx['fc2',ofcompline]
      map_fha_err = emlflxerr['fc2',ofcompline]
      if tag_exist(initmaps,'ebv') then begin
         if tag_exist(initmaps.ebv,'apply') then begin
            map_fha = emlflxcor_pp['fc2',ofcompline]
            map_fha_err = emlflxerrcor_pp['fc2',ofcompline]
         endif
      endif
      map_vha = emlvel['v%50c2',ofcompline]
      map_vha_err = emlvel['v%50c2err',ofcompline]
      map_sigha = emlvel['vsigc2',ofcompline]
      map_sigha_err = emlvel['vsigc2err',ofcompline]
      igd = where(map_fha gt 0d AND map_fha ne bad AND $
                  map_vha lt 0,ctgd)

;     Convert from flux/arcsec^2 back to straight flux
      if ctgd gt 0 then begin
         map_fha[igd] *= initdat.platescale^2d
         map_fha_err[igd] *= initdat.platescale^2d
         Rmax_kpc_ha = max(map_rkpc_ifs[igd])
;        Convert to Halpha flux using Case B
         if ofcompline eq 'Hbeta' then begin
            map_fha[igd] *= 2.86d
            map_fha_err[igd] *= 2.86d
         endif
;        Velocities in cm/s
         map_vha[igd] *= 1d5
         map_vha_err[igd] *= 1d5
         map_sigha[igd] *= 1d5
         map_sigha_err[igd] *= 1d5
      endif

;     Redshifted
      if tag_exist(initmaps.outflow_eml,'usered') then begin
         mapr_fha = emlflx['fc2',ofcompline]
         mapr_fha_err = emlflxerr['fc2',ofcompline]
         if tag_exist(initmaps,'ebv') then begin
            if tag_exist(initmaps.ebv,'apply') then begin
               mapr_fha = emlflxcor_pp['fc2',ofcompline]
               mapr_fha_err = emlflxerrcor_pp['fc2',ofcompline]
            endif
         endif
         mapr_vha = emlvel['v%50c2',ofcompline]
         mapr_vha_err = emlvel['v%50c2err',ofcompline]
         mapr_sigha = emlvel['vsigc2',ofcompline]
         mapr_sigha_err = emlvel['vsigc2err',ofcompline]
         igd = where(mapr_fha gt 0d AND mapr_fha ne bad AND $
                     mapr_vha ge 0,ctgd)
         if ctgd gt 0 then begin
            mapr_fha[igd] *= initdat.platescale^2d
            mapr_fha_err[igd] *= initdat.platescale^2d
            Rmax_kpc_ha_red = max(map_rkpc_ifs[igd])
            Rmax_kpc_ha = Rmax_kpc_ha_red > Rmax_kpc_ha ? Rmax_kpc_ha_red : Rmax_kpc_ha
            if ofcompline eq 'Hbeta' then begin
               mapr_fha[igd] *= 2.86d
               mapr_fha_err[igd] *= 2.86d
            endif
            mapr_vha[igd] *= 1d5
            mapr_vha_err[igd] *= 1d5
            mapr_sigha[igd] *= 1d5
            mapr_sigha_err[igd] *= 1d5
         endif
      endif else begin
         mapr_fha = map_fha*0d
         mapr_vha = map_fha*0d
         mapr_sigha = map_fha*0d
         mapr_fha_err = map_fha*0d
         mapr_vha_err = map_fha*0d
         mapr_sigha_err = map_fha*0d
      endelse



      ;;;;;;;;;;;;;;
      ; CRFW model ;
      ;;;;;;;;;;;;;;
  
;     Blueshifted
      map_l_ha = dblarr(dx,dy)
      map_vrad_ha = dblarr(dx,dy)
      map_sig_ha = dblarr(dx,dy)
      map_m_ha = dblarr(dx,dy)
      map_p_ha = dblarr(dx,dy)
      map_e_ha = dblarr(dx,dy)
      map_dmdt_ha = dblarr(dx,dy)
      map_dpdt_ha = dblarr(dx,dy)
      map_dedt_ha = dblarr(dx,dy)
      map_l_err_ha = dblarr(dx,dy)
      map_vrad_err_ha = dblarr(dx,dy)
      map_sig_err_ha = dblarr(dx,dy)
      map_m_err_ha = dblarr(dx,dy,2)
      map_p_err_ha = dblarr(dx,dy,2)
      map_e_err_ha = dblarr(dx,dy,2)
      map_dmdt_err_ha = dblarr(dx,dy,2)
      map_dpdt_err_ha = dblarr(dx,dy,2)
      map_dedt_err_ha = dblarr(dx,dy,2)

;     Redshifted
      mapr_l_ha = dblarr(dx,dy)
      mapr_vrad_ha = dblarr(dx,dy)
      mapr_sig_ha = dblarr(dx,dy)
      mapr_m_ha = dblarr(dx,dy)
      mapr_p_ha = dblarr(dx,dy)
      mapr_e_ha = dblarr(dx,dy)
      mapr_dmdt_ha = dblarr(dx,dy)
      mapr_dpdt_ha = dblarr(dx,dy)
      mapr_dedt_ha = dblarr(dx,dy)
      mapr_l_err_ha = dblarr(dx,dy)
      mapr_vrad_err_ha = dblarr(dx,dy)
      mapr_sig_err_ha = dblarr(dx,dy)
      mapr_m_err_ha = dblarr(dx,dy,2)
      mapr_p_err_ha = dblarr(dx,dy,2)
      mapr_e_err_ha = dblarr(dx,dy,2)
      mapr_dmdt_err_ha = dblarr(dx,dy,2)
      mapr_dpdt_err_ha = dblarr(dx,dy,2)
      mapr_dedt_err_ha = dblarr(dx,dy,2)

;     Radius of wind in cm
      Rcm_ha = 3.08567802d18 * initmaps.outflow_eml.R * 1d3
;     Radius of wind in spaxels
      Rpix_ha = (initmaps.outflow_eml.R / kpc_per_as)/initdat.platescale
  
; Solid angle subtended by pixel on a spherical wind, in units of
; radians. A solid angle element is sin(theta) * dphi * dtheta, where
; dphi and dtheta are roughly the deprojected angular sides of a
; spaxel, and where phi is the azimuthal angle (in the sky plane) and
; theta is the polar angle measured from the +z-axis (i.e., angle away
; from the line of sight).
;
; How do we compute dphi, the actual (and deprojected) azimuthal angle
; subtended by the pixel?  The size of the pixel, s, divided by the
; projected circle's radius, R, equals the angle in radians. (We use
; the projected radius instead of the spherical radius because phi is
; an azimuthal angle. E.g., if we were to compute the circumference of
; the wind at this projected radius we would use the projected radius,
; not the spherical radius.) Thus dphi = s/R = 1 / (R/s), where R/s is
; the projected wind radius in units of pixels. But R/s = r/s
; sin(theta), where r is the wind's spherical radius. So dphi = 1/(r
; sin(theta)), where r is in units of pixels. So when we multiply by
; sin(theta), we get dphi = 1/r !
  
      map_dtheta_ha = asin((map_r+0.5)/Rpix_ha) - asin((map_r-0.5)/Rpix_ha)
      map_dphi_sintheta_ha = 1/Rpix_ha
      map_domega_ha = map_dphi_sintheta_ha*map_dtheta_ha

; To deproject velocities from los to radial away from center
      map_costheta_ha = cos(asin(map_r/Rpix_ha))
    
; Extra +0.5 is to prevent map_dtheta from being undefined, and will
; introduce an error of only half a pixel. OK if Rpix >> 1.
      igood_ha = where(map_r+0.5 lt Rpix_ha AND $
                       map_fha ne bad AND $
                       map_fha gt 0 AND $
                       map_vha lt 0,ctgood_ha)
      xyind = array_indices(map_r,igood_ha)
      xigdof = xyind[0,*]
      yigdof = xyind[1,*]
      igoodr_ha = where(map_r+0.5 lt Rpix_ha AND $
                        mapr_fha ne bad AND $
                        mapr_fha gt 0 AND $
                        mapr_vha gt 0,ctgoodr_ha)

;     Compute electron density to use
;      elecdenuse = dindgen(dx,dy) + elecden_default
      elecdenuse = dblarr(dx,dy) + elecden_default
      elecdenuse_err = dblarr(dx,dy,2) + $
                       rebin(reform(elecden_err_default,1,1,2),dx,dy,2)
      elecdencomp = 'ftot'
      if ~ elecdenmap.isempty() then begin
         if elecdenmap.haskey(elecdencomp) then begin
            igds2 = where(elecdenmap[elecdencomp] ne bad,ctgds2)
            if ctgds2 gt 0 then begin
               if tag_exist(initmaps.outflow_eml,'usered') then $
                  igdof = cgsetunion(igood_ha,igoodr_ha) $
               else igdof = igood_ha
               igdelecden = cgsetintersection(igdof,igds2)
               igdelecdenxy = array_indices(elecdenuse,igdelecden)
               igdelecdenx = igdelecdenxy[0,*]
               igdelecdeny = igdelecdenxy[1,*]               
               elecdenuse[igdelecden] = 10d^mean(elecdenmap[elecdencomp,igdelecden])
               for k=0,n_elements(igdelecden)-1 do begin
                  elecdenuse_err[igdelecdenx[k],igdelecdeny[k],1] = $
                     elecdenuse[igdelecden[k]] - $
                     10d^mean(elecdenmap[elecdencomp,igdelecden[k]] - $
                              elecdenmap_errlo[elecdencomp,igdelecden[k]])
                  elecdenuse_err[igdelecdenx[k],igdelecdeny[k],0] = $
                     10d^mean(elecdenmap[elecdencomp,igdelecden[k]] + $
                              elecdenmap_errhi[elecdencomp,igdelecden[k]]) $
                     - elecdenuse[igdelecden[k]]
               endfor
               if tag_exist(initmaps.outflow_eml,'use_spaxel_elecden') then begin
                  elecdenuse[igdelecden] = 10d^elecdenmap[elecdencomp,igdelecden]
                  for k=0,n_elements(igdelecden)-1 do begin
                     elecdenuse_err[igdelecdenx[k],igdelecdeny[k],1] = $
                        elecdenuse[igdelecden[k]] - $
                        10d^(elecdenmap[elecdencomp,igdelecden[k]] - $
                             elecdenmap_errlo[elecdencomp,igdelecden[k]])
                     elecdenuse_err[igdelecdenx[k],igdelecdeny[k],0] = $
                        10d^(elecdenmap[elecdencomp,igdelecden[k]] + $
                             elecdenmap_errhi[elecdencomp,igdelecden[k]]) - $
                        elecdenuse[igdelecden[k]]
                  endfor
               endif

            endif
         endif
      endif
      


; blueshifted
      map_vrad_ha[igood_ha] = abs(map_vha[igood_ha]) / $
                                  map_costheta_ha[igood_ha]
; 1d-7 converts from ergs/s to W
      map_l_ha[igood_ha] = $
         drt_linelum(map_fha[igood_ha]*initdat.fluxunits*1d-3,ldist,/ergs)
      map_m_ha[igood_ha] = mumpsm * map_l_ha[igood_ha] / $
                           volemis / elecdenuse[igood_ha]
                                ; in units of msun
; The following formula comes from dividing dM/dt^avg_thin by M_thin
; (eq. 7 and 8 in RVS05b). Note that this is basically equivalent to
; dividing each pixel by its own dynamical time, tdyn ~ R/v.
      map_dmdt_ha[igood_ha] = map_m_ha[igood_ha] * $
                              map_vrad_ha[igood_ha] * speryr / Rcm_ha
                                ; in units of msun/yr
      map_p_ha[igood_ha] = map_m_ha[igood_ha] * msun * $
                           map_vrad_ha[igood_ha]
                                ; in units of dyne s
      map_dpdt_ha[igood_ha] = map_dmdt_ha[igood_ha] * msun / speryr * $
                              map_vrad_ha[igood_ha] * c_cms / lsun
                                ; dp/dt*c, in Lsun
      map_e_ha[igood_ha] = map_m_ha[igood_ha] * msun * $
                          (0.5d*map_vrad_ha[igood_ha]^2d + $
                           1.5d*map_sig_ha[igood_ha]^2d)
                                ; in erg
      map_dedt_ha[igood_ha] = map_dmdt_ha[igood_ha] * msun / speryr * $
                              (0.5d*map_vrad_ha[igood_ha]^2d + $
                               1.5d*map_sig_ha[igood_ha]^2d)
                                ; in erg/s
                                
      map_vrad_err_ha[igood_ha] = abs(map_vha_err[igood_ha]) / $
                                      map_costheta_ha[igood_ha]

      map_l_err_ha[igood_ha] = $
         drt_linelum(1d-3*map_fha_err[igood_ha]*initdat.fluxunits,ldist,/ergs)
         
;      map_m_err_ha[igood_ha] = mumpsm * map_l_err_ha[igood_ha] / $
;                               volemis / elecden
      for j=0,1 do begin
         xyind = array_indices(map_r,igood_ha)
         xigdof = xyind[0,*]
         yigdof = xyind[1,*]
         for k=0,ctgood_ha-1 do begin
            map_m_err_ha[xigdof[k],yigdof[k],j] = $
               map_m_ha[xigdof[k],yigdof[k]]*$
               sqrt((map_l_err_ha[xigdof[k],yigdof[k]]/$
                     map_l_ha[xigdof[k],yigdof[k]])^2d + $
                    (elecdenuse_err[xigdof[k],yigdof[k],j]/$
                     elecdenuse[xigdof[k],yigdof[k]])^2d)
            map_dmdt_err_ha[xigdof[k],yigdof[k],j] = $
               map_dmdt_ha[xigdof[k],yigdof[k]] * $
               sqrt((map_m_err_ha[xigdof[k],yigdof[k],j]/map_m_ha[xigdof[k],yigdof[k]])^2d + $
                    (map_vrad_err_ha[xigdof[k],yigdof[k]]/map_vrad_ha[xigdof[k],yigdof[k]])^2d)
            map_p_err_ha[xigdof[k],yigdof[k],j] = $
               map_p_ha[xigdof[k],yigdof[k]] * $
               sqrt((map_m_err_ha[xigdof[k],yigdof[k],j]/map_m_ha[xigdof[k],yigdof[k]])^2d + $
                    (map_vrad_err_ha[xigdof[k],yigdof[k]]/map_vrad_ha[xigdof[k],yigdof[k]])^2d)
            map_dpdt_err_ha[xigdof[k],yigdof[k],j] = $
               map_dpdt_ha[xigdof[k],yigdof[k]] * $
               sqrt((map_m_err_ha[xigdof[k],yigdof[k],j]/map_m_ha[xigdof[k],yigdof[k]])^2d + $
                    2d*(map_vrad_err_ha[xigdof[k],yigdof[k]]/map_vrad_ha[xigdof[k],yigdof[k]])^2d)
            map_e_err_ha[xigdof[k],yigdof[k],j] = $
               sqrt($
                    (map_e_ha[xigdof[k],yigdof[k]]*$
                     map_m_err_ha[xigdof[k],yigdof[k],j]/$
                     map_m_ha[xigdof[k],yigdof[k]])^2d + $
                    (msun*map_m_ha[xigdof[k],yigdof[k]]*$
                     map_vrad_ha[xigdof[k],yigdof[k]]*$
                     map_vrad_err_ha[xigdof[k],yigdof[k]])^2d + $
                    (3d*msun*map_m_ha[xigdof[k],yigdof[k]]*$
                     map_sig_ha[xigdof[k],yigdof[k]]*$
                     map_sig_err_ha[xigdof[k],yigdof[k]])^2d)
            map_dedt_err_ha[xigdof[k],yigdof[k],j] = $
               sqrt($
                    (map_dedt_ha[xigdof[k],yigdof[k]]*$
                     map_dmdt_err_ha[xigdof[k],yigdof[k],j]/$
                     map_dmdt_ha[xigdof[k],yigdof[k]])^2d + $
                    (msun/speryr*map_dmdt_ha[xigdof[k],yigdof[k]]*$
                     map_vrad_ha[xigdof[k],yigdof[k]]*$
                     map_vrad_err_ha[xigdof[k],yigdof[k]])^2d + $
                    (3d*msun/speryr*map_dmdt_ha[xigdof[k],yigdof[k]]*$
                     map_sig_ha[xigdof[k],yigdof[k]]*$
                     map_sig_err_ha[xigdof[k],yigdof[k]])^2d)
          endfor
       endfor
;      fix low errors that are too big! probably should do this earlier, but OK
;      for now
       map_m_err_ha_lo = map_m_err_ha[*,*,0]
       ibdlo = where(map_m_err_ha_lo gt map_m_ha,ctbdlo)
       if ctbdlo gt 0 then begin
          map_m_err_ha_lo[ibdlo] = map_m_ha[ibdlo]
          map_m_err_ha[*,*,0] = map_m_err_ha_lo
       endif
       map_dmdt_err_ha_lo = map_dmdt_err_ha[*,*,0]
       ibdlo = where(map_dmdt_err_ha_lo gt map_dmdt_ha,ctbdlo)
       if ctbdlo gt 0 then begin
          map_dmdt_err_ha_lo[ibdlo] = map_dmdt_ha[ibdlo]
          map_dmdt_err_ha[*,*,0] = map_dmdt_err_ha_lo
       endif
       map_p_err_ha_lo = map_p_err_ha[*,*,0]
       ibdlo = where(map_p_err_ha_lo gt map_p_ha,ctbdlo)
       if ctbdlo gt 0 then begin
          map_p_err_ha_lo[ibdlo] = map_p_ha[ibdlo]
          map_p_err_ha[*,*,0] = map_p_err_ha_lo
       endif
       map_dpdt_err_ha_lo = map_dpdt_err_ha[*,*,0]
       ibdlo = where(map_dpdt_err_ha_lo gt map_dpdt_ha,ctbdlo)
       if ctbdlo gt 0 then begin
          map_dpdt_err_ha_lo[ibdlo] = map_dpdt_ha[ibdlo]
          map_dpdt_err_ha[*,*,0] = map_dpdt_err_ha_lo
       endif
       map_e_err_ha_lo = map_e_err_ha[*,*,0]
       ibdlo = where(map_e_err_ha_lo gt map_e_ha,ctbdlo)
       if ctbdlo gt 0 then begin
          map_e_err_ha_lo[ibdlo] = map_e_ha[ibdlo]
          map_e_err_ha[*,*,0] = map_e_err_ha_lo
       endif
       map_dedt_err_ha_lo = map_dedt_err_ha[*,*,0]
       ibdlo = where(map_dedt_err_ha_lo gt map_dedt_ha,ctbdlo)
       if ctbdlo gt 0 then begin
          map_dedt_err_ha_lo[ibdlo] = map_dedt_ha[ibdlo]
          map_dedt_err_ha[*,*,0] = map_dedt_err_ha_lo
       endif

; redshifted
      if ctgoodr_ha gt 0 then begin
         mapr_vrad_ha[igoodr_ha] = abs(mapr_vha[igoodr_ha]) / $
                                      map_costheta_ha[igoodr_ha]
         mapr_l_ha[igoodr_ha] = $
            drt_linelum(1d-3*mapr_fha[igoodr_ha]*initdat.fluxunits,ldist,/ergs)
         mapr_m_ha[igoodr_ha] = mumpsm * mapr_l_ha[igoodr_ha] / $
                                volemis / elecdenuse[igoodr_ha]
         mapr_dmdt_ha[igoodr_ha] = mapr_m_ha[igoodr_ha] * $
                                   mapr_vrad_ha[igoodr_ha] * speryr / Rcm_ha
         mapr_p_ha[igoodr_ha] = mapr_m_ha[igoodr_ha] * msun * $
                                mapr_vrad_ha[igoodr_ha]
         mapr_dpdt_ha[igoodr_ha] = mapr_dmdt_ha[igoodr_ha] * msun / speryr * $
                                   mapr_vrad_ha[igoodr_ha] * c_cms / lsun
         mapr_e_ha[igoodr_ha] = mapr_m_ha[igoodr_ha] * msun * $
                                (0.5d*mapr_vrad_ha[igoodr_ha]^2d + $
                                 1.5d*mapr_sig_ha[igoodr_ha]^2d)
         mapr_dedt_ha[igoodr_ha] = mapr_dmdt_ha[igoodr_ha] * msun / speryr * $
                                   (0.5d*mapr_vrad_ha[igoodr_ha]^2d + $
                                    1.5d*mapr_sig_ha[igoodr_ha]^2d)

         mapr_vrad_err_ha[igoodr_ha] = abs(mapr_vha_err[igoodr_ha]) / $
                                         map_costheta_ha[igoodr_ha]
         mapr_l_err_ha[igoodr_ha] = $
            drt_linelum(1d-3*mapr_fha_err[igoodr_ha]*initdat.fluxunits,ldist,/ergs)
 
         for j=0,1 do begin
         xyind = array_indices(map_r,igoodr_ha)
         xigdof = xyind[0,*]
         yigdof = xyind[1,*]
         for k=0,ctgoodr_ha-1 do begin
            mapr_m_err_ha[xigdof[k],yigdof[k],j] = $
               mapr_m_ha[xigdof[k],yigdof[k]]*$
               sqrt((mapr_l_err_ha[xigdof[k],yigdof[k]]/$
                     mapr_l_ha[xigdof[k],yigdof[k]])^2d + $
                    (elecdenuse_err[xigdof[k],yigdof[k],j]/$
                     elecdenuse[xigdof[k],yigdof[k]])^2d)
            mapr_dmdt_err_ha[xigdof[k],yigdof[k],j] = $
               mapr_dmdt_ha[xigdof[k],yigdof[k]] * $
               sqrt((mapr_m_err_ha[xigdof[k],yigdof[k],j]/mapr_m_ha[xigdof[k],yigdof[k]])^2d + $
                    (mapr_vrad_err_ha[xigdof[k],yigdof[k]]/mapr_vrad_ha[xigdof[k],yigdof[k]])^2d)
            mapr_p_err_ha[xigdof[k],yigdof[k],j] = $
               mapr_p_ha[xigdof[k],yigdof[k]] * $
               sqrt((mapr_m_err_ha[xigdof[k],yigdof[k],j]/mapr_m_ha[xigdof[k],yigdof[k]])^2d + $
                    (mapr_vrad_err_ha[xigdof[k],yigdof[k]]/mapr_vrad_ha[xigdof[k],yigdof[k]])^2d)
            mapr_dpdt_err_ha[xigdof[k],yigdof[k],j] = $
               mapr_dpdt_ha[xigdof[k],yigdof[k]] * $
               sqrt((mapr_m_err_ha[xigdof[k],yigdof[k],j]/mapr_m_ha[xigdof[k],yigdof[k]])^2d + $
                    2d*(mapr_vrad_err_ha[xigdof[k],yigdof[k]]/mapr_vrad_ha[xigdof[k],yigdof[k]])^2d)
            mapr_e_err_ha[xigdof[k],yigdof[k],j] = $
               sqrt($
                    (mapr_e_ha[xigdof[k],yigdof[k]]*$
                     mapr_m_err_ha[xigdof[k],yigdof[k],j]/$
                     mapr_m_ha[xigdof[k],yigdof[k]])^2d + $
                    (msun*mapr_m_ha[xigdof[k],yigdof[k]]*$
                     mapr_vrad_ha[xigdof[k],yigdof[k]]*$
                     mapr_vrad_err_ha[xigdof[k],yigdof[k]])^2d + $
                    (3d*msun*mapr_m_ha[xigdof[k],yigdof[k]]*$
                     mapr_sig_ha[xigdof[k],yigdof[k]]*$
                     mapr_sig_err_ha[xigdof[k],yigdof[k]])^2d)
            mapr_dedt_err_ha[xigdof[k],yigdof[k],j] = $
               sqrt($
                    (mapr_dedt_ha[xigdof[k],yigdof[k]]*$
                     mapr_dmdt_err_ha[xigdof[k],yigdof[k],j]/$
                     mapr_dmdt_ha[xigdof[k],yigdof[k]])^2d + $
                    (msun/speryr*mapr_dmdt_ha[xigdof[k],yigdof[k]]*$
                     mapr_vrad_ha[xigdof[k],yigdof[k]]*$
                     mapr_vrad_err_ha[xigdof[k],yigdof[k]])^2d + $
                    (3d*msun/speryr*mapr_dmdt_ha[xigdof[k],yigdof[k]]*$
                     mapr_sig_ha[xigdof[k],yigdof[k]]*$
                     mapr_sig_err_ha[xigdof[k],yigdof[k]])^2d)
             endfor
          endfor
;      fix low errors that are too big! probably should do this earlier, but OK
;      for now
       mapr_m_err_ha_lo = mapr_m_err_ha[*,*,0]
       ibdlo = where(mapr_m_err_ha_lo gt mapr_m_ha,ctbdlo)
       if ctbdlo gt 0 then begin
          mapr_m_err_ha_lo[ibdlo] = mapr_m_ha[ibdlo]
          mapr_m_err_ha[*,*,0] = mapr_m_err_ha_lo
       endif
       mapr_dmdt_err_ha_lo = mapr_dmdt_err_ha[*,*,0]
       ibdlo = where(mapr_dmdt_err_ha_lo gt mapr_dmdt_ha,ctbdlo)
       if ctbdlo gt 0 then begin
          mapr_dmdt_err_ha_lo[ibdlo] = mapr_dmdt_ha[ibdlo]
          mapr_dmdt_err_ha[*,*,0] = mapr_dmdt_err_ha_lo
       endif
       mapr_p_err_ha_lo = mapr_p_err_ha[*,*,0]
       ibdlo = where(mapr_p_err_ha_lo gt mapr_p_ha,ctbdlo)
       if ctbdlo gt 0 then begin
          mapr_p_err_ha_lo[ibdlo] = mapr_p_ha[ibdlo]
          mapr_p_err_ha[*,*,0] = mapr_p_err_ha_lo
       endif
       mapr_dpdt_err_ha_lo = mapr_dpdt_err_ha[*,*,0]
       ibdlo = where(mapr_dpdt_err_ha_lo gt mapr_dpdt_ha,ctbdlo)
       if ctbdlo gt 0 then begin
          mapr_dpdt_err_ha_lo[ibdlo] = mapr_dpdt_ha[ibdlo]
          mapr_dpdt_err_ha[*,*,0] = mapr_dpdt_err_ha_lo
       endif
       mapr_e_err_ha_lo = mapr_e_err_ha[*,*,0]
       ibdlo = where(mapr_e_err_ha_lo gt mapr_e_ha,ctbdlo)
       if ctbdlo gt 0 then begin
          mapr_e_err_ha_lo[ibdlo] = mapr_e_ha[ibdlo]
          mapr_e_err_ha[*,*,0] = mapr_e_err_ha_lo
       endif
       mapr_dedt_err_ha_lo = mapr_dedt_err_ha[*,*,0]
       ibdlo = where(mapr_dedt_err_ha_lo gt mapr_dedt_ha,ctbdlo)
       if ctbdlo gt 0 then begin
          mapr_dedt_err_ha_lo[ibdlo] = mapr_dedt_ha[ibdlo]
          mapr_dedt_err_ha[*,*,0] = mapr_dedt_err_ha_lo
       endif
       endif


      if ctgoodr_ha gt 0 then begin
;         m_err_ha = m_ha * median([map_m_err_ha[igood_ha]/map_m_ha[igood_ha],$
;                                   mapr_m_err_ha[igoodr_ha]/mapr_m_ha[igoodr_ha]])
;         dmdt_err_ha = dmdt_ha * median([map_dmdt_err_ha[igood_ha]/map_dmdt_ha[igood_ha],$
;                                         mapr_dmdt_err_ha[igoodr_ha]/mapr_dmdt_ha[igoodr_ha]])
;         p_err_ha = p_ha * median([map_p_err_ha[igood_ha]/map_p_ha[igood_ha],$
;                                   mapr_p_err_ha[igoodr_ha]/mapr_p_ha[igoodr_ha]])
;         dpdt_err_ha = dpdt_ha * median([map_dpdt_err_ha[igood_ha]/map_dpdt_ha[igood_ha],$
;                                         mapr_dpdt_err_ha[igoodr_ha]/mapr_dpdt_ha[igoodr_ha]])
         m_ha = total(map_m_ha)+total(mapr_m_ha)
         dmdt_ha = total(map_dmdt_ha)+total(mapr_dmdt_ha)
         p_ha = total(map_p_ha)+total(mapr_p_ha)
         dpdt_ha = total(map_dpdt_ha)+total(mapr_dpdt_ha)
         e_ha = total(map_e_ha)+total(mapr_e_ha)
         dedt_ha = total(map_dedt_ha)+total(mapr_dedt_ha)
         omega_ha = 4/!DPi*(total(map_domega_ha[igood_ha])+$
                         total(map_domega_ha[igoodr_ha]))
                 
         m_err_ha = sqrt(total(total(map_m_err_ha^2d,1),1)+$
                            total(total(mapr_m_err_ha^2d,1),1))
         dmdt_err_ha = sqrt(total(total(map_dmdt_err_ha^2d,1),1)+$
                               total(total(mapr_dmdt_err_ha^2d,1),1))
         p_err_ha = sqrt(total(total(map_p_err_ha^2d,1),1)+$
                            total(total(mapr_p_err_ha^2d,1),1))
         dpdt_err_ha = sqrt(total(total(map_dpdt_err_ha^2d,1),1)+$
                               total(total(mapr_dpdt_err_ha^2d,1),1))
         e_err_ha = sqrt(total(total(map_e_err_ha^2d,1),1)+$
                            total(total(mapr_e_err_ha^2d,1),1))
         dedt_err_ha = sqrt(total(total(map_dedt_err_ha^2d,1),1)+$
                               total(total(mapr_dedt_err_ha^2d,1),1))
      endif else begin
;         m_err_ha = m_ha * median(map_m_err_ha/map_m_ha)
;         dmdt_err_ha = dmdt_ha * median(map_dmdt_err_ha/map_dmdt_ha)
;         p_err_ha = p_ha * median(map_p_err_ha/map_p_ha)
;         dpdt_err_ha = dpdt_ha * median(map_dpdt_err_ha/map_dpdt_ha)
         m_ha = 2d*total(map_m_ha)
         dmdt_ha = 2d*total(map_dmdt_ha)
         p_ha = 2d*total(map_p_ha)
         dpdt_ha = 2d*total(map_dpdt_ha)
         e_ha = 2d*total(map_e_ha)
         dedt_ha = 2d*total(map_dedt_ha)
         omega_ha = 2d/4/!DPi*total(map_domega_ha[igood_ha])
                 
         m_err_ha = 2d*sqrt(total(total(map_m_err_ha^2d,1),1))
         dmdt_err_ha = 2d*sqrt(total(total(map_dmdt_err_ha^2d,1),1))
         p_err_ha = 2d*sqrt(total(total(map_p_err_ha^2d,1),1))
         dpdt_err_ha = 2d*sqrt(total(total(map_dpdt_err_ha^2d,1),1))
         e_err_ha = 2d*sqrt(total(total(map_e_err_ha^2d,1),1))
         dedt_err_ha = 2d*sqrt(total(total(map_dedt_err_ha^2d,1),1))

      endelse
      
   endif
   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Compute fluxes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  e_of_flx = 0
  if ~ tag_exist(initdat,'noemlinfit') then begin

  fluxnorm = initdat.fluxunits ; assume fluxunits incorporates surface brightness correction if necessary
;  fluxnorm = initdat.fluxunits * initdat.platescale^2d
  emlflxsums = hash()
  emlflxsums['tot_ext'] = hash()
  emlflxsums['tot_unext_pp'] = hash()
  emlflxsums['tot_unext_med'] = hash()
  emlflxsums['comp_ext'] = hash()
  emlflxsums['comp_unext_pp'] = hash()
  emlflxsums['comp_unext_med'] = hash()
  emlflxsums['of_ext'] = hash()
  emlflxsums['of_unext_pp'] = hash()
  emlflxsums['of_unext_med'] = hash()
  foreach line,linelist.keys() do begin
     emlflxsums['tot_ext',line] = 0d
     emlflxsums['tot_unext_pp',line] = 0d
     emlflxsums['tot_unext_med',line] = 0d
     emlflxsums['comp_ext',line] = dblarr(initdat.maxncomp)
     emlflxsums['comp_unext_pp',line] = dblarr(initdat.maxncomp)
     emlflxsums['comp_unext_med',line] = dblarr(initdat.maxncomp)
     emlflxsums['of_ext',line] = 0d
     emlflxsums['of_unext_pp',line] = 0d
     emlflxsums['of_unext_med',line] = 0d
     for icomp=0,initdat.maxncomp-1 do begin
        stric = string(icomp+1,format='(I0)')
        igdflx = where(emlflx['fc'+stric,line] ne bad,ctgdflx)
        if line eq ofcompline AND $
           icomp eq 1 AND $
           tag_exist(initmaps,'outflow_eml') then begin
           if tag_exist(initmaps.outflow_eml,'usered') then begin
              igdblue = where(emlvel['v%50c'+stric,ofcompline] lt 0d AND $
                              map_rkpc_ifs lt initmaps.outflow_eml.R,ctgdblue)
              igdred = where(emlvel['v%50c'+stric,ofcompline] ne bad AND $
                             emlvel['v%50c'+stric,ofcompline] ge 0 AND $
                             map_rkpc_ifs lt initmaps.outflow_eml.R,ctgdred)
              if ctgdblue gt 0 then igdflx_of = [igdblue,igdred] $
              else igdflx_of = igdred
           endif else begin
              igdflx_of = where(emlvel['v%50c'+stric,ofcompline] lt 0d AND $
                                map_rkpc_ifs lt initmaps.outflow_eml.R,ctgdflx_of)
           endelse
        endif
        if ctgdflx gt 0 then begin
           emlflxsums['comp_ext',line,icomp] = $
              total(emlflx['fc'+stric,line,igdflx]) * fluxnorm
           if line eq ofcompline AND $
              icomp eq 1 AND $
              tag_exist(initmaps,'outflow_eml') then $
                 emlflxsums['of_ext',line] = $
                    total(emlflx['fc'+stric,line,igdflx_of]) * fluxnorm
           emlflxsums['tot_ext',line] += $
              emlflxsums['comp_ext',line,icomp]
           if tag_exist(initmaps,'ebv') then begin
              if tag_exist(initmaps.ebv,'apply') then begin
                 emlflxsums['comp_unext_pp',line,icomp] = $
                    total(emlflxcor_pp['fc'+stric,line,igdflx]) * fluxnorm
                 emlflxsums['comp_unext_med',line,icomp] = $
                    total(emlflxcor_med['fc'+stric,line,igdflx]) * fluxnorm
                 if line eq ofcompline AND $
                    icomp eq 1 AND $
                    tag_exist(initmaps,'outflow_eml') then begin
                    emlflxsums['of_unext_pp',line] = $
                       total(emlflxcor_pp['fc'+stric,line,igdflx_of]) * fluxnorm
                    emlflxsums['of_unext_med',line] = $
                       total(emlflxcor_med['fc'+stric,line,igdflx_of]) * fluxnorm
                 endif
                 emlflxsums['tot_unext_pp',line] += $
                    emlflxsums['comp_unext_pp',line,icomp]
                 emlflxsums['tot_unext_med',line] += $
                    emlflxsums['comp_unext_med',line,icomp]
              endif
           endif
        endif
     endfor
  endforeach

; Outflow flux
  if tag_exist(initmaps,'outflow_eml') then begin
     usekey = 'of_ext'
     if tag_exist(initmaps,'ebv') then $
        if tag_exist(initmaps.ebv,'apply') then $
           usekey = 'of_unext_pp'
     e_of_flx = emlflxsums[usekey,ofcompline]
     if ofcompline eq 'Hbeta' then e_of_flx *= 2.86d
  endif

  endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; TABLE AND STRUCTURE OUTPUTS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   openw,lun_stats,initdat.mapdir+initdat.label+'.stat.txt',/get_lun

;  Peak fitting
  
   if tag_exist(initdat,'decompose_qso_fit') then begin
      ; Print output of Moffat fit to QSO PSF
      printf,lun_stats,'2D Moffat fit to QSO PSF:' 
      printf,lun_stats,'  peak = ',string(qso_fitpar[1],format='(D0.4)'),$
         ' flux units'
      printf,lun_stats,'  cent (single-offset indices) = ',$
         string('[',qso_fitpar[4]+1,',',qso_fitpar[5]+1,']',$
         format='(A0,D0.2,A0,D0.2,A0)')
      printf,lun_stats,'  FWHM = ',$
         string(qso_fitpar[2],format='(D0.4)'),$
         ' spaxels'
      printf,lun_stats,'  FWHM = ',$
         string(qso_fitpar[2]*initdat.platescale,format='(D0.4)'),$
         ' arcseconds'
       printf,lun_stats,'  FWHM = ',$
         string(qso_fitpar[2]*initdat.platescale*kpc_per_as,format='(D0.4)'),$
         ' kpc'
      printf,lun_stats,'  index = ',$
         string(qso_fitpar[7],format='(D0.4)')
   endif
   if tag_exist(initmaps,'fit_empsf') then begin
      ; Print output of Moffat fit to emission line PSF
      printf,lun_stats,'2D Moffat fit to Emission Line PSF:' 
      printf,lun_stats,'  peak = ',string(empsf_fitpar[1],format='(D0.4)'),$
         ' flux units'
      printf,lun_stats,'  cent (single-offset indices) = ',$
         string('[',empsf_fitpar[4]+1,',',empsf_fitpar[5]+1,']',$
         format='(A0,D0.2,A0,D0.2,A0)')
      printf,lun_stats,'  FWHM = ',$
         string(empsf_fitpar[2],format='(D0.4)'),$
         ' spaxels'
      printf,lun_stats,'  FWHM = ',$
         string(empsf_fitpar[2]*initdat.platescale,format='(D0.4)'),$
         ' arcseconds'
      printf,lun_stats,'  FWHM = ',$
         string(empsf_fitpar[2]*initdat.platescale*kpc_per_as,format='(D0.4)'),$
         ' kpc'
      printf,lun_stats,'  index = ',$
         string(empsf_fitpar[7],format='(D0.4)')
   endif
   did_distance_intro = 0b
   if dohst then begin
      if tag_exist(initmaps.hst,'fithstpeak') then begin
         printf,lun_stats,lineofdashes
         printf,lun_stats,'Distances between assumed peak and actual (fitted) peak'
         printf,lun_stats,lineofdashes
         did_distance_intro = 1b
         printf,lun_stats,'  HST image (rotated, trimmed to IFS FOV)'
         printf,lun_stats,'    In HST pixels: ',$
                string(peakfit_hst_distance_from_nucleus_hstpix[0],'  ',$
                       peakfit_hst_distance_from_nucleus_hstpix[1],$
                       format='(D0.3,A0,D0.3)')
         printf,lun_stats,'    In kpc: ',$
                string(peakfit_hst_distance_from_nucleus_kpc[0],'  ',$
                       peakfit_hst_distance_from_nucleus_kpc[1],$
                       format='(D0.3,A0,D0.3)')
      endif
      if tag_exist(initmaps.ct,'fitifspeak') AND dohstsm then begin
         if not did_distance_intro then begin
            printf,lun_stats,lineofdashes
            printf,lun_stats,'Distances between assumed peak and actual (fitted) peak'
            printf,lun_stats,lineofdashes
         endif
         printf,lun_stats,'  HST image (rotated, trimmed, convolved, rebinned)'
         printf,lun_stats,'    In pixels: ',$
                string(peakfit_hstconv_distance_from_nucleus_pix[0],'  ',$
                       peakfit_hstconv_distance_from_nucleus_pix[1],$
                       format='(D0.3,A0,D0.3)')
         printf,lun_stats,'    In kpc: ',$
                string(peakfit_hstconv_distance_from_nucleus_kpc[0],'  ',$
                       peakfit_hstconv_distance_from_nucleus_kpc[1],$
                       format='(D0.3,A0,D0.3)')
            
      endif
   endif
   if tag_exist(initmaps,'ct') then begin
      if tag_exist(initmaps.ct,'fitifspeak') then begin
         if not did_distance_intro then begin
            printf,lun_stats,lineofdashes
            printf,lun_stats,'Distances between assumed peak and actual (fitted) peak'
            printf,lun_stats,lineofdashes
         endif
         printf,lun_stats,'  IFS continuum'
         printf,lun_stats,'    In pixels: ',$
                string(peakfit_ifs_distance_from_nucleus_pix[0],'  ',$
                       peakfit_ifs_distance_from_nucleus_pix[1],$
                       format='(D0.3,A0,D0.3)')
         printf,lun_stats,'    In kpc: ',$
                string(peakfit_ifs_distance_from_nucleus_kpc[0],'  ',$
                       peakfit_ifs_distance_from_nucleus_kpc[1],$
                       format='(D0.3,A0,D0.3)')
         if tag_exist(initdat,'decompose_qso_fit') then begin
            printf,lun_stats,'  IFS continuum, QSO only'
            printf,lun_stats,'    In pixels: ',$
                   string(peakfit_ifs_qso_distance_from_nucleus_pix[0],'  ',$
                          peakfit_ifs_qso_distance_from_nucleus_pix[1],$
                          format='(D0.3,A0,D0.3)')
            printf,lun_stats,'    In kpc: ',$
                   string(peakfit_ifs_qso_distance_from_nucleus_kpc[0],'  ',$
                          peakfit_ifs_qso_distance_from_nucleus_kpc[1],$
                          format='(D0.3,A0,D0.3)')
         endif
            
      endif
   endif


   a_vel_stats = dblarr(6)+bad
   a_fwhm_stats = dblarr(6)+bad
   a_v98_stats = dblarr(6)+bad
   e_vel_stats = dblarr(6)+bad
   e_fwhm_stats = dblarr(6)+bad
   e_v98_stats = dblarr(6)+bad
   a_vel_stats_all = dblarr(6)+bad
   a_fwhm_stats_all = dblarr(6)+bad
   a_v98_stats_all = dblarr(6)+bad
   e_vel_stats_all = dblarr(6)+bad
   e_fwhm_stats_all = dblarr(6)+bad
   e_v98_stats_all = dblarr(6)+bad

   map_of_abs = bytarr(dx,dy)
   a_of_meanxy = dblarr(2)+bad
   a_of_meanr = bad
   a_of_meanpa = bad
   a_of_meanxy_wtv50 = dblarr(2)+bad
   a_of_meanr_wtv50 = bad
   a_of_meanpa_wtv50 = bad
   a_of_meanxy_wtv98 = dblarr(2)+bad
   a_of_meanr_wtv98 = bad
   a_of_meanpa_wtv98 = bad
   a_of_meanxy_wtweq = dblarr(2)+bad
   a_of_meanr_wtweq = bad
   a_of_meanpa_wtweq = bad
   map_of_eml = bytarr(dx,dy)
   e_of_meanxy = dblarr(2)+bad
   e_of_meanr = bad
   e_of_meanpa = bad
   e_of_meanxy_wtv50 = dblarr(2)+bad
   e_of_meanr_wtv50 = bad
   e_of_meanpa_wtv50 = bad
   e_of_meanxy_wtv98 = dblarr(2)+bad
   e_of_meanr_wtv98 = bad
   e_of_meanpa_wtv98 = bad
   e_of_meanxy_wtflx = dblarr(2)+bad
   e_of_meanr_wtflx = bad
   e_of_meanpa_wtflx = bad

   statformat = '(A-5,A12,5D8.0,I5)'
   if tag_exist(initdat,'donad') then begin

      printf,lun_stats,lineofdashes
      printf,lun_stats,'ABSORPTION'
      printf,lun_stats,lineofdashes
      printf,lun_stats,'All Velocities'
      printf,lun_stats,lineofdashes
      printf,lun_stats,'#Cmp','Quantity','Mean','Median','Min','Max','StdDev','#',$
             format='(A-5,A12,5A8,A5)'
      printf,lun_stats,lineofdashes
      icomp = 0
      arr = nadabscvdfvals['v%50']
      igd = where(arr ne bad,ctgd)
      arr = arr[igd]
      printf,lun_stats,string(icomp,format='(I0)'),'Vel',mean(arr),median(arr),$
             min(arr),max(arr),stddev(arr),ctgd,format=statformat
      a_vel_stats_all = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]
      arr = nadabscvdfvals['vsig']
      igd = where(arr ne bad)
      arr = arr[igd] * 2d * sqrt(2d*alog(2d))
      printf,lun_stats,string(icomp,format='(I0)'),'FWHM',mean(arr),median(arr),$
             min(arr),max(arr),stddev(arr),ctgd,format=statformat
      a_fwhm_stats_all = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]
      arr = nadabscvdfvals['v%98']
      igd = where(arr ne bad)
      arr = arr[igd]
      printf,lun_stats,string(icomp,format='(I0)'),'V98',mean(arr),median(arr),$
             min(arr),max(arr),stddev(arr),ctgd,format=statformat
      a_v98_stats_all = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]

      if tag_exist(initmaps,'outflow_abs') then begin

         arr = nadabscvdfvals['v%50']
         igd = where(nadabscvdfvals['v%50'] lt 0 AND $
                     map_rkpc_ifs lt initmaps.outflow_abs.R,ctgd)

         printf,lun_stats,lineofdashes
         printf,lun_stats,'Blueshifted Only'
         printf,lun_stats,lineofdashes
         printf,lun_stats,'#Cmp','Quantity','Mean','Median','Min','Max','StdDev','#',$
                format='(A-5,A12,5A8,A5)'
         printf,lun_stats,lineofdashes
         icomp = 0
         arr = nadabscvdfvals['v%50']
         igd = where(nadabscvdfvals['v%50'] lt 0 AND $
                     map_rkpc_ifs lt initmaps.outflow_abs.R,ctgd)
         arr = arr[igd]
         printf,lun_stats,string(icomp,format='(I0)'),'Vel',mean(arr),median(arr),$
                min(arr),max(arr),stddev(arr),ctgd,format=statformat
         a_vel_stats = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]
         arr = nadabscvdfvals['vsig']
         arr = arr[igd] * 2d * sqrt(2d*alog(2d))
         printf,lun_stats,string(icomp,format='(I0)'),'FWHM',mean(arr),median(arr),$
                min(arr),max(arr),stddev(arr),ctgd,format=statformat
         a_fwhm_stats = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]
         arr = nadabscvdfvals['v%98']
         arr = arr[igd]
         printf,lun_stats,string(icomp,format='(I0)'),'V98',mean(arr),median(arr),$
                min(arr),max(arr),stddev(arr),ctgd,format=statformat
         a_v98_stats = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]
      
;        map of outflow spaxels
         map_of_abs = bytarr(dx,dy)
         map_of_abs[igd] = 1b

;        Outflow coordinates
         a_of_meanxy = [mean(map_x[igd])-center_axes[0],$
                        mean(map_y[igd])-center_axes[1]]*kpc_per_pix
         a_of_meanr = sqrt((mean(map_x[igd])-center_axes[0])^2d + $
                           (mean(map_y[igd])-center_axes[1])^2d)*kpc_per_pix
         a_of_meanpa = ifsf_pa(mean(map_x[igd])-center_axes[0],$
                               mean(map_y[igd])-center_axes[1])
         printf,lun_stats,icomp,'        PAof',a_of_meanpa,format='(I-0,A12,I5)'

         wtv50 = nadabscvdfvals['v%50',igd]^2d/total(nadabscvdfvals['v%50',igd]^2d)
         a_of_meanxy_wtv50 = [mean((map_x[igd]-center_axes[0])*wtv50),$
                              mean((map_y[igd]-center_axes[1])*wtv50)]*kpc_per_pix*ctgd
         a_of_meanr_wtv50 = sqrt((mean((map_x[igd]-center_axes[0])*wtv50))^2d + $
                                 (mean((map_y[igd]-center_axes[1])*wtv50))^2d)*kpc_per_pix*ctgd
         a_of_meanpa_wtv50 = ifsf_pa(mean((map_x[igd]-center_axes[0])*wtv50),$
                                     mean((map_y[igd]-center_axes[1])*wtv50))
         printf,lun_stats,icomp,'       PAv50',a_of_meanpa_wtv50,format='(I-0,A12,I5)'

         wtv98 = nadabscvdfvals['v%98',igd]^2d/total(nadabscvdfvals['v%98',igd]^2d)
         a_of_meanxy_wtv98 = [mean((map_x[igd]-center_axes[0])*wtv98),$
                              mean((map_y[igd]-center_axes[1])*wtv98)]*kpc_per_pix*ctgd
         a_of_meanr_wtv98 = sqrt((mean((map_x[igd]-center_axes[0])*wtv98))^2d + $
                                 (mean((map_y[igd]-center_axes[1])*wtv98))^2d)*kpc_per_pix*ctgd
         a_of_meanpa_wtv98 = ifsf_pa(mean((map_x[igd]-center_axes[0])*wtv98),$
                                     mean((map_y[igd]-center_axes[1])*wtv98))
         printf,lun_stats,icomp,'       PAv98',a_of_meanpa_wtv98,format='(I-0,A12,I5)'

         wtweq = map_nadabs_fitweq[igd]/total(map_nadabs_fitweq[igd])
         a_of_meanxy_wtweq = [mean((map_x[igd]-center_axes[0])*wtweq),$
                              mean((map_y[igd]-center_axes[1])*wtweq)]*kpc_per_pix*ctgd
         a_of_meanr_wtweq = sqrt((mean((map_x[igd]-center_axes[0])*wtweq))^2d + $
                                 (mean((map_y[igd]-center_axes[1])*wtweq))^2d)*kpc_per_pix*ctgd
         a_of_meanpa_wtweq = ifsf_pa(mean((map_x[igd]-center_axes[0])*wtweq),$
                                     mean((map_y[igd]-center_axes[1])*wtweq))
         printf,lun_stats,icomp,'       PAweq',a_of_meanpa_wtweq,format='(I-0,A12,I5)'

         printf,lun_stats,lineofdashes
         printf,lun_stats,'CRFW model'
         printf,lun_stats,lineofdashes
         printf,lun_stats,$
                'Rmax(obs)=',string(Rmax_kpc_nad,format='(D0.2)'),' kpc'         
         printf,lun_stats,$
                'Rmax(use)=',string(initmaps.outflow_abs.R,format='(D0.2)'),' kpc'
         printf,lun_stats,$
                '     M = ',string(m_nad,format='(E0.2)'),' M_sun'
         printf,lun_stats,$
                ' dM/dt = ',string(dmdt_nad,format='(E0.2)'),' M_sun/yr'
         printf,lun_stats,$
                '     p = ',string(p_nad,format='(E0.2)'),' dyne s'
         printf,lun_stats,$
                'dp/dt*c= ',string(dpdt_nad,format='(E0.2)'),' Lsun'
         printf,lun_stats,$
                '     E = ',string(e_nad,format='(E0.2)'),' erg'
         printf,lun_stats,$
                ' dE/dt = ',string(dedt_nad,format='(E0.2)'),' erg/s'

      endif

   endif
   if ~ tag_exist(initdat,'noemlinfit') then begin

      printf,lun_stats,lineofdashes
      printf,lun_stats,'EMISSION'
      printf,lun_stats,lineofdashes
      printf,lun_stats,'#Cmp','Quantity','Mean','Median','Min','Max','StdDev','#',$
        format='(A-5,A12,5A8,A5)'

      if ~ elecdenmap.isempty() then begin
         printf,lun_stats,lineofdashes
         printf,lun_stats,'Total'
         printf,lun_stats,lineofdashes
         arr = elecdenmap['ftot']
         igd = where(arr ne bad,ctgd)
         arr = arr[igds2]
         linstddev = stddev(arr)
         printf,lun_stats,stric,'n_e',10d^mean(arr),10d^median(arr),$
                10d^min(arr),10d^max(arr),$
                ((10d^(mean(arr)+stddev(arr))-10d^mean(arr)) + $
                (10d^mean(arr)-10d^(mean(arr)-stddev(arr))))/1d,$
                ctgd,format=statformat
      endif

      for icomp=1,initdat.maxncomp do begin
         stric = string(icomp,format='(I0)')
         printf,lun_stats,lineofdashes
         printf,lun_stats,'Component ',stric
         printf,lun_stats,lineofdashes
         arr = emlvel['v%50c'+stric,ofcompline]
         igd = where(arr ne bad,ctgd)
         arr = arr[igd]
         printf,lun_stats,stric,'v%50',mean(arr),median(arr),$
                min(arr),max(arr),stddev(arr),ctgd,format=statformat
         if tag_exist(initmaps,'outflow_eml') AND icomp eq 2 then $
            e_vel_stats_all = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]
         arr = emlvel['vsigc'+stric,ofcompline]
         igd = where(arr ne bad,ctgd)
         arr = arr[igd] * 2d * sqrt(2d*alog(2d))
         printf,lun_stats,string(icomp,format='(I0)'),'FWHM',mean(arr),median(arr),$
                min(arr),max(arr),stddev(arr),ctgd,format=statformat
         if tag_exist(initmaps,'outflow_eml') AND icomp eq 2 then $
            e_fwhm_stats_all = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]
         arr = emlvel['v%98c'+stric,ofcompline]
         igd = where(arr ne bad,ctgd)
         arr = arr[igd]
         printf,lun_stats,string(icomp,format='(I0)'),'V98',mean(arr),median(arr),$
                min(arr),max(arr),stddev(arr),ctgd,format=statformat
         if tag_exist(initmaps,'outflow_eml') AND icomp eq 2 then $
            e_v98_stats_all = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]

         if tag_exist(initmaps,'outflow_eml') AND icomp eq 2 then begin
          
            stric = string(icomp,format='(I0)')
            printf,lun_stats,lineofdashes
            printf,lun_stats,'Component ',stric,', Outflow Only'
            printf,lun_stats,lineofdashes
            arrv50 = emlvel['v%50c'+stric,ofcompline]
            arrv98 = emlvel['v%98c'+stric,ofcompline]
            if tag_exist(initmaps.outflow_eml,'usered') then begin
               igdblue = where(arrv50 lt 0 AND $
                               map_rkpc_ifs lt initmaps.outflow_eml.R,ctgdblue)
               igdred = where(arrv50 ne bad AND arrv50 ge 0 AND $
                              map_rkpc_ifs lt initmaps.outflow_eml.R,ctgdred)
               if ctgdred gt 0 then begin
                  arrv50[igdred] *= -1d
                  arrv98[igdred] *= -1d
               endif
               if ctgdblue gt 0 then igd = [igdblue,igdred] else igd = igdred
            endif else begin
               igd = where(arrv50 lt 0 AND $
                           map_rkpc_ifs lt initmaps.outflow_eml.R,ctgd)
            endelse
            arr = arrv50
            arr = arr[igd]
            printf,lun_stats,stric,'v%50',mean(arr),median(arr),$
                   min(arr),max(arr),stddev(arr),ctgd,format=statformat
            e_vel_stats = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]
            arr = emlvel['vsigc'+stric,ofcompline]
            arr = arr[igd] * 2d * sqrt(2d*alog(2d))
            printf,lun_stats,string(icomp,format='(I0)'),'FWHM',mean(arr),median(arr),$
                   min(arr),max(arr),stddev(arr),ctgd,format=statformat
            e_fwhm_stats = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]
            arr = arrv98            
            arr = arr[igd]
            printf,lun_stats,string(icomp,format='(I0)'),'V98',mean(arr),median(arr),$
                   min(arr),max(arr),stddev(arr),ctgd,format=statformat
            e_v98_stats = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]

            if ~ elecdenmap.isempty() then begin
               arr = elecdenmap['fc'+stric]
               igds2 = where(arr ne bad AND finite(arr),ctgds2)
               igds2 = cgsetintersection(igds2,igd)
               arr = arr[igds2]
               printf,lun_stats,stric,'n_e',10d^mean(arr),10d^median(arr),$
                  10d^min(arr),10d^max(arr),10d^stddev(arr),ctgd,format=statformat
            endif

;           Show blueshifted outflow only
            if tag_exist(initmaps.outflow_eml,'usered') then begin
               igd = igdblue
               ctgd = ctgdblue
            endif
            
;           Outflow coordinates
            e_of_meanxy = [mean(map_x[igd])-center_axes[0],$
                           mean(map_y[igd])-center_axes[1]]*kpc_per_pix
            e_of_meanr = sqrt((mean(map_x[igd])-center_axes[0])^2d + $
                              (mean(map_y[igd])-center_axes[1])^2d)*kpc_per_pix
            e_of_meanpa = ifsf_pa(mean(map_x[igd])-center_axes[0],$
                                  mean(map_y[igd])-center_axes[1])
            printf,lun_stats,icomp,'        PAof',e_of_meanpa,format='(I-0,A12,I5)'

            wtv50 = emlvel['v%50c'+stric,ofcompline,igd]^2d/$
                    total(emlvel['v%50c'+stric,ofcompline,igd]^2d)
            e_of_meanxy_wtv50 = [mean((map_x[igd]-center_axes[0])*wtv50),$
                                 mean((map_y[igd]-center_axes[1])*wtv50)]*kpc_per_pix*ctgd
            e_of_meanr_wtv50 = sqrt((mean((map_x[igd]-center_axes[0])*wtv50))^2d + $
                                    (mean((map_y[igd]-center_axes[1])*wtv50))^2d)*kpc_per_pix*ctgd
            e_of_meanpa_wtv50 = ifsf_pa(mean((map_x[igd]-center_axes[0])*wtv50),$
                                        mean((map_y[igd]-center_axes[1])*wtv50))
            printf,lun_stats,icomp,'       PAv50',e_of_meanpa_wtv50,format='(I-0,A12,I5)'

            wtv98 = emlvel['v%98c'+stric,ofcompline,igd]^2d/$
                    total(emlvel['v%98c'+stric,ofcompline,igd]^2d)
            e_of_meanxy_wtv98 = [mean((map_x[igd]-center_axes[0])*wtv98),$
                                 mean((map_y[igd]-center_axes[1])*wtv98)]*kpc_per_pix*ctgd
            e_of_meanr_wtv98 = sqrt((mean((map_x[igd]-center_axes[0])*wtv98))^2d + $
                                    (mean((map_y[igd]-center_axes[1])*wtv98))^2d)*kpc_per_pix*ctgd
            e_of_meanpa_wtv98 = ifsf_pa(mean((map_x[igd]-center_axes[0])*wtv98),$
                                        mean((map_y[igd]-center_axes[1])*wtv98))
            printf,lun_stats,icomp,'       PAv98',e_of_meanpa_wtv98,format='(I-0,A12,I5)'
            
            wtflx = emlflx['fc'+stric,ofcompline,igd]/$
                    total(emlflx['fc'+stric,ofcompline,igd])
            e_of_meanxy_wtflx = [mean((map_x[igd]-center_axes[0])*wtflx),$
                                 mean((map_y[igd]-center_axes[1])*wtflx)]*kpc_per_pix*ctgd
            e_of_meanr_wtflx = sqrt((mean((map_x[igd]-center_axes[0])*wtflx))^2d + $
                                    (mean((map_y[igd]-center_axes[1])*wtflx))^2d)*kpc_per_pix*ctgd
            e_of_meanpa_wtflx = ifsf_pa(mean((map_x[igd]-center_axes[0])*wtflx),$
                                        mean((map_y[igd]-center_axes[1])*wtflx))
            printf,lun_stats,icomp,'       PAflx',e_of_meanpa_wtflx,format='(I-0,A12,I5)'


;           map of outflow spaxels
            map_of_eml = bytarr(dx,dy)
            map_of_eml[igd] = 1b
            
         endif

      endfor

      if tag_exist(initmaps,'outflow_eml') then begin

         printf,lun_stats,lineofdashes
         printf,lun_stats,'CRFW model'
         printf,lun_stats,lineofdashes
         printf,lun_stats,$
                'Rmax(obs)=',string(Rmax_kpc_ha,format='(D0.2)'),' kpc'         
         printf,lun_stats,$
                'Rmax(use)=',string(initmaps.outflow_eml.R,format='(D0.2)'),' kpc'
         printf,lun_stats,$
                '      M = ',string(m_ha,format='(E0.2)'),' M_sun'
         printf,lun_stats,$
                '  dM/dt = ',string(dmdt_ha,format='(E0.2)'),' M_sun/yr'
         printf,lun_stats,$
                '      p = ',string(p_ha,format='(E0.2)'),' dyne s'
         printf,lun_stats,$
                ' dp/dt*c= ',string(dpdt_ha,format='(E0.2)'),' Lsun'
         printf,lun_stats,$
                '      E = ',string(e_ha,format='(E0.2)'),' erg'
         printf,lun_stats,$
                '  dE/dt = ',string(dedt_ha,format='(E0.2)'),' erg/s'

      endif

      linefluxes = !NULL
      if linelist.haskey('Halpha') then linefluxes = ['Halpha'] $
      else if linelist.haskey('Hbeta') then linefluxes = ['Hbeta']
      if tag_exist(initmaps,'linefluxes') then $
         linefluxes = [initmaps.linefluxes,linefluxes]
      if linefluxes ne !NULL then begin
         printf,lun_stats,lineofdashes
         printf,lun_stats,'Emission line fluxes [log(erg/s/cm^-2)]'
         printf,lun_stats,lineofdashes
         printf,lun_stats,'Line','Comp','Ext.','Unext(pp)','Unext(med)',$
                format='(5A10)'
         printf,lun_stats,lineofdashes
         foreach line,linefluxes do begin
            printf,lun_stats,line,'total',alog10(emlflxsums['tot_ext',line]),$
                   alog10(emlflxsums['tot_unext_pp',line]),$
                   alog10(emlflxsums['tot_unext_med',line]),$
                   format='(2A10,3D10.2)'
            if tag_exist(initmaps,'outflow_eml') then begin
               printf,lun_stats,line,'outflow',alog10(emlflxsums['of_ext',line]),$
                      alog10(emlflxsums['of_unext_pp',line]),$
                      alog10(emlflxsums['of_unext_med',line]),$
                      format='(2A10,3D10.2)'
            endif
         endforeach
         printf,lun_stats,lineofdashes
         printf,lun_stats,'Emission line luminosities [log(erg/s)]'
         printf,lun_stats,lineofdashes
         printf,lun_stats,'Line','Comp','Ext.','Unext(pp)','Unext(med)',$
                format='(5A10)'
         printf,lun_stats,lineofdashes
         foreach line,linefluxes do begin
            l_tot_ext = drt_linelum(emlflxsums['tot_ext',line]*1d-3,ldist,/ergs)
            l_tot_unext_pp = drt_linelum(emlflxsums['tot_unext_pp',line]*1d-3,ldist,/ergs)
            l_tot_unext_med = drt_linelum(emlflxsums['tot_unext_med',line]*1d-3,ldist,/ergs)
            printf,lun_stats,line,'total',alog10(l_tot_ext),$
                   alog10(l_tot_unext_pp),$
                   alog10(l_tot_unext_med),$
                   format='(2A10,3D10.2)'
         endforeach
      endif

   endif

   free_lun,lun_stats

   if tag_exist(initmaps,'outflow_eml') then nem = 0 else nem = -1

   windstr = {$
              nem: nem,$
;             Statistics
              a_vel_stats: a_vel_stats,$
              a_fwhm_stats: a_fwhm_stats,$
              a_v98_stats: a_v98_stats, $
              e_vel_stats: e_vel_stats,$
              e_fwhm_stats: e_fwhm_stats,$
              e_v98_stats: e_v98_stats,$
              a_vel_stats_all: a_vel_stats_all,$
              a_fwhm_stats_all: a_fwhm_stats_all,$
              a_v98_stats_all: a_v98_stats_all, $
              e_vel_stats_all: e_vel_stats_all,$
              e_fwhm_stats_all: e_fwhm_stats_all,$
              e_v98_stats_all: e_v98_stats_all,$
;             CRFW results
              a_m:m_nad,$
              a_dmdt:dmdt_nad,$
              a_p:p_nad,$
              a_dpdt:dpdt_nad,$
              a_e:e_nad,$
              a_dedt:dedt_nad,$
              a_m_err:m_err_nad,$
              a_dmdt_err:dmdt_err_nad,$
              a_p_err:p_err_nad,$
              a_dpdt_err:dpdt_err_nad,$
              a_e_err:e_err_nad,$
              a_dedt_err:dedt_err_nad,$
              a_Rmax:Rmax_kpc_nad,$
              e_m:m_ha,$
              e_dmdt:dmdt_ha,$
              e_p:p_ha,$
              e_dpdt:dpdt_ha,$
              e_e:e_ha,$
              e_dedt:dedt_ha, $
              e_m_err:m_err_ha,$
              e_dmdt_err:dmdt_err_ha,$
              e_p_err:p_err_ha,$
              e_dpdt_err:dpdt_err_ha,$
              e_e_err:e_err_ha,$
              e_dedt_err:dedt_err_ha,$
              e_Rmax:Rmax_kpc_ha,$
;             Misc
              ofcompline: ofcompline,$
              e_of_flx: e_of_flx,$
              e_total_flx: emlflxsums,$
;             Coordinates
              a_of: map_of_abs,$
              a_of_meanxy: a_of_meanxy,$
              a_of_meanr: a_of_meanr,$
              a_of_meanpa: a_of_meanpa,$
              a_of_meanxy_wtv50: a_of_meanxy_wtv50,$
              a_of_meanr_wtv50: a_of_meanr_wtv50,$
              a_of_meanpa_wtv50: a_of_meanpa_wtv50,$
              a_of_meanxy_wtv98: a_of_meanxy_wtv98,$
              a_of_meanr_wtv98: a_of_meanr_wtv98,$
              a_of_meanpa_wtv98: a_of_meanpa_wtv98,$
              a_of_meanxy_wtweq: a_of_meanxy_wtweq,$
              a_of_meanr_wtweq: a_of_meanr_wtweq,$
              a_of_meanpa_wtweq: a_of_meanpa_wtweq,$
              e_of: map_of_eml,$
              e_of_meanxy: e_of_meanxy,$
              e_of_meanr: e_of_meanr,$
              e_of_meanpa: e_of_meanpa,$
              e_of_meanxy_wtv50: e_of_meanxy_wtv50,$
              e_of_meanr_wtv50: e_of_meanr_wtv50,$
              e_of_meanpa_wtv50: e_of_meanpa_wtv50,$
              e_of_meanxy_wtv98: e_of_meanxy_wtv98,$
              e_of_meanr_wtv98: e_of_meanr_wtv98,$
              e_of_meanpa_wtv98: e_of_meanpa_wtv98,$
              e_of_meanxy_wtflx: e_of_meanxy_wtflx,$
              e_of_meanr_wtflx: e_of_meanr_wtflx,$
              e_of_meanpa_wtflx: e_of_meanpa_wtflx,$
;             Maps
              a_vel: nadabsvel, $
              a_v98: nadabsv98, $
              center_nuclei: center_nuclei,$
              center_nuclei_kpc_x: center_nuclei_kpc_x,$
              center_nuclei_kpc_y: center_nuclei_kpc_y, $
              dx: dx,$
              dy: dy,$
              e_vel: emlvel, $
              e_flx: emlflx, $
              e_flxcor_pp: emlflxcor_pp, $
              e_flxcor_med: emlflxcor_med, $
              ebv: ebv, $
              ebvmed: ebvmed,$
              errebv: errebv, $
              elecden: elecdenmap, $
              elecdenerrlo: elecdenmap_errlo, $
              elecdenerrhi: elecdenmap_errhi, $
              ibd_nadabs_fitweq: ibd_nadabs_fitweq,$
              igd_nadabs_fitweq: igd_nadabs_fitweq,$
              kpc_per_as: kpc_per_as,$
              lr: lr, $
              lrerrlo: lrerrlo, $
              lrerrhi: lrerrhi, $
              map_r: map_r,$
              map_rkpc_ifs: map_rkpc_ifs,$
              stel_vel: stel_vel,$
              stel_ebv: stel_ebv,$
              xarr_kpc: xran_kpc,$
              xran_kpc: xran_kpc,$
              yarr_kpc: yran_kpc,$
              yran_kpc: yran_kpc $
             }
   save,windstr,file=initdat.mapdir+initdat.label+'.xdr'

   save,emlvel,file=initdat.mapdir+initdat.label+'.emlvel.xdr'

   if tag_exist(initmaps,'ebv') then begin
      if tag_exist(initmaps.ebv,'calc') AND $
         tag_exist(initmaps.ebv,'apply') then begin
            save,emlflxcor_pp,file=initdat.mapdir+initdat.label+'.emlflxcor_pp.xdr'
            save,emlflxcor_med,file=initdat.mapdir+initdat.label+'.emlflxcor_med.xdr'
      endif
   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; OTHER PLOTS (GALAXY-SPECIFIC)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   plotinfo = {dx: dx,$
               dy: dy,$
               ldist: ldist,$
               map_r: map_r,$
               map_rkpc_ifs: map_rkpc_ifs,$
               map_rkpc_hst: map_rkpc_hst,$
               kpc_per_as: kpc_per_as,$
               kpc_per_pix: kpc_per_pix,$
               xran_kpc: xran_kpc,$
               yran_kpc: yran_kpc,$
               xarr_kpc: xarr_kpc,$
               yarr_kpc: yarr_kpc,$
               carr: carr,$
               center_nuclei: center_nuclei,$
               center_nuclei_kpc_x: center_nuclei_kpc_x,$
               center_nuclei_kpc_y: center_nuclei_kpc_y,$
               aspectrat: aspectrat}
;               xsec_endpoints: xsec_endpoints}

   if tag_exist(initmaps,'fcn_oplots') then begin
      if tag_exist(initmaps,'tags_oplots') then begin
         tags = initmaps.tags_oplots
         args_oplots = create_struct(tags[0], scope_varfetch(tags[0]))
         for i=1,n_elements(tags)-1 do $
            if isa(scope_varfetch(tags[i])) then $
               args_oplots = $
                  struct_addtags(args_oplots,$
                                 create_struct(tags[i],scope_varfetch(tags[i])))
         call_procedure,initmaps.fcn_oplots,initdat,initmaps,plotinfo,$
                        _extra=args_oplots
      endif else begin
         call_procedure,initmaps.fcn_oplots,initdat,initmaps,plotinfo
      endelse
   endif

end
