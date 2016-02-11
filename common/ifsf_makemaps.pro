; docformat = 'rst'
;
;+
;
; This procedure makes maps of various quantities. Contains three helper routines: 
; IFSF_PLOTRANGE, IFSF_PLOTCOMPASS, IFSF_PLOTAXESNUC.
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
;                       *cont_radprof.eps plots
;      2015sep21, DSNR, big changes to dereddening procedures, and other changes
;                       to line plotting procedures
;      2016jan24, DSNR, added 'diskline' and 'ofparline' tags to INITMAPS
;      2016feb04, DSNR, fixed treatment of HST image sizes
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
pro ifsf_plotaxesnuc,xran_kpc,yran_kpc,xnuc,ynuc,nolab=nolab
   COMPILE_OPT IDL2, HIDDEN
   if not keyword_set(nolab) then begin
      cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
      cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
   endif else begin
      cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save,xtickn=replicate(' ',60)
      cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save,ytickn=replicate(' ',60)   
   endelse
   cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
   cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
   cgoplot,xnuc,ynuc,psym=1
end
;
pro ifsf_plotcompass,xarr,yarr,carr=carr
   COMPILE_OPT IDL2, HIDDEN
   if ~ keyword_set(carr) then carr='Black'
   cgarrow,xarr[0],yarr[0],xarr[1],yarr[1],/data,/solid,color=carr
   cgarrow,xarr[0],yarr[0],xarr[2],yarr[2],/data,/solid,color=carr
   cgtext,xarr[3],yarr[3],'N',color=carr,align=0.5
   cgtext,xarr[4],yarr[4],'E',color=carr,align=0.5
end
;
function ifsf_plotrange,auto=auto,$
                        rline=rline,matline=matline,$
                        rcomp=rcomp,matcomp=matcomp,$
                        rquant=rquant,matquant=matquant,$
                        rncbdiv=rncbdiv,rlo=rlo,rhi=rhi,$
                        mapgd=mapgd,divinit=divinit,ncbdivmax=ncbdivmax
;
   if keyword_set(auto) then doauto=1 else doauto=0
   if ~ doauto then begin
      if keyword_set(rline) AND keyword_set(matline) AND $
         keyword_set(rquant) AND keyword_set(matquant) AND $
         keyword_set(rncbdiv) AND keyword_set(rlo) AND keyword_set(rhi) then begin
         if keyword_set(rcomp) AND keyword_set(matcomp) then $
            ithislinecomp = where(rline eq matline AND $
                                  rcomp eq matcomp AND $
                                  rquant eq matquant,ctthisline) $
         else $
            ithislinecomp = where(rline eq matline AND $
                                  rquant eq matquant,ctthisline)
         if ctthisline eq 1 then begin
            zran = [rlo[ithisline],rhi[ithisline]]
            ncbdiv = rncbdiv[ithisline]
            ncbdiv = ncbdiv[0]
         endif else doauto=1
      endif else doauto=1
   endif
   if doauto AND $
      keyword_set(mapgd) AND $
      keyword_set(divinit) AND $
      keyword_set(ncbdivmax) then begin
      zran = [min(mapgd),max(mapgd)]
      divarr = ifsf_cbdiv(zran,divinit,ncbdivmax)
      ncbdiv = divarr[0]
   endif else begin
      print,'IFSF_PLOTRANGE: Proper keywords not specified.'
      stop
   endelse
   return,[zran,zran[1]-zran[0],ncbdiv]
end
;
;
;-------------------------------------------------------------------------------
;
;
pro ifsf_makemaps,initproc

   fwhm2sig = 2d*sqrt(2d*alog(2d))
   plotquantum = 2.5 ; in inches
   bad = 1d99
   c_kms = 299792.458d
   ncbdivmax = 7
   maxnadabscomp = 3
   maxnademcomp = 3

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
   y = 0.9d                      ; Na ionization fraction
   a = -5.69d                    ; Na abundance
   b = -0.95d                    ; Na depletion
   volemis = 2.63d-25            ; volume emissivity of Ha = product of recomb. coeff. 
                                 ; and photon energy; units erg cm^3 s^-1
   elecden = 10                  ; electron density, cm^-3

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
   hasrangepar=0
   hasratpar_comp=0
   hasratpar_cvdf=0
   if tag_exist(initmaps,'center_axes') then $
      center_axes = initmaps.center_axes
   if tag_exist(initmaps,'center_nuclei') then $
      center_nuclei = initmaps.center_nuclei
   if tag_exist(initmaps,'comp_lr_args') then begin
      comp_lr_args = initmaps.comp_lr_args
      hasratpar_comp=1
   endif
   if tag_exist(initmaps,'cvdf_lr_args') AND $
      tag_exist(initmaps,'cvdf_lr_ftags') AND $
      tag_exist(initmaps,'cvdf_lr_ftitles') then begin
      cvdf_lr_args = initmaps.cvdf_lr_args
      cvdf_lr_ftags = initmaps.cvdf_lr_ftags
      cvdf_lr_ftitles = initmaps.cvdf_lr_ftitles
      hasratpar_cvdf=1
   endif
   if tag_exist(initmaps,'rangefile') then begin
      rangefile = initmaps.rangefile
      hasrangepar=1
   endif

;  Get linelist
   if ~ tag_exist(initmaps,'noemlinfit') then $
      linelist = ifsf_linelist(initdat.lines)
   if tag_exist(initdat,'donad') then $
      nadlinelist = ifsf_linelist(['NaD1','NaD2','HeI5876'])

;  Get range file
;
;  plot types, in order; used for correlating with input ranges (array 
;  rangequant)
   plottypes = ['flux','velocity','sigma']
   hasrangefile=0
   if hasrangepar then begin
      if file_test(rangefile) then begin
         readcol,rangefile,rangeline,rangecomp,rangequant,rangelo,rangehi,$
         rangencbdiv,format='(A,I,A,D,D,I)',/silent
         hasrangefile=1
      endif else print,'IFSF_MAKEMAPS: Range file not found.'
   endif
   if keyword_set(rangefile) then begin
      if file_test(rangefile) then begin
         readcol,rangefile,rangeline,rangecomp,rangequant,rangelo,rangehi,$
            rangencbdiv,format='(A,I,A,D,D,I)',/silent
         hasrangefile=1
      endif else print,'IFSF_MAKEMAPS: Range file not found.'
   endif

;  Restore line maps
   if ~ tag_exist(initmaps,'noemlinfit') then begin
      if not tag_exist(initdat,'outlines') then outlines = linelist->keys() $
      else outlines = initdat.outlines
      restore,file=initdat.outdir+initdat.label+'.lin.xdr'
      restore,file=initdat.outdir+initdat.label+'.tlin.xdr'

      size_tmp = size(linmaps[outlines[0]])
      dx = size_tmp[1]
      dy = size_tmp[2]
      if center_axes[0] eq -1 then center_axes = [double(dx)/2d,double(dy)/2d]
      if center_nuclei[0] eq -1 then center_nuclei = center_axes
   endif
   
;  Restore continuum parameters
   if tag_exist(initdat,'decompose_ppxf_fit') OR $
      tag_exist(initdat,'decompose_qso_fit') then begin
      restore,file=initdat.outdir+initdat.label+'.cont.xdr'
   endif

;  Get NaD parameters
   if tag_exist(initdat,'donad') then begin
      restore,file=initdat.outdir+initdat.label+'.nadspec.xdr'
      restore,file=initdat.outdir+initdat.label+'.nadfit.xdr'
      if tag_exist(initmaps,'noemlinfit') then begin
         size_tmp = size(nadcube.weq)
         dx = size_tmp[1]
         dy = size_tmp[2]
         if center_axes[0] eq -1 then center_axes = [double(dx)/2d,double(dy)/2d]
         if center_nuclei[0] eq -1 then center_nuclei = center_axes
      endif
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
      print,'IFSF_MAKEMAPS: No emission line or absorption line data specified.'
      print,'               Aborting.'
      goto,badinput
   endif
   
;  Figure aspect ratio multiplier
   if tag_exist(initmaps,'aspectrat') then aspectrat = initmaps.aspectrat $
   else aspectrat = 1d

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Compute some things
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  Luminosity and angular size distances
   ldist = lumdist(initdat.zsys_gas,H0=73,Omega_m=0.27,Lambda0=0.73,/silent)
   kpc_per_as = ldist/(1+initdat.zsys_gas^2)*1000d/206265d
   kpc_per_pix = initdat.platescale * kpc_per_as


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Load and process continuum data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  Data cube
   if not tag_exist(initdat,'datext') then datext=1 else datext=initdat.datext
   if not tag_exist(initdat,'varext') then varext=2 else varext=initdat.varext
   if not tag_exist(initdat,'dqext') then dqext=3 else dqext=initdat.dqext
   ctcube = ifsf_readcube(initdat.infile,/quiet,oned=oned,$
                          datext=datext,varext=varext,dqext=dqext)
   if tag_exist(initmaps,'fluxfactor') then begin
      ctcube.dat *= initmaps.fluxfactor
      ctcube.var *= (initmaps.fluxfactor)^2d
   endif

;  HST data
   dohst=0
   dohstbl=0
   dohstrd=0
   dohstsm=0
   dohstcol=0
   dohstcolsm=0
   if tag_exist(initmaps,'hst') AND tag_exist(initmaps,'hstbl') then begin
      dohstbl=1
      hstbl = readfits(initmaps.hstbl.file,header,/silent,/ext)
      hst_big_ifsfov = dblarr(4,2)
      if tag_exist(initmaps.hstbl,'platescale') then $
         hstpsbl = initmaps.hstbl.platescale $
      else hstpsbl = 0.05d
      if tag_exist(initmaps.hst,'subim_sm') AND $
         tag_exist(initmaps.hstbl,'sclargs_sm') then begin
         hst_sm_ifsfov = dblarr(4,2)
         bhst_sm = ifsf_hstsubim(hstbl,[initmaps.hst.subim_sm,$
                                 initmaps.hst.subim_sm],$
                                 [dx,dy],initdat.platescale,$
                                 initdat.positionangle,center_nuclei,$
                                 initmaps.hst.refcoords,$
                                 initmaps.hstbl.scllim,$
                                 sclargs=initmaps.hstbl.sclargs_sm,$
                                 ifsbounds=hst_sm_ifsfov,hstps=hstpsbl)
      endif
      bhst_big = ifsf_hstsubim(hstbl,[initmaps.hst.subim_big,$
                               initmaps.hst.subim_big],$
                               [dx,dy],initdat.platescale,$
                               initdat.positionangle,center_nuclei,$
                               initmaps.hst.refcoords,$
                               initmaps.hstbl.scllim,$
                               sclargs=initmaps.hstbl.sclargs_big,$
                               ifsbounds=hst_big_ifsfov,hstps=hstpsbl)
      bhst_fov = ifsf_hstsubim(hstbl,[0,0],[dx,dy],initdat.platescale,$
                               initdat.positionangle,center_nuclei,$
                               initmaps.hst.refcoords,$
                               initmaps.hstbl.scllim,$
                               sclargs=initmaps.hstbl.sclargs_fov,$
                               /fov,hstps=hstpsbl)
      bhst_fov_ns = ifsf_hstsubim(hstbl,[0,0],[dx,dy],initdat.platescale,$
                                  initdat.positionangle,center_nuclei,$
                                  initmaps.hst.refcoords,[0,0],/noscl,/fov,$
                                  hstps=hstpsbl)
      if tag_exist(initmaps,'hstblsm') then begin
         dohstsm=1

;        For F05189, mask central pixels before smoothing
         if initdat.label eq 'f05189' then begin
            size_tmp = size(hstbl)
            map_x_tmp = rebin(dindgen(size_tmp[1]),size_tmp[1],size_tmp[2])
            map_y_tmp = rebin(transpose(dindgen(size_tmp[2])),$
                              size_tmp[1],size_tmp[2])
            map_rkpc_tmp = sqrt((map_x_tmp - (initmaps.hst.refcoords[0]+$
                                 initmaps.hstbl.nucoffset[0]-1))^2d + $
                                (map_y_tmp - (initmaps.hst.refcoords[1]+$
                                 initmaps.hstbl.nucoffset[1]-1))^2d) $
                           * initmaps.hstbl.platescale * kpc_per_as
            ipsf = where(map_rkpc_tmp le 0.15d)
            ipsf_bkgd = where(map_rkpc_tmp gt 0.15d AND map_rkpc_tmp le 0.25d)
            hstbl_tmp = hstbl
            hstbl_tmp[ipsf] = median(hstbl[ipsf_bkgd])
            hstblsm = filter_image(hstbl_tmp,fwhm=initmaps.hst.smoothfwhm,/all)
         endif else begin
            hstblsm = filter_image(hstbl,fwhm=initmaps.hst.smoothfwhm,/all)
         endelse
         
;         print,max(hstbl[3250:3270,2700:2720]),max(hstblsm[3250:3270,2700:2720])
         
         bhst_fov_sm = ifsf_hstsubim(hstblsm,[0,0],[dx,dy],$
                                     initdat.platescale,$
                                     initdat.positionangle,center_nuclei,$
                                     initmaps.hst.refcoords,$
                                     initmaps.hstblsm.scllim,$
                                     sclargs=initmaps.hstblsm.sclargs,$
                                     /fov,hstps=hstpsbl)
         bhst_fov_sm_ns= ifsf_hstsubim(hstblsm,[0,0],[dx,dy],$
                                       initdat.platescale,$
                                       initdat.positionangle,center_nuclei,$
                                       initmaps.hst.refcoords,[0,0],/noscl,$
                                       /fov,hstps=hstpsbl)
         bhst_fov_sm_ns_rb = congrid(bhst_fov_sm_ns,dx,dy,/interp,/center)        
      endif
   endif      
   if tag_exist(initmaps,'hst') AND tag_exist(initmaps,'hstrd') then begin
      dohstrd=1
      hstrd = readfits(initmaps.hstrd.file,header,/silent,/ext)
      hst_big_ifsfov = dblarr(4,2)
      if tag_exist(initmaps.hstrd,'platescale') then $
         hstpsrd = initmaps.hstrd.platescale $
      else hstpsrd = 0.05d
      if tag_exist(initmaps.hst,'subim_sm') AND $
         tag_exist(initmaps.hstrd,'sclargs_sm') then begin
         hst_sm_ifsfov = dblarr(4,2)
         rhst_sm = ifsf_hstsubim(hstrd,[initmaps.hst.subim_sm,$
                                        initmaps.hst.subim_sm],$
                                 [dx,dy],initdat.platescale,$
                                 initdat.positionangle,center_nuclei,$
                                 initmaps.hst.refcoords,$
                                 initmaps.hstrd.scllim,$
                                 sclargs=initmaps.hstrd.sclargs_sm,$
                                 ifsbounds=hst_sm_ifsfov,hstps=hstpsrd)
      endif
      rhst_big = ifsf_hstsubim(hstrd,[initmaps.hst.subim_big,$
                               initmaps.hst.subim_big],$
                               [dx,dy],initdat.platescale,$
                               initdat.positionangle,center_nuclei,$
                               initmaps.hst.refcoords,$
                               initmaps.hstrd.scllim,$
                               sclargs=initmaps.hstrd.sclargs_big,$
                               ifsbounds=hst_big_ifsfov,hstps=hstpsrd)
      rhst_fov = ifsf_hstsubim(hstrd,[0,0],[dx,dy],initdat.platescale,$
                               initdat.positionangle,center_nuclei,$
;                               0,center_nuclei,$
                               initmaps.hst.refcoords,$
                               initmaps.hstrd.scllim,$
                               sclargs=initmaps.hstrd.sclargs_fov,$
                               /fov,hstps=hstpsrd)
      rhst_fov_ns = ifsf_hstsubim(hstrd,[0,0],[dx,dy],initdat.platescale,$
                                  initdat.positionangle,center_nuclei,$
;                                  0,center_nuclei,$
                                  initmaps.hst.refcoords,[0,0],/noscl,/fov,$
                                  hstps=hstpsrd)
      if tag_exist(initmaps,'hstrdsm') then begin
         dohstsm=1

;        For F05189, mask central pixels before smoothing
         if initdat.label eq 'f05189' then begin
            size_tmp = size(hstrd)
            map_x_tmp = rebin(dindgen(size_tmp[1]),size_tmp[1],size_tmp[2])
            map_y_tmp = rebin(transpose(dindgen(size_tmp[2])),$
                              size_tmp[1],size_tmp[2])
            map_rkpc_tmp = sqrt((map_x_tmp - (initmaps.hst.refcoords[0]+$
                                 initmaps.hstrd.nucoffset[0]-1))^2d + $
                                (map_y_tmp - (initmaps.hst.refcoords[1]+$
                                 initmaps.hstrd.nucoffset[1]-1))^2d) $
                           * initmaps.hstrd.platescale * kpc_per_as
            ipsf = where(map_rkpc_tmp le 0.15d)
            ipsf_bkgd = where(map_rkpc_tmp gt 0.15d AND map_rkpc_tmp le 0.25d)
            hstrd_tmp = hstrd
            hstrd_tmp[ipsf] = median(hstrd[ipsf_bkgd])
            hstrdsm = filter_image(hstrd_tmp,fwhm=initmaps.hst.smoothfwhm,/all)
         endif else begin
            hstrdsm = filter_image(hstrd,fwhm=initmaps.hst.smoothfwhm,/all)
         endelse

;         print,max(hstrd[3250:3270,2700:2720]),max(hstrdsm[3250:3270,2700:2720])

         rhst_fov_sm = ifsf_hstsubim(hstrdsm,[0,0],[dx,dy],$
                                     initdat.platescale,$
                                     initdat.positionangle,center_nuclei,$
                                     initmaps.hst.refcoords,$
                                     initmaps.hstrdsm.scllim,$
                                     sclargs=initmaps.hstrdsm.sclargs,$
                                     /fov,hstps=hstpsrd)
         rhst_fov_sm_ns= ifsf_hstsubim(hstrdsm,[0,0],[dx,dy],$
                                       initdat.platescale,$
                                       initdat.positionangle,center_nuclei,$
                                       initmaps.hst.refcoords,[0,0],/noscl,$
                                       /fov,hstps=hstpsrd)
         rhst_fov_sm_ns_rb = congrid(rhst_fov_sm_ns,dx,dy,/interp,/center)
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
      zprd = -2.5d*alog10(initmaps.hstrd.photflam) - 2.408d - $
             5d*alog10(initmaps.hstrd.photplam)
      zpbl = -2.5d*alog10(initmaps.hstbl.photflam) - 2.408d - $
             5d*alog10(initmaps.hstbl.photplam)

;     Take a bunch of random samples of HST image
      uplim = 0.1 ; this gets rid of cosmic rays and stars ...
      size_hst = size(hstrd)
      pxhst = round(size_hst[1]/10)
      pyhst = round(size_hst[2]/10)
      sdev = dblarr(4)
      hsttmp = hstbl[3*pxhst:4*pxhst,3*pyhst:4*pyhst]
      sdev[0] = stddev(hsttmp[where(hsttmp ne 0 AND hsttmp le uplim)])
      hsttmp = hstbl[3*pxhst:4*pxhst,6*pyhst:7*pyhst]
      sdev[1] = stddev(hsttmp[where(hsttmp ne 0 AND hsttmp le uplim)])
      hsttmp = hstbl[6*pxhst:7*pxhst,3*pyhst:4*pyhst]
      sdev[2] = stddev(hsttmp[where(hsttmp ne 0 AND hsttmp le uplim)])
      hsttmp = hstbl[6*pxhst:7*pxhst,6*pyhst:7*pyhst]
      sdev[3] = stddev(hsttmp[where(hsttmp ne 0 AND hsttmp le uplim)])
      sdevbl = median(sdev)
      hsttmp = hstrd[3*pxhst:4*pxhst,3*pyhst:4*pyhst]
      sdev[0] = stddev(hsttmp[where(hsttmp ne 0 AND hsttmp le uplim)])
      hsttmp = hstrd[3*pxhst:4*pxhst,6*pyhst:7*pyhst]
      sdev[1] = stddev(hsttmp[where(hsttmp ne 0 AND hsttmp le uplim)])
      hsttmp = hstrd[6*pxhst:7*pxhst,3*pyhst:4*pyhst]
      sdev[2] = stddev(hsttmp[where(hsttmp ne 0 AND hsttmp le uplim)])
      hsttmp = hstrd[6*pxhst:7*pxhst,6*pyhst:7*pyhst]
      sdev[3] = stddev(hsttmp[where(hsttmp ne 0 AND hsttmp le uplim)])
      sdevrd = median(sdev)
;     Find bad pixels
      colsigthr = 3d
      ibdcol = where(hstrd le colsigthr*sdevrd OR $
                     hstbl le colsigthr*sdevbl)
      hstcol = -2.5d*alog10(hstbl/hstrd) + zpbl - zprd
      hstcol[ibdcol] = 1d99
;     Extract and scale
      if tag_exist(initmaps.hst,'subim_sm') then begin
         chst_sm = ifsf_hstsubim(hstcol,[initmaps.hst.subim_sm,$
                                         initmaps.hst.subim_sm],$
                                 [dx,dy],initdat.platescale,$
                                 initdat.positionangle,center_nuclei,$
                                 initmaps.hst.refcoords,$
                                 initmaps.hstcol.scllim,$
                                 sclargs=initmaps.hstcol.sclargs,hstps=hstpsbl)
      endif
      chst_big = ifsf_hstsubim(hstcol,[initmaps.hst.subim_big,$
                               initmaps.hst.subim_big],$
                               [dx,dy],initdat.platescale,$
                               initdat.positionangle,center_nuclei,$
                               initmaps.hst.refcoords,$
                               initmaps.hstcol.scllim,$
                               sclargs=initmaps.hstcol.sclargs,hstps=hstpsbl)
      chst_fov = ifsf_hstsubim(hstcol,[0,0],[dx,dy],initdat.platescale,$
                               initdat.positionangle,center_nuclei,$
                               initmaps.hst.refcoords,$
                               initmaps.hstcol.scllim,$
                               sclargs=initmaps.hstcol.sclargs,$
                               /fov,hstps=hstpsbl)
;     Extract unscaled color image
      chst_fov_ns = ifsf_hstsubim(hstcol,[0,0],[dx,dy],initdat.platescale,$
                                  initdat.positionangle,center_nuclei,$
                                  initmaps.hst.refcoords,$
                                  initmaps.hstcol.scllim,/noscl,/fov,$
                                  hstps=hstpsbl)
   endif
   if tag_exist(initmaps,'hst') AND tag_exist(initmaps,'hstcolsm') then begin
      dohstcolsm=1
;     Take a bunch of random samples of HST image
      uplim = 0.1 ; this gets rid of cosmic rays and stars ...
      size_hst = size(hstrd)
      pxhst = round(size_hst[1]/10)
      pyhst = round(size_hst[2]/10)
      sdev = dblarr(4)
      hsttmp = hstblsm[3*pxhst:4*pxhst,3*pyhst:4*pyhst]
      sdev[0] = stddev(hsttmp[where(hsttmp ne 0 AND hsttmp le uplim)])
      hsttmp = hstblsm[3*pxhst:4*pxhst,6*pyhst:7*pyhst]
      sdev[1] = stddev(hsttmp[where(hsttmp ne 0 AND hsttmp le uplim)])
      hsttmp = hstblsm[6*pxhst:7*pxhst,3*pyhst:4*pyhst]
      sdev[2] = stddev(hsttmp[where(hsttmp ne 0 AND hsttmp le uplim)])
      hsttmp = hstblsm[6*pxhst:7*pxhst,6*pyhst:7*pyhst]
      sdev[3] = stddev(hsttmp[where(hsttmp ne 0 AND hsttmp le uplim)])
      sdevbl = median(sdev)
      hsttmp = hstrdsm[3*pxhst:4*pxhst,3*pyhst:4*pyhst]
      sdev[0] = stddev(hsttmp[where(hsttmp ne 0 AND hsttmp le uplim)])
      hsttmp = hstrdsm[3*pxhst:4*pxhst,6*pyhst:7*pyhst]
      sdev[1] = stddev(hsttmp[where(hsttmp ne 0 AND hsttmp le uplim)])
      hsttmp = hstrdsm[6*pxhst:7*pxhst,3*pyhst:4*pyhst]
      sdev[2] = stddev(hsttmp[where(hsttmp ne 0 AND hsttmp le uplim)])
      hsttmp = hstrdsm[6*pxhst:7*pxhst,6*pyhst:7*pyhst]
      sdev[3] = stddev(hsttmp[where(hsttmp ne 0 AND hsttmp le uplim)])
      sdevrd = median(sdev)
;     Find bad pixels
      colsigthr = 3d
      ibdcol = where(hstrdsm le colsigthr*sdevrd OR $
                     hstblsm le colsigthr*sdevbl)
      hstcolsm = -2.5d*alog10(hstblsm/hstrdsm) + zpbl - zprd
      hstcolsm[ibdcol] = 1d99
;     Extract and scale
      if tag_exist(initmaps.hst,'subim_sm') then begin
         cshst_sm = ifsf_hstsubim(hstcolsm,[initmaps.hst.subim_sm,$
                                            initmaps.hst.subim_sm],$
                                  [dx,dy],initdat.platescale,$
                                  initdat.positionangle,center_nuclei,$
                                  initmaps.hst.refcoords,$
                                  initmaps.hstcolsm.scllim,$
                                  sclargs=initmaps.hstcolsm.sclargs,hstps=hstpsbl)
      endif
      cshst_big = ifsf_hstsubim(hstcolsm,[initmaps.hst.subim_big,$
                                initmaps.hst.subim_big],$
                                [dx,dy],initdat.platescale,$
                                initdat.positionangle,center_nuclei,$
                                initmaps.hst.refcoords,$
                                initmaps.hstcolsm.scllim,$
                                sclargs=initmaps.hstcolsm.sclargs,hstps=hstpsbl)
      cshst_fov_s = ifsf_hstsubim(hstcolsm,[0,0],[dx,dy],initdat.platescale,$
                                initdat.positionangle,center_nuclei,$
                                initmaps.hst.refcoords,$
                                initmaps.hstcolsm.scllim,$
                                sclargs=initmaps.hstcolsm.sclargs,$
                                /fov,hstps=hstpsbl)
;     Extract unscaled color image and convert to same pixel scale as IFS data
      cshst_fov_ns = ifsf_hstsubim(hstcolsm,[0,0],[dx,dy],initdat.platescale,$
                                   initdat.positionangle,center_nuclei,$
                                   initmaps.hst.refcoords,[0,0],/noscl,/fov,$
                                   hstps=hstpsbl)
      cshst_fov_rb = congrid(cshst_fov_ns,dx,dy,/interp,/center)
   endif
   hstrd=0
   hstbl=0
   hstcol=0
   hstrdsm=0
   hstblsm=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Fit QSO PSF
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if tag_exist(initdat,'decompose_qso_fit') then begin
      qso_map = total(contcube.qso,3)
      maxqso_map = max(qso_map)
;      qso_err = stddev(contcube.qso,dim=3,/double)
      qso_err = sqrt(median(ctcube.var,dim=3,/double))
      qso_map /= maxqso_map
      qso_err /= max(median(ctcube.dat,dim=3,/double))
      
;     2D Moffat fit to continuum flux vs. radius
      est=[0d,1d,2d,2d,center_nuclei[0]-1d,center_nuclei[1]-1d,0d,2.5d]
      qso_fit = mpfit2dpeak(qso_map,qso_fitpar,/moffat,/circular,est=est,$
                            error=qso_err)
   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Compute some more things
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

;  Radii in kpc
;  GMOS FOV
   map_x = rebin(dindgen(dx)+1,dx,dy)
   map_y = rebin(transpose(dindgen(dy)+1),dx,dy)
   map_r = sqrt((map_x - center_axes[0])^2d + (map_y - center_axes[1])^2d)
   map_rkpc_ifs = map_r * kpc_per_pix
   if tag_exist(initdat,'decompose_qso_fit') then begin
      map_rnuc = sqrt((map_x - qso_fitpar[4]+1)^2d + $
                      (map_y - qso_fitpar[5]+1)^2d)
      map_rnuckpc_ifs = map_rnuc * kpc_per_pix
      psf1d_x = dindgen(101)/100d*max(map_rnuckpc_ifs)
      psf1d_y = alog10(moffat(psf1d_x,[qso_fitpar[1],$
                                       0d,$
                                       qso_fitpar[2]*kpc_per_pix,$
                                       qso_fitpar[7]]))      
   endif
;  HST FOV
   if (dohstrd OR dohstbl) then begin
      if dohstbl then size_subim = size(bhst_fov) $
      else size_subim = size(rhst_fov)
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
         if dohstrd then begin
            hstplatescale = initmaps.hstbl.platescale
            kpc_per_hstpix = hstplatescale * kpc_per_as
            if initmaps.hstbl.platescale ne initmaps.hstrd.platescale then begin
               print,'WARNING: HST blue and red plate scales differ;'
               print,'         using blue platescale for radius calculations.'
            endif
         endif
;        Radius of each HST pixel from axis [0,0] point, in kpc
         map_rkpc_hst = map_r_hst * initmaps.hstbl.platescale * kpc_per_as
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
   endif

;  Sort emission line components if requested
   if tag_exist(initmaps,'fcnsortcomp') AND $
      tag_exist(initmaps,'sortlines') AND $
      tag_exist(initmaps,'sorttype') then begin
      if tag_exist(initmaps,'argssortcomp') then $
         linmaps = call_function(initmaps.fcnsortcomp,dx,dy,linmaps,$
                                 initdat.linetie,initmaps.sortlines,$
                                 initmaps.sorttype,$
                                 _extra=initmaps.argssortcomp) $
      else $
         linmaps = call_function(initmaps.fcnsortcomp,dx,dy,linmaps,$
                                 initdat.linetie,initmaps.sortlines,$
                                 initmaps.sorttype)
   endif

;  Apply a sigma cut on a line-by-line basis. Previous sigma cuts may have been 
;  less strict (see, e.g., IFSF_SORTCOMP).
   if tag_exist(initmaps,'sigthresh') then begin
      foreach line,outlines do begin
         flux = linmaps[line,*,*,*,0]
         fluxerr = linmaps[line,*,*,*,1]
         wave = linmaps[line,*,*,*,2]
         sig = linmaps[line,*,*,*,3]
         fluxpk = linmaps[line,*,*,*,4]
         tflux = tlinmaps[line,*,*,*,0]
         tfluxerr = tlinmaps[line,*,*,*,1]
         ibd = where(flux gt 0d AND $
                     flux ne bad AND $
                     flux lt fluxerr*initmaps.sigthresh,ctbd)
         if ctbd gt 0 then begin
            flux[ibd] = bad
            fluxerr[ibd] = bad
            wave[ibd] = bad
            sig[ibd] = bad
            fluxpk[ibd] = bad
            linmaps[line,*,*,*,0] = flux
            linmaps[line,*,*,*,1] = fluxerr
            linmaps[line,*,*,*,2] = wave
            linmaps[line,*,*,*,3] = sig
            linmaps[line,*,*,*,4] = fluxpk
         endif
         tibd = where(tflux gt 0d AND $
                      tflux ne bad AND $
                      tflux lt tfluxerr*initmaps.sigthresh,tctbd)
         if tctbd gt 0 then begin
            tflux[ibd] = bad
            tfluxerr[ibd] = bad
            tlinmaps[line,*,*,*,0] = tflux
            tlinmaps[line,*,*,*,1] = tfluxerr
         endif
      endforeach
   endif
   
;  Emission line maps and parameters
   
   if ~ tag_exist(initdat,'noemlinfit') then begin
      linspecmaps = hash()  ; line maps, possibly dereddened
      linspecpars = hash() ; CVDFs, possibly dereddened
      linspecpars_arr = hash() ; CVDFs in array form, possibly dereddened
      ofpars = hash() ; outflow CVDF, possibly dereddened
;     Set reference line for rotation curve, for defining outflows, by resorting
;     OUTLINES so that reference line comes first.
      if ~ tag_exist(initmaps,'diskline') then diskline='Halpha' $
      else diskline = initmaps.diskline
      idiskline = outlines.where(diskline,count=ctdiskline)
      if ctdiskline gt 0 then outlines.move,idiskline,0 $
      else print,'IFSF_MAKEMAPS: Disk line not found; using first element in ',$
                 'OUTLINES list.'
      foreach line,outlines do begin
;        Line maps
         linspecmaps[line] = $
            ifsf_cmplinspecmaps(linmaps[line,*,*,*,4],$
                                linmaps[line,*,*,*,2],$
                                linmaps[line,*,*,*,3],$
                                initdat.maxncomp,linelist[line],$
                                initdat.zsys_gas)
                        
;        Cumulative velocity distribution functions, computed using extincted
;        fluxes and, optionally, fluxes dereddened using component extinctions                                
         linpararr=1b
         ofpars_line=0b
         sigthresh=0b
         ofthresh=0b
         ofignore=0b
         diffthresh=0b
         if tag_exist(initmaps,'compof') then begin
            ofpars_line=1b
            if tag_exist(initmaps,'compof_sigthresh') then $
               sigthresh=initmaps.compof_sigthresh
            if tag_exist(initmaps,'compof_ofthresh') then $
               ofthresh=initmaps.compof_ofthresh
            if tag_exist(initmaps,'compof_diffthresh') then $
               diffthresh=initmaps.compof_diffthresh
            if tag_exist(initmaps,'compof_ignore') then $
               ofignore=initmaps.compof_ignore
         endif
         if line eq diskline then diskrot=0b $
         else diskrot=linspecpars[diskline].vpk
         linspecpars[line] = $
            ifsf_cmplinspecpars(linspecmaps[line],linpararr=linpararr,$
                                ofpars=ofpars_line,$
                                sigthresh=sigthresh,ofthresh=ofthresh,$
                                ofignore=ofignore,diffthresh=diffthresh,$
                                diskrot=diskrot)
         if tag_exist(initmaps,'compof') then ofpars[line] = ofpars_line
         linspecpars_arr[line] = linpararr
      endforeach         
      linspecpars_tags = tag_names(linspecpars[outlines[0]])
      if tag_exist(initmaps,'compof') then $
         ofpars_tags = tag_names(ofpars[outlines[0]])
   endif
        
;  Line ratios by component, with errors
   if hasratpar_comp then linrats = ifsf_lineratios(linmaps,linelist)

;  Line ratios summed over all components, with errors 
   if hasratpar_comp then tlinrats = ifsf_lineratios(tlinmaps,linelist)

;  Line ratios by CVDF, without errors
   if hasratpar_cvdf then $
      linrats_cvdf = ifsf_lineratios(linspecpars_arr,linelist,/noerr)

; Apply median filter to E(B-V) maps?
   if tag_exist(initmaps,'ebv_medfilt') then begin
      for i=0,4 do begin
         tmpmap = $
            filter_image(linrats_cvdf['ebv',*,*,i],$
                         median=initmaps.ebv_medfilt,$
                         /all_pixels)
         ibd = where(tmpmap eq bad,ctbd)
         if ctbd gt 0 then tmpmap[ibd] = bad
         linrats_cvdf['ebv',*,*,i] = tmpmap
      endfor
   endif

;  Unextincted emission line maps and parameters
   if ~ tag_exist(initmaps,'noemlinfit') AND $
      (tag_exist(initmaps,'applyebv') OR $
       tag_exist(initmaps,'applyebv_tot') OR $
       tag_exist(initmaps,'applyebv_single')) then begin
      linspecmaps_ext = linspecmaps ; extincted line maps
      linspecmaps = hash()
      linspecpars_ext = linspecpars ; extincted CVDFs
      linspecpars = hash()
      linspecpars_ext_arr = linspecpars_arr ; extincted CVDFs in array form
      linspecpars_arr = hash()
      ofpars_ext = ofpars ; extincted outflow CVDF
      ofpars = hash()
;     applyebv: Apply extinction calculated in each component to that component.
;        Set tag to array with number of elements equal to maxncomp; each 
;        element specifies whether or not extinction applied to that component.
;     applyebv_tot: Apply extinction calculated using total line flux to the
;        entire line. Set tag to byte flag.
      if tag_exist(initmaps,'applyebv_single') then begin
         tlinrats_tmp = tlinrats['ebv',*,*,0]
         igd = where(tlinrats_tmp ne bad)
         medebv = median(tlinrats_tmp[igd])
         tlinrats_tmp = dblarr(dx,dy)+bad
         tlinrats_tmp[igd] = medebv
         doebv_tmp = rebin(tlinrats_tmp,dx,dy,initdat.maxncomp)
      endif else if tag_exist(initmaps,'applyebv_tot') then $
         doebv_tmp = rebin(tlinrats['ebv',*,*,0],dx,dy,initdat.maxncomp) $
      else $
         doebv_tmp = $
            linrats['ebv'] * $
            double(rebin(reform(initmaps.applyebv,1,1,initdat.maxncomp),$
                         dx,dy,initdat.maxncomp))
      foreach line,outlines do begin
;        apply E(B-V) only to those components that are specified
;        and tied to Halpha or to the same line as Halpha
         if (initdat.linetie)[line] eq 'Halpha' OR $
            (initdat.linetie)[line] eq (initdat.linetie)['Halpha'] then $
               doebv = doebv_tmp $
         else doebv=0b
;        Line maps
         linspecmaps[line] = $
            ifsf_cmplinspecmaps(linmaps[line,*,*,*,4],$
                                linmaps[line,*,*,*,2],$
                                linmaps[line,*,*,*,3],$
                                initdat.maxncomp,linelist[line],$
                                initdat.zsys_gas,ebv=doebv)

         linpararr=1b
         ofpars_line=0b
         sigthresh=0b
         ofthresh=0b
         ofignore=0b
         diffthresh=0b
         if tag_exist(initmaps,'compof') then begin
            ofpars_line=1b
            if tag_exist(initmaps,'compof_sigthresh') then $
               sigthresh=initmaps.compof_sigthresh
            if tag_exist(initmaps,'compof_ofthresh') then $
               ofthresh=initmaps.compof_ofthresh
            if tag_exist(initmaps,'compof_diffthresh') then $
               diffthresh=initmaps.compof_diffthresh
            if tag_exist(initmaps,'compof_ignore') then $
               ofignore=initmaps.compof_ignore
         endif
         if line eq diskline then diskrot=0b $
         else diskrot=linspecpars[diskline].vpk
         linspecpars[line] = $
            ifsf_cmplinspecpars(linspecmaps[line],linpararr=linpararr,$
                                ofpars=ofpars_line,$
                                sigthresh=sigthresh,ofthresh=ofthresh,$
                                ofignore=ofignore,diffthresh=diffthresh,$
                                diskrot=diskrot)
         if tag_exist(initmaps,'compof') then ofpars[line] = ofpars_line
         linspecpars_arr[line] = linpararr
      endforeach         
   endif

   if tag_exist(initdat,'donad') then begin
;     Compute velocities and column densities of NaD model fits
      nadabscf = dblarr(dx,dy,maxnadabscomp+1)+bad
      errnadabscf = dblarr(dx,dy,maxnadabscomp+1,2)+bad
      nadabstau = dblarr(dx,dy,maxnadabscomp+1)+bad
      errnadabstau = dblarr(dx,dy,maxnadabscomp+1,2)+bad
      nadabsvel = dblarr(dx,dy,maxnadabscomp+1)+bad
      errnadabsvel = dblarr(dx,dy,maxnadabscomp+1)+bad
      nademvel = dblarr(dx,dy,maxnademcomp+1)+bad
      errnademvel = dblarr(dx,dy,maxnademcomp+1)+bad
      nadabssig = dblarr(dx,dy,maxnadabscomp+1)+bad
      errnadabssig = dblarr(dx,dy,maxnadabscomp+1)+bad
      nademsig = dblarr(dx,dy,maxnademcomp+1)+bad
      errnademsig = dblarr(dx,dy,maxnademcomp+1)+bad
      nadabsv98 = dblarr(dx,dy,maxnadabscomp+1)+bad
      errnadabsv98 = dblarr(dx,dy,maxnadabscomp+1)+bad
      nademv98 = dblarr(dx,dy,maxnademcomp+1)+bad
      errnademv98 = dblarr(dx,dy,maxnademcomp+1)+bad
      nadabsnh = dblarr(dx,dy)+bad
      nadabscnh = dblarr(dx,dy,maxnadabscomp+1)+bad
      errnadabsnh = dblarr(dx,dy,2)+bad
      errnadabscnh = dblarr(dx,dy,maxnadabscomp+1,2)+bad
      nadabsncomp = intarr(dx,dy)+bad
      nademncomp = intarr(dx,dy)+bad
      for i=0,dx-1 do begin
         for j=0,dy-1 do begin
            igd = where(nadfit.waveabs[i,j,*] ne bad AND $
                        nadfit.waveabs[i,j,*] ne 0,ctgd)
            if ctgd gt 0 then begin
               nnai = total(nadfit.tau[i,j,igd]*nadfit.sigmaabs[i,j,igd]) / $
                      (1.497d-15/sqrt(2d)*nadlinelist['NaD1']*0.3180d)
               tauerrlo = nadfit.tau[i,j,igd] - $
                          10d^(alog10(nadfit.tau[i,j,igd]) - $
                               nadfit.tauerr[i,j,igd,0])
               tauerrhi = 10d^(alog10(nadfit.tau[i,j,igd]) + $
                               nadfit.tauerr[i,j,igd,1]) - nadfit.tau[i,j,igd]
               nnaierrlo = $
                  nnai*sqrt(total((tauerrlo/$
                                   nadfit.tau[i,j,igd])^2d + $
                                  (nadfit.sigmaabserr[i,j,igd,0]/$
                                   nadfit.sigmaabs[i,j,igd])^2d))
               nnaierrhi = $
                  nnai*sqrt(total((tauerrhi/$
                                   nadfit.tau[i,j,igd])^2d + $
                                  (nadfit.sigmaabserr[i,j,igd,1]/$
                                   nadfit.sigmaabs[i,j,igd])^2d))
               nadabsnh[i,j] = nnai/(1-0.9d)/10^(-5.69d - 0.95d)
               errnadabsnh[i,j,*] = [nnaierrlo,nnaierrhi]/$
                                    (1-0.9d)/10^(-5.69d - 0.95d)
               nadabsncomp[i,j] = ctgd
;              For now, average lo/hi errors in wavelength and sigma
               tmptau=nadfit.tau[i,j,igd]
               tmptauerrlo = nadfit.tauerr[i,j,*,0]
               tmptauerrhi = nadfit.tauerr[i,j,*,1]
               tmpcf=nadfit.cf[i,j,igd]
               tmpcferrlo = nadfit.cferr[i,j,*,0]
               tmpcferrhi = nadfit.cferr[i,j,*,1]
               tmpwaveabs = nadfit.waveabs[i,j,igd]
               tmpwaveabserr = mean(nadfit.waveabserr[i,j,*,*],dim=4)
               tmpwaveabserr = tmpwaveabserr[igd]
               tmpsigabs=nadfit.sigmaabs[i,j,igd]
               tmpsigabserr = mean(nadfit.sigmaabserr[i,j,*,*],dim=4)
               tmpsigabserr = tmpsigabserr[igd]
            endif
;           Sort absorption line wavelengths. In output arrays,
;           first element of third dimension holds data for spaxels with only
;           1 component. Next elements hold velocities for spaxels with more than
;           1 comp, in order of increasing blueshift. Formula for computing error
;           in velocity results from computing derivative in Wolfram Alpha w.r.t.
;           lambda and rearranging on paper.
            if ctgd eq 1 then begin
               nadabscnh[i,j,0]=nadabsnh[i,j]
               errnadabscnh[i,j,0,*]=errnadabsnh[i,j,*]
               zdiff = tmpwaveabs/nadlinelist['NaD1']-1d - initdat.zsys_gas
               nadabsvel[i,j,0] = c_kms * ((zdiff+1d)^2d - 1d) / $
                                  ((zdiff+1d)^2d + 1d)
               errnadabsvel[i,j,0] = c_kms * (4d/nadlinelist['NaD1']*(zdiff+1d)/$
                                     ((zdiff+1d)^2d + 1d)^2d) * tmpwaveabserr
               nadabssig[i,j,0] = tmpsigabs
               errnadabssig[i,j,0] = tmpsigabserr
               nadabsv98[i,j,0] = nadabsvel[i,j,0]-2d*nadabssig[i,j,0]
               errnadabsv98[i,j,0] = $
                  sqrt(errnadabsvel[i,j,0]^2d + 4d*errnadabssig[i,j,0]^2d)
               nadabscf[i,j,0] = tmpcf
               errnadabscf[i,j,0,0:1] = [tmpcferrlo[igd],tmpcferrhi[igd]]
               nadabstau[i,j,0] = tmptau
               errnadabstau[i,j,0,0:1] = [tmptauerrlo[igd],tmptauerrhi[igd]]
            endif else if ctgd gt 1 then begin
               sortgd = sort(tmpwaveabs)
               nnai = reverse(tmptau[sortgd]*tmpsigabs[sortgd]) / $
                      (1.497d-15/sqrt(2d)*nadlinelist['NaD1']*0.3180d)
               nadabscnh[i,j,1:ctgd] = nnai/(1-0.9d)/10^(-5.69d - 0.95d)
               nnaierrlo = $
                  nnai*sqrt(reverse((tmptauerrlo[sortgd]/tmptau[sortgd])^2d + $
                                    (tmpsigabserr[sortgd]/tmpsigabs[sortgd])^2d))
               nnaierrhi = $
                  nnai*sqrt(reverse((tmptauerrhi[sortgd]/tmptau[sortgd])^2d + $
                                    (tmpsigabserr[sortgd]/tmpsigabs[sortgd])^2d))
               errnadabscnh[i,j,1:ctgd,0] = nnaierrlo / $
                                            (1-0.9d)/10^(-5.69d - 0.95d)
               errnadabscnh[i,j,1:ctgd,1] = nnaierrhi / $
                                            (1-0.9d)/10^(-5.69d - 0.95d)
               zdiff = reverse(tmpwaveabs[sortgd])/nadlinelist['NaD1']-1d - $
                       initdat.zsys_gas
               nadabsvel[i,j,1:ctgd] = c_kms * ((zdiff+1d)^2d - 1d) / $
                                       ((zdiff+1d)^2d + 1d)
               errnadabsvel[i,j,1:ctgd] = $
                  c_kms * (4d/nadlinelist['NaD1']*(zdiff+1d)/$
                  ((zdiff+1d)^2d + 1d)^2d) * reverse(tmpwaveabserr[sortgd])
               nadabssig[i,j,1:ctgd] = reverse(tmpsigabs[sortgd])
               errnadabssig[i,j,1:ctgd] = reverse(tmpsigabserr[sortgd])
               nadabsv98[i,j,1:ctgd] = nadabsvel[i,j,1:ctgd]-$
                                       2d*nadabssig[i,j,1:ctgd]
               errnadabsv98[i,j,1:ctgd] = $
                  sqrt(errnadabsvel[i,j,1:ctgd]^2d + $
                  4d*errnadabssig[i,j,1:ctgd]^2d)
               nadabscf[i,j,1:ctgd] = reverse(tmpcf[sortgd])
               errnadabscf[i,j,1:ctgd,0] = reverse(tmpcferrlo[sortgd])
               errnadabscf[i,j,1:ctgd,1] = reverse(tmpcferrhi[sortgd])
               nadabstau[i,j,1:ctgd] = reverse(tmptau[sortgd])
               errnadabstau[i,j,1:ctgd,0] = reverse(tmptauerrlo[sortgd])
               errnadabstau[i,j,1:ctgd,1] = reverse(tmptauerrhi[sortgd])
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
               zdiff = tmpwaveem/nadlinelist['NaD1']-1d - initdat.zsys_gas
               nademvel[i,j,0] = c_kms * ((zdiff+1d)^2d - 1d) / $
                                 ((zdiff+1d)^2d + 1d)
               errnademvel[i,j,0] = c_kms * (4d/nadlinelist['NaD1']*(zdiff+1d)/$
                                    ((zdiff+1d)^2d + 1d)^2d) * tmpwaveemerr
               nademsig[i,j,0] = tmpsigem
               errnadabssig[i,j,0] = tmpsigemerr
               nademv98[i,j,0] = nademvel[i,j,0]+2d*nademsig[i,j,0]
               errnademv98[i,j,0] = $
                  sqrt(errnademvel[i,j,0]^2d + 4d*errnademsig[i,j,0]^2d)
            endif else if ctgd gt 1 then begin
               sortgd = sort(tmpwaveem)
               zdiff = tmpwaveem[sortgd]/nadlinelist['NaD1']-1d -initdat.zsys_gas
               nademvel[i,j,1:ctgd] = c_kms * ((zdiff+1d)^2d - 1d) / $
                                      ((zdiff+1d)^2d + 1d)               
               errnademvel[i,j,1:ctgd] = $
                  c_kms * (4d/nadlinelist['NaD1']*(zdiff+1d)/$
                  ((zdiff+1d)^2d + 1d)^2d) * $
                  tmpwaveemerr[sortgd]
               nademsig[i,j,1:ctgd] = tmpsigem[sortgd]
               errnademsig[i,j,1:ctgd] = tmpsigemerr[sortgd]
               nademv98[i,j,1:ctgd] = nademvel[i,j,1:ctgd]+2d*nademsig[i,j,1:ctgd]
               errnademv98[i,j,1:ctgd] = $
                  sqrt(errnademvel[i,j,1:ctgd]^2d + $
                  4d*errnademsig[i,j,1:ctgd]^2d)                                       
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
      cvd_nad_maps = ifsf_cmplinspecmaps(nadfit.tau,nadfit.waveabs,$
                                         nadfit.sigmaabs,initnad.maxncomp,$
                                         nadlinelist['NaD1'],$
                                         initdat.zsys_gas)
      cvd_nad_pars = ifsf_cmplinspecpars(cvd_nad_maps)

      
   endif

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Fit PSF to Emission Line Map
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if tag_exist(initmaps,'fit_empsf') then begin
      linmap_tmp = linspecmaps[initmaps.fit_empsf.line]
      vel = initmaps.fit_empsf.vel
      ivel = value_locate(linmap_tmp.vel,vel)
      empsf_map = linmap_tmp.flux[*,*,ivel]
      maxempsf_map = max(empsf_map,imax)
      empsf_map /= maxempsf_map

;     Use error in total flux for error in line
      empsf_err = tlinmaps[initmaps.fit_empsf.line,*,*,0,1]
      empsf_err /= empsf_err[imax]

;     2D Moffat fit to continuum flux vs. radius
      parinfo = REPLICATE({fixed:0b},8)
      parinfo[0].fixed = 1b
      est=[0d,1d,2d,2d,center_nuclei[0]-1d,center_nuclei[1]-1d,0d,2.5d]
      empsf_fit = $
         mpfit2dpeak(empsf_map,empsf_fitpar,/moffat,/circular,est=est,$
                     parinfo=parinfo,error=empsf_err)

      map_rempsf = sqrt((map_x - empsf_fitpar[4]+1)^2d + $
                        (map_y - empsf_fitpar[5]+1)^2d)
      map_rempsfkpc_ifs = map_rempsf * kpc_per_pix
      empsf1d_x = dindgen(101)/100d*max(map_rempsfkpc_ifs)
      empsf1d_y = alog10(moffat(empsf1d_x,[empsf_fitpar[1],$
                                           0d,$
                                           empsf_fitpar[2]*kpc_per_pix,$
                                           empsf_fitpar[7]]))      

   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Continuum plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if dohstrd OR dohstbl then begin
      npx = 2
      npy = 2
   endif else begin
      npx = 1
      npy = 1
   endelse
   if dohstrd OR dohstbl then begin
      if dohstrd AND dohstbl then begin
         cap1 = textoidl(initmaps.hstbl.label+'+'+initmaps.hstrd.label)
         cap2 = textoidl(initmaps.hstbl.label+'+'+initmaps.hstrd.label)
         cap3 = textoidl(initmaps.hstbl.label+'+'+initmaps.hstrd.label)
         cap4 = textoidl('IFS cont.')
      endif else if dohstrd then begin
         cap1 = textoidl(initmaps.hstrd.label)
         cap2 = textoidl(initmaps.hstrd.label)
         cap3 = textoidl(initmaps.hstrd.label)
         cap4 = textoidl('IFS cont.')
      endif else begin
         cap1 = textoidl(initmaps.hstbl.label)
         cap2 = textoidl(initmaps.hstbl.label)
         cap3 = textoidl(initmaps.hstbl.label)
         cap4 = textoidl('IFS cont.')
      endelse
   endif else cap4 = textoidl('IFS cont.')
;  arrays for positions for zoom box
   posbox1x = dblarr(2)
   posbox1y = dblarr(2)
   posbox2x = dblarr(2)
   posbox2y = dblarr(2)

   cgps_open,initdat.mapdir+initdat.label+'cont.eps',charsize=1,/encap,$
             /inches,xs=plotquantum*npx,ys=plotquantum*npy*aspectrat,/qui
   pos = cglayout([npx,npy],ixmar=[2d,2d],iymar=[2d,2d],$
                  oxmar=[1,0],oymar=[0,0],xgap=0,ygap=0,unit=!D.X_PX_CM/3.0)

   if (dohstrd OR dohstbl) then begin
      i = 0
;     HST continuum, large scale
      if dohstbl then size_subim = size(bhst_big) $
      else size_subim = size(rhst_big)
      if dohstrd AND dohstbl then begin
         mapscl = bytarr(3,size_subim[1],size_subim[2])
         mapscl[0,*,*] = rhst_big
         mapscl[2,*,*] = bhst_big
         mapscl[1,*,*] = byte((double(rhst_big)+double(bhst_big))/2d)
         mapscl = rebin(mapscl,3,size_subim[1]*10,size_subim[2]*10,/sample)
      endif else begin
         mapscl = bytarr(size_subim[1],size_subim[2])
         if dohstrd then mapscl = rhst_big
         if dohstbl then mapscl = bhst_big
         mapscl = rebin(mapscl,size_subim[1]*10,size_subim[2]*10,/sample)
      endelse
      cgloadct,65,/reverse
      cgimage,mapscl,$
              /keep,pos=pos[*,i],opos=truepos,$
              noerase=i ne 0,missing_value=bad,missing_index=255,$
              missing_color='white'
      cgplot,[0],xsty=5,ysty=5,xran=[0,size_subim[1]],$
             yran=[0,size_subim[2]],position=truepos,$
             /nodata,/noerase,title=cap1
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

      i = 1
;     HST continuum, IFS FOV
      if dohstbl then size_subim = size(bhst_fov) $
      else size_subim = size(rhst_fov)
      if dohstbl AND dohstrd then begin
         mapscl = bytarr(3,size_subim[1],size_subim[2])
         mapscl[0,*,*] = rhst_fov
         mapscl[2,*,*] = bhst_fov
         mapscl[1,*,*] = byte((double(rhst_fov)+double(bhst_fov))/2d)
         ctmap = (rhst_fov_ns+bhst_fov_ns)/2d
         mapscl = rebin(mapscl,3,size_subim[1]*10,size_subim[2]*10,/sample)
      endif else begin
         mapscl = bytarr(size_subim[1],size_subim[2])
         if dohstrd then begin
            mapscl = rhst_fov
            ctmap = rhst_fov_ns
         endif else if dohstbl then begin
            mapscl = bhst_fov
            ctmap = bhst_fov_ns
         endif
         mapscl = rebin(mapscl,size_subim[1]*10,size_subim[2]*10,/sample)
      endelse

      
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos[*,i],opos=truepos,$
              noerase=i ne 0,missing_value=bad,missing_index=255,$
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
;        Fitted peak coordinate in HST pixels, in single-offset coordinates
         peakfit_hstpix = [a[4]+xhst_sub[0]+1,a[5]+yhst_sub[0]+1]
         peakfit_hst_distance_from_nucleus_hstpix = peakfit_hstpix - $
                                                    center_nuclei_hst[*,0]
         peakfit_hst_distance_from_nucleus_kpc = $
            peakfit_hst_distance_from_nucleus_hstpix * kpc_per_hstpix
         size_hstpix = size(ctmap)
         cgplot,[0],xsty=5,ysty=5,xran=[1,size_hstpix[1]],$
                yran=[1,size_hstpix[2]],position=truepos,$
                /nodata,/noerase,title=cap2
         cgoplot,peakfit_hstpix[0],peakfit_hstpix[1],psym=1,color='Red'         
      endif else begin
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase,title=cap2         
      endelse
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
      cgtext,xran_kpc[0]+(xran_kpc[1]-xran_kpc[0])*0.05,$
             yran_kpc[1]-(yran_kpc[1]-yran_kpc[0])*0.1,$
             'IFS FOV',/data,color='white'
      ifsf_plotcompass,xarr_kpc,yarr_kpc,carr=carr
      posbox1x[1] = truepos[0]
      posbox1y[1] = truepos[3]
      posbox2x[1] = truepos[0]
      posbox2y[1] = truepos[1]

      i = 2
;     smoothed HST continuum, IFS FOV
      if dohstsm then begin
;;        3-color image
;         mapscl = bytarr(3,size_subim[1],size_subim[2])
;         if dohstrd then mapscl[0,*,*] = rhst_fov_sm
;         if dohstbl then mapscl[2,*,*] = bhst_fov_sm
;         if dohstrd AND dohstbl then $
;            mapscl[1,*,*] = byte((double(rhst_fov_sm)+double(bhst_fov_sm))/2d)
;        Flux image
         if dohstbl AND dohstrd then $
            ctmap = (double(rhst_fov_sm_ns_rb)+double(bhst_fov_sm_ns_rb))/2d $
         else if dohstbl then ctmap = bhst_fov_sm_ns_rb $
         else ctmap = rhst_fov_sm_ns_rb
         ctmap /= max(ctmap)
         zran = initmaps.ct.scllim
         dzran = zran[1]-zran[0]
         mapscl = cgimgscl(rebin(ctmap,dx*20,dy*20,/sample),$
                           minval=zran[0],max=zran[1],$
                           stretch=initmaps.ct.stretch)
         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos[*,i],opos=truepos,$
                 noerase=i ne 0,missing_value=bad,missing_index=255,$
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
;           Fitted peak coordinate in IFS pixels, in single-offset coordinates
            peakfit_pix = [a[4]+x_sub[0]+1,a[5]+y_sub[0]+1]
            peakfit_hstconv_distance_from_nucleus_pix = peakfit_pix - $
                                                        center_nuclei[*,0]
            peakfit_hstconv_distance_from_nucleus_kpc = $
                peakfit_hstconv_distance_from_nucleus_pix * kpc_per_pix
            cgplot,[0],xsty=5,ysty=5,xran=[1,dx],$
                   yran=[1,dy],position=truepos,$
                   /nodata,/noerase,title=cap3
            cgoplot,peakfit_pix[0],peakfit_pix[1],psym=1,color='Red'         
         endif else begin
            cgplot,[0],xsty=5,ysty=5,position=truepos,$
                   /nodata,/noerase,title=cap3        
         endelse
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
         cgtext,xran_kpc[0]+(xran_kpc[1]-xran_kpc[0])*0.05,$
                yran_kpc[1]-(yran_kpc[1]-yran_kpc[0])*0.1,$
                'IFS FOV, conv.',/data,color='white'
      endif
 
      i=3
 
   endif else i=0

   if tag_exist(initmaps,'ct') then begin
      ictlo = value_locate(ctcube.wave,initmaps.ct.sumrange[0])
      icthi = value_locate(ctcube.wave,initmaps.ct.sumrange[1])
      zran = initmaps.ct.scllim
      dzran = zran[1]-zran[0]
      ctmap = total(ctcube.dat[*,*,ictlo:icthi],3)
      ctmap /= max(ctmap)
      ctmap_save = ctmap
      mapscl = cgimgscl(rebin(ctmap,dx*20,dy*20,/sample),$
                        minval=zran[0],max=zran[1],stretch=initmaps.ct.stretch)
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos[*,i],opos=truepos,$
              noerase=i ne 0,missing_value=bad,missing_index=255,$
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
;           Fitted peak coordinate in IFS pixels, in single-offset coordinates
            peakfit_pix = [a[4]+x_sub[0]+1,a[5]+y_sub[0]+1]
            peakfit_pix_ifs = peakfit_pix ; save for later
            peakfit_ifs_distance_from_nucleus_pix = peakfit_pix - $
                                                    center_nuclei[*,0]
            peakfit_ifs_distance_from_nucleus_kpc = $
                peakfit_ifs_distance_from_nucleus_pix * kpc_per_pix
            cgplot,[0],xsty=5,ysty=5,xran=[1,dx],$
                   yran=[1,dy],position=truepos,$
                   /nodata,/noerase,title=cap4
            cgoplot,peakfit_pix[0],peakfit_pix[1],psym=1,color='Red'         
      endif else begin
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase,title=cap4    
      endelse
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
   endif

   cgplot,posbox1x,posbox1y,color='Red',$
          xsty=5,ysty=5,/noerase,xran=[0,1],yran=[0,1],pos=[0,0,1,1]
   cgoplot,posbox2x,posbox2y,color='Red'

   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Continuum color plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if (dohstrd AND dohstbl) then begin

      npx = 2
      npy = 2
      cap1 = textoidl(initmaps.hstbl.label+'-'+initmaps.hstrd.label)
      cap2 = textoidl(initmaps.hstbl.label+'-'+initmaps.hstrd.label)
      cap3 = textoidl(initmaps.hstbl.label+'-'+initmaps.hstrd.label)
      cap4 = textoidl('IFS cont.')
      cbform='(D0.1)'

      cgps_open,initdat.mapdir+initdat.label+'color.eps',charsize=1,/encap,$
                /inches,xs=plotquantum*npx,ys=plotquantum*npy,/qui
      pos = cglayout([npx,npy],ixmar=[2,2],iymar=[2,2],oxmar=[0,0],oymar=[0,0],$
                     xgap=0,ygap=0,unit=!D.X_PX_CM/3.0)

      i = 0
      size_subim = size(chst_big)
      cgloadct,65
      cgimage,rebin(chst_big,size_subim[1]*10,size_subim[2]*10,/sample),$
              /keep,pos=pos[*,i],opos=truepos,$
              noerase=i ne 0,missing_value=bad,missing_index=255,$
              missing_color='white'
      cgplot,[0],xsty=5,ysty=5,xran=[0,size_subim[1]],$
             yran=[0,size_subim[2]],position=truepos,$
             /nodata,/noerase,title=cap1,color='Black'
      cgoplot,[hst_big_ifsfov[*,0],hst_big_ifsfov[0,0]],$
              [hst_big_ifsfov[*,1],hst_big_ifsfov[0,1]],color='Black'
      imsize = string(fix(initmaps.hst.subim_big*kpc_per_as),format='(I0)')
      cgtext,size_subim[1]*0.05,size_subim[2]*0.9,$
             textoidl(imsize+'\times'+imsize+' kpc'),color='black'
      posbox1x[0] = truepos[0]+(truepos[2]-truepos[0])*$
                    hst_big_ifsfov[0,0]/size_subim[1]
      posbox1y[0] = truepos[1]+(truepos[3]-truepos[1])*$
                    hst_big_ifsfov[3,1]/size_subim[2]
      posbox2x[0] = posbox1x[0]
      posbox2y[0] = truepos[1]+(truepos[3]-truepos[1])*$
                    hst_big_ifsfov[0,1]/size_subim[2]


      i = 1
      size_subim = size(chst_fov)
      cgloadct,65
      cgimage,rebin(chst_fov,size_subim[1]*10,size_subim[2]*10,/sample),$
              /keep,pos=pos[*,i],opos=truepos,$
              noerase=i ne 0,missing_value=bad,missing_index=255,$
              missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
             /nodata,/noerase,title=cap2
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
      cgtext,-2,1.8,'IFS FOV',/data,color='black'
      ifsf_plotcompass,xarr_kpc,yarr_kpc
      zran=initmaps.hstcol.scllim
      dzran = zran[1]-zran[0]
      ncbdiv = initmaps.hstcol.ncbdiv
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
                         (dzran - zran[1]),format=cbform)
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
                 ticknames=ticknames,/ver,/right,charsize=0.6
      posbox1x[1] = truepos[0]
      posbox1y[1] = truepos[3]
      posbox2x[1] = truepos[0]
      posbox2y[1] = truepos[1]

      i = 2
      size_subim = size(cshst_fov_s)
;     smoothed HST continuum, IFS FOV
      if dohstcolsm then begin
         cgloadct,65
         cgimage,rebin(cshst_fov_s,size_subim[1]*10,size_subim[2]*10,/sample),$
                 /keep,pos=pos[*,i],opos=truepos,$
                 noerase=i ne 0,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase,title=cap3
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
         cgtext,-2,1.8,'IFS FOV, conv.',/data,color='black'
         zran=initmaps.hstcolsm.scllim
         dzran = zran[1]-zran[0]
         ncbdiv = initmaps.hstcolsm.ncbdiv
         ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
                            (dzran - zran[1]),format=cbform)
         cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
         cgcolorbar,position=cbpos,divisions=ncbdiv,$
                    ticknames=ticknames,/ver,/right,charsize=0.6
      endif
 
      i=3
 
      if tag_exist(initmaps,'ct') then begin
         ictlo = value_locate(ctcube.wave,initmaps.ct.sumrange[0])
         icthi = value_locate(ctcube.wave,initmaps.ct.sumrange[1])
         zran = initmaps.ct.scllim
         dzran = zran[1]-zran[0]
         ctmap = total(ctcube.dat[*,*,ictlo:icthi],3)
         ctmap /= max(ctmap)
         mapscl = cgimgscl(rebin(ctmap,dx*20,dy*20,/sample),$
                           minval=zran[0],max=zran[1],stretch=initmaps.ct.stretch)
         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos[*,i],opos=truepos,$
                 noerase=i ne 0,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase,title=cap4
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
      endif

      cgplot,posbox1x,posbox1y,color='Black',$
             xsty=5,ysty=5,/noerase,xran=[0,1],yran=[0,1],pos=[0,0,1,1]
      cgoplot,posbox2x,posbox2y,color='Black'

      cgps_close

   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Continuum radial profiles
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if tag_exist(initmaps,'ct') then begin
      
      npy = 2
      npx = 1
      if tag_exist(initdat,'decompose_qso_fit') then npx = 3

      cgps_open,initdat.mapdir+initdat.label+'cont_radprof.eps',charsize=1,/encap,$
                /inches,xs=plotquantum*2*npx,ys=plotquantum*2*npy,/qui

      pos = cglayout([npx,npy],ixmar=[2,2],iymar=[2,2],oxmar=[5,5],oymar=[10,5],$
                     xgap=5,ygap=5,unit=!D.X_PX_CM/3.0)
      
;     Total flux
      cgplot,map_rkpc_ifs,alog10(ctmap_save),yran=[-4,0],$
             xran=[0,max(map_rkpc_ifs)],/xsty,/ysty,psym=16,symsize=1d,$
             pos=pos[*,0],aspect=1d,title='Host Cont. + QSO PSF',$
             xtit = 'Radius (kpc)',ytit = 'log I/I$\downmax$
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

      mapscl = cgimgscl(rebin(ctmap_save,dx*20,dy*20,/sample),$
                        minval=zran[0],max=zran[1],stretch=initmaps.ct.stretch)
      cgloadct,65,/reverse
      if tag_exist(initdat,'decompose_qso_fit') then posuse=pos[*,3] else $
         posuse=pos[*,1]
      cgimage,mapscl,/keep,pos=posuse,opos=truepos,$
              /noerase,missing_value=bad,missing_index=255,$
              missing_color='white'
      if tag_exist(initmaps.ct,'fitifspeak') AND $
         tag_exist(initmaps.ct,'fitifspeakwin_kpc') then begin
          cgplot,[0],xsty=5,ysty=5,xran=[1,dx],$
                 yran=[1,dy],position=truepos,$
                 /nodata,/noerase
          cgoplot,peakfit_pix_ifs[0],peakfit_pix_ifs[1],psym=1,color='Red'         
      endif else begin
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase 
      endelse
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y

      if tag_exist(initdat,'decompose_qso_fit') then begin
         qso_map = total(contcube.qso,3)
         qso_map /= max(qso_map)
         cgplot,map_rkpc_ifs,alog10(qso_map),yran=[-4,0],$
                xran=[0,max(map_rkpc_ifs)],/xsty,/ysty,psym=16,symsize=1d,$
                pos=pos[*,1],/noerase,aspect=1d,title='QSO PSF'
         if tag_exist(initdat,'decompose_qso_fit') then begin
            cgoplot,psf1d_x,psf1d_y,color='Red'
         endif else if tag_exist(initmaps,'fit_empsf') then begin
            cgoplot,empsf1d_x,empsf1d_y,color='Red'
         endif else if tag_exist(initmaps,'ctradprof_psffwhm') then begin
            x = dindgen(101)/100d*max(map_rkpc_ifs)
            fwhm=initmaps.ctradprof_psffwhm * kpc_per_as
;           Gaussian
            y = alog10(gaussian(x,[1d,0d,fwhm/2.35]))
            cgoplot,x,y,color='Black'
;           Moffat, index = 1.5
            y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/1.5d)-1),1.5d]))
            cgoplot,x,y,color='Red',/linesty
;           Moffat, index = 2.5
            y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/2.5d)-1),2.5d]))
            cgoplot,x,y,color='Red'
;           Moffat, index = 5
            y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/5d)-1),5d]))
            cgoplot,x,y,color='Blue'
         endif

         mapscl = cgimgscl(rebin(qso_map,dx*20,dy*20,/sample),$
                           minval=zran[0],max=zran[1],stretch=initmaps.ct.stretch)
         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos[*,4],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         if tag_exist(initmaps.ct,'fitifspeak') AND $
            tag_exist(initmaps.ct,'fitifspeakwin_kpc') then begin
            peakfit_pix = [qso_fitpar[4]+1,qso_fitpar[5]+1]
            peakfit_ifs_qso_distance_from_nucleus_pix = peakfit_pix - $
                                                        center_nuclei[*,0]
            peakfit_ifs_qso_distance_from_nucleus_kpc = $
                peakfit_ifs_qso_distance_from_nucleus_pix * kpc_per_pix
            cgplot,[0],xsty=5,ysty=5,xran=[1,dx],$
                   yran=[1,dy],position=truepos,$
                   /nodata,/noerase
            cgoplot,peakfit_pix[0],peakfit_pix[1],psym=1,color='Red'         
         endif else begin
            cgplot,[0],xsty=5,ysty=5,position=truepos,$
                   /nodata,/noerase 
         endelse
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y

         host_map = total(contcube.host,3)
         host_map /= max(host_map)
         cgplot,map_rkpc_ifs,alog10(host_map),yran=[-4,0],$
                xran=[0,max(map_rkpc_ifs)],/xsty,/ysty,psym=16,symsize=1d,$
                pos=pos[*,2],/noerase,aspect=1d,title='Host Continuum'
         if tag_exist(initdat,'decompose_qso_fit') then begin
            cgoplot,psf1d_x,psf1d_y,color='Red'
         endif else if tag_exist(initmaps,'fit_empsf') then begin
            cgoplot,empsf1d_x,empsf1d_y,color='Red'
         endif else if tag_exist(initmaps,'ctradprof_psffwhm') then begin
            x = dindgen(101)/100d*max(map_rkpc_ifs)
            fwhm=initmaps.ctradprof_psffwhm * kpc_per_as
;           Gaussian
            y = alog10(gaussian(x,[1d,0d,fwhm/2.35]))
            cgoplot,x,y,color='Black'
;           Moffat, index = 1.5
            y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/1.5d)-1),1.5d]))
            cgoplot,x,y,color='Red',/linesty
;           Moffat, index = 2.5
            y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/2.5d)-1),2.5d]))
            cgoplot,x,y,color='Red'
;           Moffat, index = 5
            y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/5d)-1),5d]))
            cgoplot,x,y,color='Blue'
         endif
         
         mapscl = cgimgscl(rebin(host_map,dx*20,dy*20,/sample),$
                           minval=zran[0],max=zran[1],stretch=initmaps.ct.stretch)
         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos[*,5],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
      endif

      cgps_close
   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Continuum components: Percentage of summed continuum fit by additive polynomial
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;   cgps_open,initdat.mapdir+initdat.label+'cont_polypct.eps',charsize=1,/encap,$
;         /inches,xs=plotquantum*2,ys=plotquantum*2,/qui
;
;   map = contcube.cont_fit_poly_tot_pct
;   ibd = where(map eq bad)
;   map[ibd] = 0l
;      
;   zran = [0,1]
;   dzran = zran[1]-zran[0]
;   mapscl = cgimgscl(rebin(map,dx*20,dy*20,/sample),$
;                     minval=zran[0],maxval=zran[1],ncolors=5)
;   cgloadct,8,/brewer,ncolors=5
;   cgimage,mapscl,/keep,opos=truepos,mar=0.2d
;   cgplot,[0],xsty=5,ysty=5,position=truepos,$
;          /nodata,/noerase
;   ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;   cbpos=[truepos[2],truepos[1],truepos[2]+0.02,truepos[3]]
;   ticknames = ['0.0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1.0']
;   cgDCBar,ncolors=5,position=cbpos,labels=ticknames,/ver,/right,$
;           charsize=0.6
;   cgps_close
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Continuum components: Percentage of summed continuum fit by additive polynomial
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;   cgps_open,initdat.mapdir+initdat.label+'NaDcont_polypct.eps',charsize=1,/encap,$
;         /inches,xs=plotquantum*2,ys=plotquantum*2,/qui
;
;   map = contcube.cont_fit_poly_nad_pct
;   ibd = where(map eq bad)
;   map[ibd] = 0l
;      
;   zran = [0,1]
;   dzran = zran[1]-zran[0]
;   mapscl = cgimgscl(rebin(map,dx*20,dy*20,/sample),$
;                     minval=zran[0],maxval=zran[1],ncolors=5)
;   cgloadct,8,/brewer,ncolors=5
;   cgimage,mapscl,/keep,opos=truepos,mar=0.2d
;   cgplot,[0],xsty=5,ysty=5,position=truepos,$
;          /nodata,/noerase
;   ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;   cbpos=[truepos[2],truepos[1],truepos[2]+0.02,truepos[3]]
;   ticknames = ['0.0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1.0']
;   cgDCBar,ncolors=5,position=cbpos,labels=ticknames,/ver,/right,$
;           charsize=0.6
;   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plots of individual emission lines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if ~ tag_exist(initdat,'noemlinfit') then begin

;     Size of plot grid
      npx = initdat.maxncomp
      npy = 3

;     Linemap indices to plot below: flux, wavelength (converted to velocity), sigma
      ilinmap = [0,2,3]

;     Loop through emission lines
      foreach line,outlines do begin

;        Get syntax of linelabel right; otherwise call to DEVICE chokes
         linelab=line
         ilb = strpos(linelab,'[')
         if ilb ne -1 then $
            linelab = strmid(linelab,0,ilb)+'\'+strmid(linelab,ilb)
         irb = strpos(linelab,']')
         if irb ne -1 then $
            linelab = strmid(linelab,0,irb)+'\'+strmid(linelab,irb)

         cgps_open,initdat.mapdir+initdat.label+linelab+'_comp.eps',charsize=1,/encap,$
                   /inches,xs=plotquantum*npx,ys=plotquantum*npy*aspectrat,/qui
         pos = cglayout([npx,npy],ixmar=[3,3],iymar=[3,3],oxmar=[0,0],oymar=[0,0],$
                        xgap=0,ygap=0,unit=!D.X_PX_CM/3.0)

;        loop through plot types
         for j=0,2 do begin

;           Set up colorbar labeling
            if j eq 0 then cbform = '(D0.1)' else cbform = '(I0)'


;           Set up ranges for all components at once
            if tag_exist(initmaps,'setplotrange_allcomp') then begin
        
               mapallcomp = linmaps[line,*,*,0:initdat.maxncomp-1,ilinmap[j]]
               ibd = where(mapallcomp eq bad,ctbd)
               igd = where(mapallcomp ne bad AND mapallcomp ne 0,ctgd)
               if j eq 1 then begin
;                 redshift with respect to galaxy systemic
                  zdiff = mapallcomp[igd]/linelist[line]-1d - initdat.zsys_gas
;                 relativistic velocity shift;
;                 see http://hyperphysics.phy-astr.gsu.edu/hbase/relativ/reldop2.html
                  mapallcomp[igd] = c_kms * ((zdiff+1d)^2d - 1d) / ((zdiff+1d)^2d + 1d)
               endif
               zran = [min(mapallcomp[igd]),max(mapallcomp[igd])]
               if j eq 0 then begin
                  zmax_flux = zran[1]
                  zran=[0,1]
                  dzran = 1
                  ncbdiv = 5
               endif else begin
                  divarr = ifsf_cbdiv(zran,100d,ncbdivmax)
                  ncbdiv = divarr[0]
                  dzran = zran[1]-zran[0]
               endelse

            endif

;           Loop through velocity components
            for i=0,initdat.maxncomp-1 do begin

;              Plot index
               iplot = j*npx+i

;              Get map and scale
               map = linmaps[line,*,*,i,ilinmap[j]]
               fluxmap = linmaps[line,*,*,i,0]
               ibd = where(map eq bad AND ~ finite(map),ctbd)
               inan = where(~finite(map),ctnan)
               izero = where(fluxmap eq 0)
               igd = where(map ne bad AND fluxmap ne 0 AND finite(map),ctgd)

               map[ibd] = bad
               map[inan] = bad
               map[izero] = bad

               if ctgd gt 0 then begin
            
                  if j eq 0 AND tag_exist(initmaps,'fluxfactor') then begin
                     map[igd] *= initmaps.fluxfactor
                  endif else if j eq 1 then begin
                     zdiff = map[igd]/linelist[line]-1d - initdat.zsys_gas
                     map[igd] = c_kms * ((zdiff+1d)^2d - 1d) / ((zdiff+1d)^2d + 1d)
                  endif

;                 Set up ranges for each component separately

;                 Check for manual range first ...
                  hasrange = 0
                  if hasrangefile then begin
                     ithisline = where(line eq rangeline AND $
                                       i+1 eq rangecomp AND $
                                       plottypes[j] eq rangequant,ctthisline)
                     if ctthisline eq 1 then begin
                        zran = [rangelo[ithisline],rangehi[ithisline]]
                        dzran = zran[1]-zran[0]
                        ncbdiv = rangencbdiv[ithisline]
                        ncbdiv = ncbdiv[0]
                        hasrange = 1
                     endif
                  endif
;                 otherwise set it automagically.
                  if ~tag_exist(initmaps,'setplotrange_allcomp') AND $
                     ~hasrange then begin
                     zran = [min(map[igd]),max(map[igd])]
                     if j eq 0 then begin
                        zmax_flux = zran[1]
                        zran=[0,1]
                        dzran = 1
                        ncbdiv = 5
                     endif else begin
                        divarr = ifsf_cbdiv(zran,100d,ncbdivmax)
                        ncbdiv = divarr[0]
                        dzran = zran[1]-zran[0]
                     endelse
                  endif
 
                  if j eq 0 then map[igd] = map[igd]/zmax_flux
                  if ctnan gt 0 then map[inan] = bad
                  mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
                                  min=zran[0],max=zran[1])
   
;                 Plot image
                  if j eq 0 then begin
                     cgloadct,65,/reverse
                     title='c'+string(i+1,format='(I0)')+' flux'
                     title += ' ('+string(zmax_flux,format='(E0.2)')+')'
                  endif
                  if j eq 1 then begin
                     cgloadct,74,/reverse
                     title='c'+string(i+1,format='(I0)')+' velocity'
                  endif
                  if j eq 2 then begin
                     cgloadct,65,/reverse
                     title='c'+string(i+1,format='(I0)')+' sigma'
                  endif
                  cgimage,mapscl,/keep,pos=pos[*,iplot],opos=truepos,$
                          noerase=iplot ne 0,missing_value=bad,missing_index=255,$
                          missing_color='white'
                  cgplot,[0],xsty=5,ysty=5,position=truepos,$
                         /nodata,/noerase,title=title
                  ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;                 Colorbar
                  cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
                  ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
                                     (dzran - zran[1]),format=cbform)
                  cgcolorbar,position=cbpos,divisions=ncbdiv,$
                             ticknames=ticknames,/ver,/right,charsize=0.6

               endif
 
            endfor

         endfor
 
         cgps_close

      endforeach

   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; 1D velocity cross-sections [CVDF]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  hash to hold endpoints for plotting cross-sections on maps
   xsec_endpoints = 0b
   if tag_exist(initmaps,'xsec') then begin

      xsec = initmaps.xsec
      nsec = n_elements(xsec.line)
      xsec_endpoints = dblarr(2,2,nsec)

      for i=0,nsec-1 do begin
;        Get map and scale
         itag = where(strcmp(linspecpars_tags,xsec.tag[i],/fold_case) eq 1)
         map = linspecpars[xsec.line[i]].(itag)
         ibd = where(map eq bad OR ~finite(map),ctbd)
         igd = where(map ne bad AND finite(map),ctgd)

         if ctgd gt 0 then begin
            hasrange = 0
            if hasrangefile then begin
            ithisline = where(xsec.line[i] eq rangeline AND $
                              xsec.tag[i] eq rangequant,ctthisline)
               if ctthisline eq 1 then begin
                  zran = [rangelo[ithisline],rangehi[ithisline]]
                  dzran = zran[1]-zran[0]
                  ncbdiv = rangencbdiv[ithisline]                   
                  ncbdiv = ncbdiv[0]
                  hasrange = 1
               endif
            endif
         endif
         if ~hasrange then begin
            zran = [min(map[igd]),max(map[igd])]
            divarr = ifsf_cbdiv(zran,100d,ncbdivmax)
            ncbdiv = divarr[0]
            dzran = zran[1]-zran[0]
         endif
 
         if ctbd gt 0 then map[ibd] = bad

         if tag_exist(xsec,'nearest') then nearest=1b else nearest=0b
         slit = ifsfit_xsec(map,[xsec.xcenter[i],xsec.ycenter[i]],$
                            xsec.length[i],xsec.angle[i],bad=bad,$
                            nearest=nearest,ends=ends)
         xsec_endpoints[*,*,i] = ends
         islitgd = where(slit[*,1] ne bad,ctslitgd)
         if ctslitgd gt 0 then begin
            rslit = slit[islitgd,0]
            mslit = slit[islitgd,1]
         endif else begin
            rslit=0d
            mslit=0d
         endelse

         file = initdat.mapdir+initdat.label+xsec.line[i]+'_cvdf_xsec_'+$
                xsec.tag[i]+'_'+string(xsec.angle[i],format='(I0)')+'deg.eps'
         cgps_open,file,charsize=1,/encap,/inches,xs=7.5d,ys=7.5d,/qui
         xran = 1.05d*[min(rslit),max(rslit)]*kpc_per_pix
         cgplot,rslit*kpc_per_pix,mslit,/xsty,/ysty,$
                xran=xran,yran=zran,xtit='Position (kpc)',$
                ytit=xsec.title[i],psym=16,symsize=2
;         cgoplot,rslit*kpc_per_as,mslit
         cgoplot,[0d,0d],zran,/linesty
         cgoplot,xran,[0d,0d],/linesty
         cgps_close

      endfor

   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plots of individual emission lines [CVDF]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if ~ tag_exist(initdat,'noemlinfit') then begin

;     Quantities to plot
      if tag_exist(initmaps,'cvdf_ftags') then ftags = initmaps.cvdf_ftags $
      else ftags = ['ftot','fpk','fv50','fv84','fv98']
      if tag_exist(initmaps,'cvdf_vtags') then vtags = initmaps.cvdf_vtags $
      else vtags = ['sig','vpk','v50','v84','v98']
      if tag_exist(initmaps,'cvdf_ftitles') then ftitles = initmaps.cvdf_ftitles $
      else ftitles = ['F$\downtot$','F$\downpk$','F$\downv50$',$
                    'F$\downv84$','F$\downv98$']
      if tag_exist(initmaps,'cvdf_vtitles') then vtitles = initmaps.cvdf_vtitles $
      else vtitles = ['$\sigma$','v$\downpeak$','v$\down50$',$
                      'v$\down84$','v$\down98$']

;     Size of plot grid
      npx = max([n_elements(vtags),n_elements(ftags)])
      if tag_exist(initmaps,'doemlinradprof') then npy = 3 else npy = 2
      
;     Loop through emission lines
      foreach line,outlines do begin

;        Get syntax of linelabel right; otherwise call to DEVICE chokes
         linelab=line
         ilb = strpos(linelab,'[')
         if ilb ne -1 then $
            linelab = strmid(linelab,0,ilb)+'\'+strmid(linelab,ilb)
         irb = strpos(linelab,']')
         if irb ne -1 then $
            linelab = strmid(linelab,0,irb)+'\'+strmid(linelab,irb)

         cgps_open,initdat.mapdir+initdat.label+linelab+'_cvdf.eps',$
                   charsize=1,/encap,$
                   /inches,xs=plotquantum*npx,ys=plotquantum*npy*aspectrat,/qui
         pos = cglayout([npx,npy],ixmar=[3,3],iymar=[3,3],$
                        oxmar=[0,0],oymar=[0,0],$
                        xgap=0,ygap=0,unit=!D.X_PX_CM/3.0)

         ichan=0
;        loop through plot types
         for j=0,npx-1 do begin

;           FLUXES

            cbform = '(D0.1)' ; colorbar syntax
            iplot = j ; plot index
            if j lt n_elements(ftags) then begin

;              Get map and scale
;              Channel map option
               ftag_use = ftags[j]
               if ftags[j] eq 'ch' then begin
                  if tag_exist(initmaps,'cvdf_channels') then begin
                     ftag_use = 'ch'+string(ichan,format='(I0)')
                     linmap_tmp = linspecmaps[line]
                     vel = initmaps.cvdf_channels[ichan]
                     ivel = $
                        value_locate(linmap_tmp.vel,vel)
                     map = linmap_tmp.flux[*,*,ivel]
                     ichan++
                  endif else begin
                     map = 0b
                  endelse
;              vel% map option
               endif else begin
                  itag = where(strcmp(linspecpars_tags,ftags[j],/fold_case) eq 1)
                  map = linspecpars[line].(itag)
               endelse
               ibd = where(map eq bad AND ~ finite(map),ctbd)
               inan = where(~finite(map),ctnan)
               igd = where(map ne bad AND map ne 0 AND finite(map),ctgd)
               ctgdlin = ctgd ; saving for later ...
            
               if ctgd gt 0 then begin
               
                  if tag_exist(initmaps,'fluxfactor') then $
                     map[igd] *= initmaps.fluxfactor
                                 
                  zran=[0,1]
                  dzran = 1
                  ncbdiv = 5
                  zmax_flux = max(map[igd])

                  if hasrangefile then begin
                     ithisline = where(line eq rangeline AND $
                                       ftag_use eq rangequant,ctthisline)
                     if ctthisline eq 1 then zmax_flux = rangelo[ithisline]
                  endif
               
                  map[igd] = map[igd]/zmax_flux[0]
                  if ctnan gt 0 then map[inan] = bad
                  izero = where(map eq 0d,ctzero)
                  if ctzero gt 0 then map[izero] = bad
                  mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
                                  min=zran[0],max=zran[1])
   
;                 Save for radial profiles
                  maplin = map
                  maplin[igd] = alog10(maplin[igd])
                  if ctzero gt 0 then maplin[izero] = bad
   
;                 Plot image
                  cgloadct,65,/reverse
                  title=ftitles[j]
                  title += ' ('+string(zmax_flux,format='(E0.2)')+')'
                  cgimage,mapscl,/keep,pos=pos[*,iplot],opos=truepos,$
                          noerase=iplot ne 0,missing_value=bad,missing_index=255,$
                          missing_color='white'
                  cgplot,[0],xsty=5,ysty=5,position=truepos,$
                         /nodata,/noerase,title=title
                  ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,$
                                   center_nuclei_kpc_y
;                 Colorbar
                  cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
                  ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
                                     (dzran - zran[1]),format=cbform)
                  cgcolorbar,position=cbpos,divisions=ncbdiv,$
                             ticknames=ticknames,/ver,/right,charsize=0.6
               endif
            endif else ctgdlin=0

;           VELOCITIES
                          
            cbform = '(I0)' ; colorbar syntax

            if j lt n_elements(vtags) then begin


               iplot = npx+j ; plot index

;              Get map and scale
               itag = where(strcmp(linspecpars_tags,vtags[j],/fold_case) eq 1)
               map = linspecpars[line].(itag)
               ibd = where(map eq bad AND ~ finite(map),ctbd)
               inan = where(~finite(map),ctnan)
               igd = where(map ne bad AND map ne 0 AND finite(map),ctgd)

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
 
                  if ctnan gt 0 then map[inan] = bad
                  mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
                                  min=zran[0],max=zran[1])
   
;                 Plot image
                  if vtags[j] eq 'sig' then cgloadct,65,/reverse $
                  else cgloadct,74,/reverse
                  title=vtitles[j]
                  cgimage,mapscl,/keep,pos=pos[*,iplot],opos=truepos,$
                          noerase=iplot ne 0,missing_value=bad,missing_index=255,$
                          missing_color='white'
                  cgplot,[0],xsty=5,ysty=5,position=truepos,$
                         /nodata,/noerase,title=title,xran=[0,dx],yran=[0,dy]
;                 Velocity contours
                  if tag_exist(initmaps,'contourlevels') then begin
                     key = line+'_'+vtags[j]
;                    Not sure why levels aren't being labeled
                     if initmaps.contourlevels->haskey(key) then begin
;                        nlevels = n_elements(initmaps.contourlevels[key])
                        cgcontour,map,dindgen(dx)+0.5,dindgen(dy)+0.5,$
                                  /overplot,color=0,c_linesty=2,c_thick=4,$
                                  levels=initmaps.contourlevels[key],$
                                  max=initmaps.contourmax[key] ;,$
;                                  c_labels=dblarr(nlevels)+1b,$
;                                  c_charsize=0.75
                     endif
                  endif
;                 Cross section
                  if xsec_endpoints[0] ne 0b then $
                     for k=0,n_elements(initmaps.xsec.line)-1 do $
                        if initmaps.xsec.line[k] eq line AND $
                           initmaps.xsec.tag[k] eq vtags[j] then $
                           cgoplot,xsec_endpoints[*,0,k]+0.5d,$
                                   xsec_endpoints[*,1,k]+0.5d,$
                                   thick=4,linesty=0
                  ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,$
                                   center_nuclei_kpc_y
;                 Colorbar
                  cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
                  ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
                                     (dzran - zran[1]),format=cbform)
                  cgcolorbar,position=cbpos,divisions=ncbdiv,$
                             ticknames=ticknames,/ver,/right,charsize=0.6

               endif               
            endif

;           Radial profiles

            if ctgdlin gt 0 AND tag_exist(initmaps,'doemlinradprof') then begin
;
; Moffat index values chosen to match turbulence theory (5), IRAF default (2.5),
; and wingy profile (1.5). These are the same chosen by Trujillo et al. 2001.
                          
               iplot = 2*npx+j ; plot index
               yran = [-4d,0d]

               title = ftitles[j]+' vs. R'
               cgplot,map_rkpc_ifs,maplin,/xsty,/ysty,yran=yran,$
                      xran=[0,max(map_rkpc_ifs)],/noerase,$
                      title=title,pos=pos[*,iplot],psym=16,symsize=1d,$
                      aspect=1.0d
               if tag_exist(initdat,'decompose_qso_fit') then begin
                  cgoplot,psf1d_x,psf1d_y,color='Red'
               endif else if tag_exist(initmaps,'fit_empsf') then begin
                  cgoplot,empsf1d_x,empsf1d_y,color='Red'
               endif else if tag_exist(initmaps,'emlinradprof_psffwhm') then begin
                  x = dindgen(101)/100d*max(map_rkpc_ifs)
                  fwhm=initmaps.emlinradprof_psffwhm * kpc_per_as
;                 Gaussian
                  y = alog10(gaussian(x,[1d,0d,fwhm/2.35]))
                  cgoplot,x,y,color='Black'
;                 Moffat, index = 1.5
                  y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/1.5d)-1),1.5d]))
                  cgoplot,x,y,color='Red',/linesty
;                 Moffat, index = 2.5
                  y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/2.5d)-1),2.5d]))
                  cgoplot,x,y,color='Red'
;                 Moffat, index = 5
                  y = alog10(moffat(x,[1d,0d,fwhm/2d/sqrt(2^(1/5d)-1),5d]))
                  cgoplot,x,y,color='Blue'
               endif

            endif
 
         endfor
 
         cgps_close

      endforeach

   endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plots of line ratios
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if hasratpar_comp then begin

      keys = comp_lr_args->keys()

;     Loop through plot files
      foreach plt,keys do begin

;        Size of plot grid
         arrsiz = size(comp_lr_args[plt])
         npx = arrsiz[1]
         if arrsiz[0] gt 1 then begin
            npy = arrsiz[2]
            nplots = npx*npy
         endif else begin
            npy = 1
            nplots = npx
         endelse

         cgps_open,initdat.mapdir+initdat.label+plt+'_comp.eps',charsize=1,/encap,$
            /inches,xs=plotquantum*npx,ys=plotquantum*npy*aspectrat,/qui
         pos = cglayout([npx,npy],ixmar=[3,3],iymar=[3,3],oxmar=[0,0],oymar=[0,0],$
            xgap=0,ygap=0,unit=!D.X_PX_CM/3.0)
         cbform = '(D0.2)'

;        Loop through plot panels
         for i=0,nplots-1 do begin

            tmpstr = strsplit(comp_lr_args[plt,i],'_',/extract)
            vcomp = fix(tmpstr[0])
            ptype = tmpstr[1]

            map = linrats[ptype[0],*,*,vcomp-1]
            ibd = where(map eq bad AND ~ finite(map),ctbd)
            inan = where(~finite(map),ctnan)
            igd = where(map ne bad AND finite(map),ctgd)

;           Line ratio maps
            if n_elements(tmpstr) eq 2 AND ctgd gt 0 then begin
               
               hasrange = 0
               if hasrangefile then begin
                  ithisline = where(ptype eq rangeline AND $
                     vcomp eq rangecomp,ctthisline)
                  if ctthisline eq 1 then begin
                     zran = [rangelo[ithisline],rangehi[ithisline]]
                     dzran = zran[1]-zran[0]
                     ncbdiv = rangencbdiv[ithisline]
                     ncbdiv = ncbdiv[0]
                     hasrange = 1
                  endif
               endif
               if ~tag_exist(initmaps,'setplotrange_allcomp') AND $
                  ~hasrange then begin
                  zran = [min(map[igd]),max(map[igd])]
                  dzran = zran[1] - zran[0]
                  ncbdiv = ifsf_cbdiv(zran,0.5,7)
               endif

               title='c'+string(vcomp,format='(I0)')+' '
               if ptype eq 'n2ha' then title+=textoidl('[NII]/H\alpha')
               if ptype eq 'o1ha' then title+=textoidl('[OI]/H\alpha')
               if ptype eq 'o3hb' then title+=textoidl('[OIII]/H\beta')
               if ptype eq 'ebv' then title+=textoidl('E(B-V)')

               if ctnan gt 0 then map[inan] = bad
               mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
                               min=zran[0],max=zran[1])

               cgloadct,65,/reverse
               cgimage,mapscl,/keep,pos=pos[*,i],opos=truepos,$
                       noerase=i ne 0,missing_value=bad,missing_index=255,$
                       missing_color='white'
               cgplot,[0],xsty=5,ysty=5,position=truepos,$
                      /nodata,/noerase,title=title
               ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
               cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
               ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
                                  (dzran - zran[1]),format=cbform)
               cgcolorbar,position=cbpos,divisions=ncbdiv,$
                  ticknames=ticknames,/ver,/right,charsize=0.6

;           VO plots
            endif else if n_elements(tmpstr) eq 4 AND ctgd gt 0 then begin
               
               map2 = linrats[tmpstr[3],*,*,vcomp-1]
               ibd2 = where(map2 eq bad AND ~ finite(map2),ctbd2)
               inan2 = where(~finite(map2),ctnan2)
               igd2 = where(map2 ne bad AND finite(map2),ctgd2)


               if ctgd2 gt 0 then begin

                  title='c'+string(vcomp,format='(I0)')
                  ptype = ptype+'_vs_'+tmpstr[3]
                  if ptype eq 'n2ha_vs_o3hb' then begin
                     xran = [-1.99d,0.99d]
                     yran = [-1.19d,1.49d]
                     xkew1 = 0.05d*indgen(110)-5d
                     ykew1 = 0.61d / (xkew1-0.47d)+1.19d
                     xkew2 = xkew1
                     ykew2 = xkew1-xkew1-99d
                     xtit = textoidl('[NII]/H\alpha')
                     ytit = textoidl('[OIII]/H\beta')
                  endif
                  if ptype eq 'o1ha_vs_o3hb' then begin
                     xran = [-1.99d,0.49d]
                     yran = [-1.19d,1.49d]
                     xkew1 = 0.05d*indgen(110)-5d
                     ykew1 = 0.73d / (xkew1+0.59d)+1.33d
                     xkew2 = xkew1
                     ykew2 = xkew1-xkew1-99d
                     xtit = textoidl('[OI]/H\alpha')
                     ytit = textoidl('[OIII]/H\beta')
                  endif

                  if ctnan gt 0 then map[inan] = bad
                  if ctnan2 gt 0 then map2[inan2] = bad

                  cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,pos=pos[*,i],$
                     /nodata,noerase=i ne 0,title=title,xtit=xtit,ytit=ytit,$
                     aspect=1d
                  cgoplot,map,map2,psym=1
                  cgoplot,xkew1,ykew1
                  cgoplot,xkew2,ykew2
                  
               endif
               
            endif

         endfor

         cgps_close

      endforeach

   endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plots of line ratios (by CVDF)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if hasratpar_cvdf then begin

      ftags = cvdf_lr_ftags
      titles = cvdf_lr_ftitles
      keys = cvdf_lr_args->keys()

;     Loop through plot files
      foreach plt,keys do begin

;        Size of plot grid
         arrsiz = size(cvdf_lr_args[plt])
         npx = arrsiz[1]
         if arrsiz[0] gt 1 then begin
            npy = arrsiz[2]
            nplots = npx*npy
         endif else begin
            npy = 1
            nplots = npx
         endelse

         cgps_open,initdat.mapdir+initdat.label+plt+'_cvdf.eps',charsize=1,/encap,$
            /inches,xs=plotquantum*npx,ys=plotquantum*npy*aspectrat,/qui
         pos = cglayout([npx,npy],ixmar=[3,3],iymar=[3,3],oxmar=[0,0],oymar=[0,0],$
            xgap=0,ygap=0,unit=!D.X_PX_CM/3.0)
         cbform = '(D0.2)'

;        Loop through plot panels
         for i=0,nplots-1 do begin

            tmpstr = strsplit(cvdf_lr_args[plt,i],'_',/extract)
            ftype = tmpstr[0]
            vcomp = where(ftype eq ftags)
            ptype = tmpstr[1]

            map = linrats_cvdf[ptype[0],*,*,vcomp]
            ibd = where(map eq bad AND ~ finite(map),ctbd)
            inan = where(~finite(map),ctnan)
            igd = where(map ne bad AND finite(map),ctgd)

;           Line ratio maps
            if n_elements(tmpstr) eq 2 AND ctgd gt 0 then begin
               
               hasrange = 0
               if hasrangefile then begin
                  ithisline = where(ptype eq rangeline AND $
                     ftype eq rangequant,ctthisline)
                  if ctthisline eq 1 then begin
                     zran = [rangelo[ithisline],rangehi[ithisline]]
                     dzran = zran[1]-zran[0]
                     ncbdiv = rangencbdiv[ithisline]
                     ncbdiv = ncbdiv[0]
                     hasrange = 1
                  endif
               endif
               if ~tag_exist(initmaps,'setplotrange_allcomp') AND $
                  ~hasrange then begin
                  zran = [min(map[igd]),max(map[igd])]
                  dzran = zran[1] - zran[0]
                  ncbdiv = ifsf_cbdiv(zran,0.5,7)
               endif

               title = titles[vcomp]
               if ptype eq 'n2ha' then title+=textoidl('[NII]/H\alpha')
               if ptype eq 'o1ha' then title+=textoidl('[OI]/H\alpha')
               if ptype eq 'o3hb' then title+=textoidl('[OIII]/H\beta')
               if ptype eq 'ebv' then title+=textoidl('E(B-V)')

               if ctnan gt 0 then map[inan] = bad
               mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
                               min=zran[0],max=zran[1])

               cgloadct,65,/reverse
               cgimage,mapscl,/keep,pos=pos[*,i],opos=truepos,$
                       noerase=i ne 0,missing_value=bad,missing_index=255,$
                       missing_color='white'
               cgplot,[0],xsty=5,ysty=5,position=truepos,$
                      /nodata,/noerase,title=title
               ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
               cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
               ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
                                  (dzran - zran[1]),format=cbform)
               cgcolorbar,position=cbpos,divisions=ncbdiv,$
                  ticknames=ticknames,/ver,/right,charsize=0.6

;           VO plots
            endif else if n_elements(tmpstr) eq 4 AND ctgd gt 0 then begin
               
               map2 = linrats_cvdf[tmpstr[3],*,*,vcomp]
               ibd2 = where(map2 eq bad AND ~ finite(map2),ctbd2)
               inan2 = where(~finite(map2),ctnan2)
               igd2 = where(map2 ne bad AND finite(map2),ctgd2)


               if ctgd2 gt 0 then begin
                  
                  title = titles[vcomp]
                     
                  ptype = ptype+'_vs_'+tmpstr[3]
                  if ptype eq 'n2ha_vs_o3hb' then begin
                     xran = [-1.99d,0.99d]
                     yran = [-1.19d,1.49d]
                     xkew1 = 0.05d*indgen(110)-5d
                     ykew1 = 0.61d / (xkew1-0.47d)+1.19d
                     xkew2 = xkew1
                     ykew2 = xkew1-xkew1-99d
                     xtit = textoidl('[NII]/H\alpha')
                     ytit = textoidl('[OIII]/H\beta')
                  endif
                  if ptype eq 'o1ha_vs_o3hb' then begin
                     xran = [-1.99d,0.49d]
                     yran = [-1.19d,1.49d]
                     xkew1 = 0.05*indgen(85)-5
                     ykew1 = 0.73d / (xkew1+0.59d)+1.33d
                     xkew2 = xkew1
                     ykew2 = xkew1-xkew1-99d
                     xtit = textoidl('[OI]/H\alpha')
                     ytit = textoidl('[OIII]/H\beta')
                  endif

                  if ctnan gt 0 then map[inan] = bad
                  if ctnan2 gt 0 then map2[inan2] = bad

                  cgplot,[0],/xsty,/ysty,xran=xran,yran=yran,pos=pos[*,i],$
                     /nodata,noerase=i ne 0,title=title,xtit=xtit,ytit=ytit,$
                     aspect=1d
                  cgoplot,map,map2,psym=1
                  cgoplot,xkew1,ykew1
                  cgoplot,xkew2,ykew2
                  
               endif
               
            endif

         endfor

         cgps_close

      endforeach

   endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD EQUIVALENT WIDTH
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if tag_exist(initdat,'donad') then begin

      if maxnademncomp_act gt 0 then ny=3 else ny=1

      cgps_open,initdat.mapdir+initdat.label+'NaDempweq.eps',charsize=1,/encap,$
         /inches,xs=plotquantum*2,ys=plotquantum*ny,/qui

      pos = cglayout([2,ny],ixmar=[2,2],iymar=[2,2],oxmar=[0,0],oymar=[0,0],$
         xgap=0,ygap=0,unit=!D.X_PX_CM/3.0)
      cbform = '(I0)'

;
;     ABSORPTION
;
      map = nadcube.weq[*,*,0]
      igd = where(map ge initmaps.nadabsweq_snrthresh*nadcube.weq[*,*,1] AND $
                  map gt 0d AND map ne bad)
      ibd = where(map lt initmaps.nadabsweq_snrthresh*nadcube.weq[*,*,1] OR $
                  map eq 0d OR map eq bad)

;     Set up range
;     Check for manual range first ...
      hasrange = 0
      if hasrangefile then begin
         ithisline = where(rangeline eq 'NaDabs' AND $
                           rangequant eq 'weq',ctthisline)
         if ctthisline eq 1 then begin
            zran = [rangelo[ithisline],rangehi[ithisline]]
            dzran = zran[1]-zran[0]
            ncbdiv = rangencbdiv[ithisline]
            ncbdiv = ncbdiv[0]
            hasrange = 1
         endif
      endif
;     otherwise set it automagically.
      if ~hasrange then begin
         zran = [min(map[igd]),max(map[igd])]
         divarr = ifsf_cbdiv(zran,2d,ncbdivmax)
         ncbdiv = divarr[0]
         dzran = zran[1]-zran[0]
      endif

;     replace bad points
      map[ibd] = bad

;     Save some things for later
      zran_empweqabs = zran
      map_empweqabs = map

      mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
                      min=zran[0],max=zran[1])
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos[*,0],opos=truepos,$
         noerase=i ne 0,missing_value=bad,missing_index=255,$
         missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
;         /nodata,/noerase,title='W$\down eq$(NaD abs, $\angstrom$)'
         /nodata,/noerase,title=textoidl('W_{eq}^{abs}')+' ($\angstrom$)'
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6

;     Error
      err = nadcube.weq[*,*,1]
      map[igd] /= err[igd]

;     Set up range
;     Check for manual range first ...
      hasrange = 0
      if hasrangefile then begin
         ithisline = where(rangeline eq 'NaDabs' AND $
            rangequant eq 'weqsnr',ctthisline)
         if ctthisline eq 1 then begin
            zran = [rangelo[ithisline],rangehi[ithisline]]
            dzran = zran[1]-zran[0]
            ncbdiv = rangencbdiv[ithisline]
            ncbdiv = ncbdiv[0]
            hasrange = 1
         endif
      endif
;     otherwise set it automagically.
      if ~hasrange then begin
         zran = [min(map[igd]),max(map[igd])]
         divarr = ifsf_cbdiv(zran,10d,ncbdivmax)
         ncbdiv = divarr[0]
         dzran = zran[1]-zran[0]
      endif

      mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
         min=zran[0],max=zran[1])
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos[*,1],opos=truepos,$
         noerase=i ne 0,missing_value=bad,missing_index=255,$
         missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
         /nodata,/noerase,title=textoidl('W_{eq}^{abs}/\deltaW')
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6

;
;     EMISSION
;

      if ny gt 1 then begin

      map = nadcube.weq[*,*,2]
      igd = where(abs(map) ge initmaps.nademweq_snrthresh*nadcube.weq[*,*,3] AND $
                  abs(map) gt 0d AND map ne bad)
      ibd = where(abs(map) lt initmaps.nademweq_snrthresh*nadcube.weq[*,*,3] OR $
                  map eq 0d OR map eq bad)

;     Set up range
;     Check for manual range first ...
      hasrange = 0
      if hasrangefile then begin
         ithisline = where(rangeline eq 'NaDem' AND $
            rangequant eq 'weq',ctthisline)
         if ctthisline eq 1 then begin
            zran = [rangelo[ithisline],rangehi[ithisline]]
            dzran = zran[1]-zran[0]
            ncbdiv = rangencbdiv[ithisline]
            ncbdiv = ncbdiv[0]
            hasrange = 1
         endif
      endif
;     otherwise set it automagically.
      if ~hasrange then begin
         zran = [min(map[igd]),max(map[igd])]
         divarr = ifsf_cbdiv(zran,2d,ncbdivmax)
         ncbdiv = divarr[0]
         dzran = zran[1]-zran[0]
      endif

;     replace bad points
      map[ibd] = bad

;     Save some things for later
      zran_empweqem = zran
      map_empweqem = map

      mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
         min=zran[0],max=zran[1])
      cgloadct,65
      cgimage,mapscl,/keep,pos=pos[*,2],opos=truepos,$
         noerase=i ne 0,missing_value=bad,missing_index=255,$
         missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
         /nodata,/noerase,title=textoidl('W_{eq}^{em}')+' ($\angstrom$)'
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6

;     Error
      err = nadcube.weq[*,*,3]
      map[igd] /= -err[igd]

;     Set up range
;     Check for manual range first ...
      hasrange = 0
      if hasrangefile then begin
         ithisline = where(rangeline eq 'NaDem' AND $
            rangequant eq 'weqsnr',ctthisline)
         if ctthisline eq 1 then begin
            zran = [rangelo[ithisline],rangehi[ithisline]]
            dzran = zran[1]-zran[0]
            ncbdiv = rangencbdiv[ithisline]
            ncbdiv = ncbdiv[0]
            hasrange = 1
         endif
      endif
;     otherwise set it automagically.
      if ~hasrange then begin
         zran = [min(map[igd]),max(map[igd])]
         divarr = ifsf_cbdiv(zran,2d,ncbdivmax)
         ncbdiv = divarr[0]
         dzran = zran[1]-zran[0]
      endif

      mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
         min=zran[0],max=zran[1])
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos[*,3],opos=truepos,$
         noerase=i ne 0,missing_value=bad,missing_index=255,$
         missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
         /nodata,/noerase,title=textoidl('W_{eq}^{em}/\deltaW')
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6



;     Flux
      cbform = '(D0.2)'
      map = nadcube.emflux[*,*,0]
      igd = where(abs(nadcube.weq[*,*,2]) ge $
                  initmaps.nademweq_snrthresh*nadcube.weq[*,*,3] AND $
                  abs(map) gt 0d AND map ne bad)
      ibd = where(abs(nadcube.weq[*,*,2]) lt $
                  initmaps.nademweq_snrthresh*nadcube.weq[*,*,3] OR $
                  map eq 0d OR map eq bad)

      if tag_exist(initmaps,'fluxfactor') then $
         map[igd] *= initmaps.fluxfactor
         
;     Set up range
;     Check for manual range first ...
      hasrange = 0
      if hasrangefile then begin
         ithisline = where(rangeline eq 'NaDem' AND $
            rangequant eq 'flux',ctthisline)
         if ctthisline eq 1 then begin
            zran = [rangelo[ithisline],rangehi[ithisline]]
            dzran = zran[1]-zran[0]
            ncbdiv = rangencbdiv[ithisline]
            ncbdiv = ncbdiv[0]
            hasrange = 1
         endif
      endif
;     otherwise set it automagically.
      if ~hasrange then begin
         zran = [min(map[igd]),max(map[igd])]
         divarr = ifsf_cbdiv(zran,2d,ncbdivmax)
         ncbdiv = divarr[0]
         dzran = zran[1]-zran[0]
      endif
      
;     replace bad points
      map[ibd] = bad

;     Save some things for later
      zran_empfluxem = zran
      map_empfluxem = map

      mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
         min=zran[0],max=zran[1])
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos[*,4],opos=truepos,$
         noerase=i ne 0,missing_value=bad,missing_index=255,$
         missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
         /nodata,/noerase,$
         title=textoidl('I^{em} / 10^{-16}')
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6

      endif

      cgps_close

   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD EQUIVALENT WIDTH (FITS)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if tag_exist(initdat,'donad') then begin

      if maxnademncomp_act gt 0 then ny=3 else ny=1

      cgps_open,initdat.mapdir+initdat.label+'NaDfitweq.eps',charsize=1,/encap,$
         /inches,xs=plotquantum*2,ys=plotquantum*ny,/qui

      pos = cglayout([2,ny],ixmar=[3,3],iymar=[3,3],oxmar=[0,0],oymar=[0,0],$
         xgap=0,ygap=0,unit=!D.X_PX_CM/3.0)
      cbform = '(I0)'

;
;     ABSORPTION
;
      map = nadfit.weqabs[*,*,0]
      maperrlo = nadfit.weqabserr[*,*,0]
      maperrhi = nadfit.weqabserr[*,*,1]
      maperravg = (maperrlo + maperrhi) / 2d
      maperrlo_emp = nadcube.weq[*,*,1]
      maperrhi_emp = nadcube.weq[*,*,1]
      igd = where(map gt 0d AND $
                  map ne bad AND $
                  map ge initmaps.nadabsweq_snrthresh*maperravg)
      ibd = where(map eq 0d OR $
                  map eq bad OR $
                  map lt initmaps.nadabsweq_snrthresh*maperravg)
      igd_nadfitabsweq = igd
      ibd_nadfitabsweq = ibd

;     Set up range
;     Check for manual range first ...
      hasrange = 0
      if hasrangefile then begin
         ithisline = where(rangeline eq 'NaDabs' AND $
                           rangequant eq 'weqfit',ctthisline)
         if ctthisline eq 1 then begin
            zran = [rangelo[ithisline],rangehi[ithisline]]
            dzran = zran[1]-zran[0]
            ncbdiv = rangencbdiv[ithisline]
            ncbdiv = ncbdiv[0]
            hasrange = 1
         endif
      endif
;     otherwise set it automagically.
      if ~hasrange then begin
         zran = [min(map[igd]),max(map[igd])]
         divarr = ifsf_cbdiv(zran,2d,ncbdivmax)
         ncbdiv = divarr[0]
         dzran = zran[1]-zran[0]
      endif
      zran_fitweqabs = zran

;     replace bad points
      map[ibd] = bad

      mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
                      min=zran[0],max=zran[1])
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos[*,0],opos=truepos,$
         missing_value=bad,missing_index=255,missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
         /nodata,/noerase,title=textoidl('W_{eq}^{abs}')+' ($\angstrom$)'
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
      ifsf_plotcompass,xarr_kpc,yarr_kpc
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6

      ibderr = where(map eq 0d OR map eq bad OR $
                     map_empweqabs eq 0d OR map_empweqabs eq bad)
      maperrlo[ibderr] = 0d
      maperrhi[ibderr] = 0d
      maperrlo_emp[ibderr] = 0d
      maperrhi_emp[ibderr] = 0d

      igd_empweqabs = where(map_empweqabs ne bad)
      xran = [0d,max([map[igd],map_empweqabs[igd_empweqabs]])*1.1d]
      yran = xran
      cgplot,map_empweqabs,map,/xsty,/ysty,$
         xran=xran,yran=yran,psym=3,pos=pos[*,1],$
         xtit=textoidl('W_{eq}^{abs} (empirical)'),$
         ytit=textoidl('W_{eq}^{abs} (fit)'),$
         /noerase,aspect=1d,$
         chars=0.8,err_width=0,err_color='Gray',/err_clip,$
         err_xlow=maperrlo_emp,err_xhigh=maperrhi_emp,$
         err_ylow=maperrlo,err_yhigh=maperrhi
      cgoplot,map_empweqabs,map,psym=16,symsize=0.4
      cgoplot,xran,xran
      inz = where(map ne 0d AND map ne bad AND $
                  map_empweqabs ne 0d AND map_empweqabs ne bad)
      weq_rms = sqrt(mean((map[inz]-map_empweqabs[inz])^2d))
      cgoplot,xran-weq_rms,xran,linesty=3
      cgoplot,xran+weq_rms,xran,linesty=3
;
;     EMISSION
;

      if ny gt 1 then begin

      map = nadfit.weqem[*,*,0]
      maperrlo = nadfit.weqemerr[*,*,0]
      maperrhi = nadfit.weqemerr[*,*,1]
      maperravg = (maperrlo + maperrhi)/2d
      maperrlo_emp = nadcube.weq[*,*,3]
      maperrhi_emp = nadcube.weq[*,*,3]
;      igd = where(abs(map) gt 0d AND map ne bad)
;      ibd = where(map eq 0d OR map eq bad)
      igd = where(abs(map) gt 0d AND $
                  map ne bad)
      ibd = where(map eq 0d OR $
                  map eq bad)
      igd_nadfitemweq = igd
      ibd_nadfitemweq = ibd

;     Set up range
;     Check for manual range first ...
      hasrange = 0
      if hasrangefile then begin
         ithisline = where(rangeline eq 'NaDem' AND $
            rangequant eq 'weq',ctthisline)
         if ctthisline eq 1 then begin
            zran = [rangelo[ithisline],rangehi[ithisline]]
            dzran = zran[1]-zran[0]
            ncbdiv = rangencbdiv[ithisline]
            ncbdiv = ncbdiv[0]
            hasrange = 1
         endif
      endif
;     otherwise set it automagically.
      if ~hasrange then begin
         zran = [min(map[igd]),max(map[igd])]
         divarr = ifsf_cbdiv(zran,2d,ncbdivmax)
         ncbdiv = divarr[0]
         dzran = zran[1]-zran[0]
      endif

;     replace bad points
      map[ibd] = 0d

      mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
         min=zran[0],max=zran[1])
      cgloadct,65
      cgimage,mapscl,/keep,pos=pos[*,2],opos=truepos,$
         /noerase,missing_value=bad,missing_index=255,$
         missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
         /nodata,/noerase,title=textoidl('W_{eq}^{em}')+' ($\angstrom$)'
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6

      ibderr = where(map eq 0d OR map eq bad OR $
                     map_empweqem eq 0d OR map_empweqem eq bad)
      maperrlo[ibderr] = 0d
      maperrhi[ibderr] = 0d
      maperrlo_emp[ibderr] = 0d
      maperrhi_emp[ibderr] = 0d

      igd_empweqem = where(map_empweqem ne bad)
      xran = [min([map[igd],map_empweqem[igd_empweqem]])*1.1d,0d]
      yran = xran
      cgplot,map_empweqem,map,/xsty,/ysty,xran=xran,yran=yran,psym=3,$
         pos=pos[*,3],/noerase,aspect=1d,$
         xtit=textoidl('W_{eq}^{em} (empirical)'),$
         ytit=textoidl('W_{eq}^{em} (fit)'),$    
         chars=0.8d,err_color='Gray',err_width=0,/err_clip,$
         err_xlow=maperrlo_emp,err_xhigh=maperrhi_emp,$
         err_ylow=maperrlo,err_yhigh=maperrhi
      cgoplot,map_empweqem,map,psym=16,symsize=0.4d
      cgoplot,xran,xran
      inz = where(map ne 0d AND map ne bad AND $
                  map_empweqem ne 0d AND map_empweqem ne bad)
      weq_rms = sqrt(mean((map[inz]-map_empweqem[inz])^2d))
      cgoplot,xran-weq_rms,xran,linesty=3
      cgoplot,xran+weq_rms,xran,linesty=3


;     Flux
      cbform = '(D0.2)'
      map = nadfit.totfluxem[*,*,0]
      maperrlo = nadfit.totfluxemerr[*,*,0]
      maperrhi = nadfit.totfluxemerr[*,*,1]
      maperrlo_emp = nadcube.emflux[*,*,1]
      maperrhi_emp = nadcube.emflux[*,*,1]
;      igd = where(abs(map) gt 0d AND map ne bad)
;      ibd = where(map eq 0d OR map eq bad)
      igd = igd_nadfitemweq
      ibd = ibd_nadfitemweq

      if tag_exist(initmaps,'fluxfactor') then $
         map[igd] *= initmaps.fluxfactor

;     Set up range
;     Check for manual range first ...
      hasrange = 0
      if hasrangefile then begin
         ithisline = where(rangeline eq 'NaDem' AND $
            rangequant eq 'flux',ctthisline)
         if ctthisline eq 1 then begin
            zran = [rangelo[ithisline],rangehi[ithisline]]
            dzran = zran[1]-zran[0]
            ncbdiv = rangencbdiv[ithisline]
            ncbdiv = ncbdiv[0]
            hasrange = 1
         endif
      endif
;     otherwise set it automagically.
      if ~hasrange then begin
         zran = [0,max(map[igd])]
         divarr = ifsf_cbdiv(zran,initmaps.nademflux_cbint,ncbdivmax)
         ncbdiv = divarr[0]
         dzran = zran[1]-zran[0]
      endif
      
;     replace bad points
      map[ibd] = bad

      mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
         min=zran[0],max=zran[1])
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos[*,4],opos=truepos,$
         /noerase,missing_value=bad,missing_index=255,$
         missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
         /nodata,/noerase,title='I$\upem$ / 10$\up-16$'
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6

      igderr = where(map ne 0d AND map ne bad AND $
                     map_empfluxem ne 0d AND map_empfluxem ne bad)
      ibderr = where(map eq 0d OR map eq bad OR $
                     map_empfluxem eq 0d OR map_empfluxem eq bad)
      maperrlo[ibderr] = 0d
      maperrhi[ibderr] = 0d
      maperrlo_emp[ibderr] = 0d
      maperrhi_emp[ibderr] = 0d
      if tag_exist(initmaps,'fluxfactor') then begin
         maperrlo[igderr] *= initmaps.fluxfactor
         maperrhi[igderr] *= initmaps.fluxfactor
         maperrlo_emp[igderr] *= initmaps.fluxfactor
         maperrhi_emp[igderr] *= initmaps.fluxfactor
      endif

      igd_empfluxem = where(map_empfluxem ne bad)
      xyran = [0d,max([map[igd],map_empfluxem[igd_empfluxem]])*1.1d]
      cgplot,map_empfluxem,map,/xsty,/ysty,xran=xyran,yran=xyran,psym=3,$
         pos=pos[*,5],xtit='I$\upem$ / 10$\up-16$ (empirical)',$
         ytit='I$\upem$ / 10$\up-16$ (fit)',/noerase,aspect=1d,$
         chars=0.8d,err_color='Gray',err_width=0,/err_clip,$
         err_xlow=maperrlo_emp,err_xhigh=maperrhi_emp,$
         err_ylow=maperrlo,err_yhigh=maperrhi
      cgoplot,map_empfluxem,map,psym=16,symsize=0.4d
      cgoplot,xyran,xyran
      inz = where(map ne 0d AND map ne bad AND $
                  map_empfluxem ne 0d AND map_empfluxem ne bad)
      weq_rms = sqrt(mean((map[inz]-map_empfluxem[inz])^2d))
      cgoplot,xyran-weq_rms,xyran,linesty=3
      cgoplot,xyran+weq_rms,xyran,linesty=3

      endif

      cgps_close

   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD components map (fits)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if tag_exist(initdat,'donad') then begin

;      cgps_open,initdat.mapdir+initdat.label+'NaDfitncomp.eps',charsize=1.5,$
;         /encap,$
;         /inches,xs=plotquantum*2,ys=plotquantum*2,/qui
;      map = nadabsncomp
;      ibd = where(map eq bad)
;      map[ibd] = 0l
;      
;      zran = [0,2]
;      ncbdiv = 2
;      dzran = zran[1]-zran[0]
;      cbform = '(I0)'
;
;      mapscl = cgimgscl(rebin(map,dx*20,dy*20,/sample),$
;                        minval=zran[0],maxval=zran[1],ncolors=3)
;      cgloadct,8,/brewer,ncolors=3
;      cgimage,mapscl,opos=truepos,mar=0.1d
;      cgplot,[0],xsty=5,ysty=5,position=truepos,$
;         /nodata,/noerase,title=textoidl('N_{comp}^{abs}')
;      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;;      cbpos=[pos[2,0],pos[1,0],pos[2,0]+0.02,pos[3,0]]
;      cbpos=[truepos[2],truepos[1],truepos[2]+0.02,truepos[3]]
;      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
;         (dzran - zran[1]),format=cbform)
;      cgDCBar,ncolors=3,position=cbpos,labels=ticknames,/ver,/right,$
;         charsize=1.5
;         
;      cgps_close

      cgps_open,initdat.mapdir+initdat.label+'NaDfitncomp.eps',charsize=1.5,$
         /encap,$
         /inches,xs=plotquantum*2,ys=plotquantum*2,/qui
      pos = cglayout([1,1],ixmar=[0,0],iymar=[0,0],oxmar=[3,5],oymar=[6,1],$
                     xgap=0,ygap=0,aspect=1,unit=!D.X_PX_CM/3.0)
      map = nadabsncomp
;      ibd = where(map eq bad)
      ibd = ibd_nadfitabsweq
      map[ibd] = 0l
      
      zran = [0,2]
      ncbdiv = 2
      dzran = zran[1]-zran[0]
      cbform = '(I0)'

      mapscl = cgimgscl(rebin(map,dx*20,dy*20,/sample),$
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
      
   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD "empirical" velocity maps
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if tag_exist(initdat,'donad') then begin

      if maxnademncomp_act gt 0 then ny=3 else ny=1

      cgps_open,initdat.mapdir+initdat.label+'NaDempvel.eps',charsize=1,/encap,$
         /inches,xs=plotquantum*3,ys=plotquantum*ny,/qui

      pos = cglayout([3,ny],ixmar=[3,3],iymar=[3,3],oxmar=[0,0],oymar=[0,0],$
         xgap=0,ygap=0,unit=!D.X_PX_CM/3.0)
      cbform = '(I0)'

;
;     ABSORPTION
;
      map = nadcube.vel[*,*,0]
      igd = where(nadcube.weq[*,*,0] ge initmaps.nadabsweq_snrthresh*nadcube.weq[*,*,1] $
                  AND map gt 0d AND map ne bad)
      ibd = where(nadcube.weq[*,*,0] lt initmaps.nadabsweq_snrthresh*nadcube.weq[*,*,1] $
                  OR map eq 0d OR map eq bad)

;     Set up range
;     Check for manual range first ...
      hasrange = 0
      if hasrangefile then begin
         ithisline = where(rangeline eq 'NaDabs' AND $
            rangequant eq 'empvelwid',ctthisline)
         if ctthisline eq 1 then begin
            zran = [rangelo[ithisline],rangehi[ithisline]]
            dzran = zran[1]-zran[0]
            ncbdiv = rangencbdiv[ithisline]
            ncbdiv = ncbdiv[0]
            hasrange = 1
         endif
      endif
;     otherwise set it automagically.
      if ~hasrange then begin
         zran = [min(map[igd]),max(map[igd])]
         divarr = ifsf_cbdiv(zran,200d,ncbdivmax)
         ncbdiv = divarr[0]
         dzran = zran[1]-zran[0]
      endif

;     replace bad points with "bad"
      map[ibd] = bad

      mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
         min=zran[0],max=zran[1])
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos[*,0],opos=truepos,$
         noerase=i ne 0,missing_value=bad,missing_index=255,$
         missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
         /nodata,/noerase,title='$\Delta$v(NaD abs, km/s)'
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6

      map = nadcube.vel[*,*,1]
      igd = where(nadcube.weq[*,*,0] ge initmaps.nadabsweq_snrthresh*nadcube.weq[*,*,1] $
                  AND map ne bad)
      ibd = where(nadcube.weq[*,*,0] lt initmaps.nadabsweq_snrthresh*nadcube.weq[*,*,1] $
                  OR map eq bad)

;     Set up range
;     Check for manual range first ...
      hasrange = 0
      if hasrangefile then begin
         ithisline = where(rangeline eq 'NaDabs' AND $
            rangequant eq 'empvelavg',ctthisline)
         if ctthisline eq 1 then begin
            zran = [rangelo[ithisline],rangehi[ithisline]]
            dzran = zran[1]-zran[0]
            ncbdiv = rangencbdiv[ithisline]
            ncbdiv = ncbdiv[0]
            hasrange = 1
         endif
      endif
;     otherwise set it automagically.
      if ~hasrange then begin
         zran = [min(map[igd]),max(map[igd])]
         divarr = ifsf_cbdiv(zran,100d,ncbdivmax)
         ncbdiv = divarr[0]
         dzran = zran[1]-zran[0]
      endif

;     replace bad points with "bad"
      map[ibd] = bad

      mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
         min=zran[0],max=zran[1])
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos[*,1],opos=truepos,$
         noerase=i ne 0,missing_value=bad,missing_index=255,$
         missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
         /nodata,/noerase,title='<v>(NaD abs, km/s)'
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6


      map = nadcube.vel[*,*,2]
      igd = where(nadcube.weq[*,*,0] ge initmaps.nadabsweq_snrthresh*nadcube.weq[*,*,1] $
                  AND map ne bad)
      ibd = where(nadcube.weq[*,*,0] lt initmaps.nadabsweq_snrthresh*nadcube.weq[*,*,1] $
                  OR map eq bad)

;     Set up range
;     Check for manual range first ...
      hasrange = 0
      if hasrangefile then begin
         ithisline = where(rangeline eq 'NaDabs' AND $
            rangequant eq 'empvelmax',ctthisline)
         if ctthisline eq 1 then begin
            zran = [rangelo[ithisline],rangehi[ithisline]]
            dzran = zran[1]-zran[0]
            ncbdiv = rangencbdiv[ithisline]
            ncbdiv = ncbdiv[0]
            hasrange = 1
         endif
      endif
;     otherwise set it automagically.
      if ~hasrange then begin
         zran = [min(map[igd]),max(map[igd])]
         divarr = ifsf_cbdiv(zran,100d,ncbdivmax)
         ncbdiv = divarr[0]
         dzran = zran[1]-zran[0]
      endif

;     replace bad points with "bad"
      map[ibd] = bad

      mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
         min=zran[0],max=zran[1])
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos[*,2],opos=truepos,$
         noerase=i ne 0,missing_value=bad,missing_index=255,$
         missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
         /nodata,/noerase,title='v$\downmax$(NaD abs, km/s)'
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6

;
;     EMISSION
;

      if ny gt 1 then begin

      map = nadcube.vel[*,*,3]
      igd = where(abs(nadcube.weq[*,*,2]) ge $
                  initmaps.nademweq_snrthresh*nadcube.weq[*,*,3] $
                  AND map gt 0d AND map ne bad)
      ibd = where(abs(nadcube.weq[*,*,2]) lt $
                  initmaps.nademweq_snrthresh*nadcube.weq[*,*,3] $
                  OR map eq 0d OR map eq bad)

;     Set up range
;     Check for manual range first ...
      hasrange = 0
      if hasrangefile then begin
         ithisline = where(rangeline eq 'NaDem' AND $
            rangequant eq 'empvelwid',ctthisline)
         if ctthisline eq 1 then begin
            zran = [rangelo[ithisline],rangehi[ithisline]]
            dzran = zran[1]-zran[0]
            ncbdiv = rangencbdiv[ithisline]
            ncbdiv = ncbdiv[0]
            hasrange = 1
         endif
      endif
;     otherwise set it automagically.
      if ~hasrange then begin
         zran = [min(map[igd]),max(map[igd])]
         divarr = ifsf_cbdiv(zran,200d,ncbdivmax)
         ncbdiv = divarr[0]
         dzran = zran[1]-zran[0]
      endif

;     replace bad points with "bad"
      map[ibd] = bad

      mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
         min=zran[0],max=zran[1])
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos[*,3],opos=truepos,$
         noerase=i ne 0,missing_value=bad,missing_index=255,$
         missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
         /nodata,/noerase,title='$\Delta$v(NaD em, km/s)'
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6

      map = nadcube.vel[*,*,4]
      igd = where(abs(nadcube.weq[*,*,2]) ge $
                  initmaps.nademweq_snrthresh*nadcube.weq[*,*,3] $
                  AND map ne bad)
      ibd = where(abs(nadcube.weq[*,*,2]) lt $
                  initmaps.nademweq_snrthresh*nadcube.weq[*,*,3] $
                  OR map eq bad)

;     Set up range
;     Check for manual range first ...
      hasrange = 0
      if hasrangefile then begin
         ithisline = where(rangeline eq 'NaDem' AND $
            rangequant eq 'empvelavg',ctthisline)
         if ctthisline eq 1 then begin
            zran = [rangelo[ithisline],rangehi[ithisline]]
            dzran = zran[1]-zran[0]
            ncbdiv = rangencbdiv[ithisline]
            ncbdiv = ncbdiv[0]
            hasrange = 1
         endif
      endif
;     otherwise set it automagically.
      if ~hasrange then begin
         zran = [min(map[igd]),max(map[igd])]
         divarr = ifsf_cbdiv(zran,100d,ncbdivmax)
         ncbdiv = divarr[0]
         dzran = zran[1]-zran[0]
      endif

;     replace bad points with "bad"
      map[ibd] = bad

      mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
         min=zran[0],max=zran[1])
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos[*,4],opos=truepos,$
         noerase=i ne 0,missing_value=bad,missing_index=255,$
         missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
         /nodata,/noerase,title='<v>(NaD em, km/s)'
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6

      map = nadcube.vel[*,*,5]
      igd = where(abs(nadcube.weq[*,*,2]) ge $
                  initmaps.nademweq_snrthresh*nadcube.weq[*,*,3] $
                  AND map ne bad)
      ibd = where(abs(nadcube.weq[*,*,2]) lt $
                  initmaps.nademweq_snrthresh*nadcube.weq[*,*,3] $
                  OR map eq bad)

;     Set up range
;     Check for manual range first ...
      hasrange = 0
      if hasrangefile then begin
         ithisline = where(rangeline eq 'NaDem' AND $
            rangequant eq 'empvelmax',ctthisline)
         if ctthisline eq 1 then begin
            zran = [rangelo[ithisline],rangehi[ithisline]]
            dzran = zran[1]-zran[0]
            ncbdiv = rangencbdiv[ithisline]
            ncbdiv = ncbdiv[0]
            hasrange = 1
         endif
      endif
;     otherwise set it automagically.
      if ~hasrange then begin
         zran = [min(map[igd]),max(map[igd])]
         divarr = ifsf_cbdiv(zran,100d,ncbdivmax)
         ncbdiv = divarr[0]
         dzran = zran[1]-zran[0]
      endif

;     replace bad points with "bad"
      map[ibd] = bad

      mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
         min=zran[0],max=zran[1])
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos[*,5],opos=truepos,$
         noerase=i ne 0,missing_value=bad,missing_index=255,$
         missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
         /nodata,/noerase,title='v$\downmax$(NaD em, km/s)'
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6

      endif

      cgps_close

   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; NaD fitted velocity maps
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if tag_exist(initdat,'donad') then begin

      nx = 3
      nyabs=0
      if maxnadabsncomp_act gt 0 then begin
         if donadabsonecomp then nyabs=1
         if donadabsmulticomp then nyabs+=maxnadabsncomp_act
      endif
      nyem=0
      if maxnademncomp_act gt 0 then begin
         if donademonecomp then nyem=1
         if donademmulticomp then nyem+=maxnademncomp_act
      endif
      ny = nyabs+nyem

      cgps_open,initdat.mapdir+initdat.label+'NaDfitvel.eps',charsize=1,/encap,$
                /inches,xs=plotquantum*nx,ys=plotquantum*ny,/qui

      pos = cglayout([nx,ny],ixmar=[2,3],iymar=[2,2],oxmar=[2,0],oymar=[0,2],$
                     xgap=0,ygap=0,unit=!D.X_PX_CM/3.0)
      cbform = '(I0)'
      abslab = ['Abs (1 Comp)']
      if nyabs gt 1 then $
         for i=1,nyabs-1 do $
            abslab = [abslab,string('Abs (',i,' of ',nyabs-1,' Comp)',$
                                    format=('(A0,I0,A0,I0,A0)'))]
      emlab = ['Em (1 Comp)']
      if nyem gt 1 then $
         for i=1,nyem-1 do $
            emlab = [emlab,string('Em (',i,' of ',nyem-1,' Comp)',$
                                  format=('(A0,I0,A0,I0,A0)'))]

;
;     ABSORPTION
;
      for i=0,nyabs-1 do begin
;        sigma
         map = nadabssig[*,*,i]
;         igd = where(map ne bad)
;         ibd = where(map eq bad)
         igd_thiscomp = where(map ne bad)
         ibd_thiscomp = where(map eq bad)
         igd = cgsetintersection(igd_nadfitabsweq,igd_thiscomp)
         ibd = cgsetunion(ibd_nadfitabsweq,ibd_thiscomp)
         map[ibd] = bad
         plotdat = ifsf_plotrange(/auto,mapgd=map[igd],divinit=100d,$
                                  ncbdivmax=ncbdivmax)
         mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
                         min=plotdat[0],max=plotdat[1])
         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos[*,0+i*nx],opos=truepos,$
                 noerase=i ne 0,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
         if i eq 0 then begin
            cgtext,(xran_kpc[0]+xran_kpc[1])/2d,$
                   yran_kpc[1]+0.1*(yran_kpc[1]-yran_kpc[0]),$
                  '$\sigma$ (km/s)',chars=1.25,align=0.5
            ifsf_plotcompass,xarr_kpc,yarr_kpc
         endif
         cgtext,xran_kpc[0]-0.17*(xran_kpc[1]-xran_kpc[0]),$
                (yran_kpc[0]+yran_kpc[1])/2d,$
                abslab[i],align=0.5,orient=90d,chars=1.25
         cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,/ver,/right,charsize=0.6
         
;        v50
         map = nadabsvel[*,*,i]
;         igd = where(map ne bad)
;         ibd = where(map eq bad)
         igd_thiscomp = where(map ne bad)
         ibd_thiscomp = where(map eq bad)
         igd = cgsetintersection(igd_nadfitabsweq,igd_thiscomp)
         ibd = cgsetunion(ibd_nadfitabsweq,ibd_thiscomp)
         map[ibd] = bad
         plotdat = ifsf_plotrange(/auto,mapgd=map[igd],divinit=200d,$
                                  ncbdivmax=ncbdivmax)
         mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
                         min=plotdat[0],max=plotdat[1])
         cgloadct,22,/brewer,/reverse
;         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos[*,1+i*nx],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
         if i eq 0 then $
            cgtext,(xran_kpc[0]+xran_kpc[1])/2d,$
                   yran_kpc[1]+0.1*(yran_kpc[1]-yran_kpc[0]),$
                  'v$\down50$ (km/s)',chars=1.25,align=0.5
         cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,/ver,/right,charsize=0.6
;        v98
         map = nadabsv98[*,*,i]
;         igd = where(map ne bad)
;         ibd = where(map eq bad)
         igd_thiscomp = where(map ne bad)
         ibd_thiscomp = where(map eq bad)
         igd = cgsetintersection(igd_nadfitabsweq,igd_thiscomp)
         ibd = cgsetunion(ibd_nadfitabsweq,ibd_thiscomp)
         map[ibd] = bad
         plotdat = ifsf_plotrange(/auto,mapgd=map[igd],divinit=200d,$
                                  ncbdivmax=ncbdivmax)
         mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
            min=plotdat[0],max=plotdat[1])
         cgloadct,22,/brewer,/reverse
;         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos[*,2+i*nx],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
         if i eq 0 then $
            cgtext,(xran_kpc[0]+xran_kpc[1])/2d,$
                   yran_kpc[1]+0.1*(yran_kpc[1]-yran_kpc[0]),$
                  'v$\down98$ (km/s)',chars=1.25,align=0.5
         cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,/ver,/right,charsize=0.6
;
      endfor
;
;     EMISSION
;
      for i=0,nyem-1 do begin
;        sigma
         map = nademsig[*,*,i]
;         igd = where(map ne bad)
;         ibd = where(map eq bad)
         igd_thiscomp = where(map ne bad)
         ibd_thiscomp = where(map eq bad)
         igd = cgsetintersection(igd_nadfitemweq,igd_thiscomp)
         ibd = cgsetunion(ibd_nadfitemweq,ibd_thiscomp)
         map[ibd] = bad
         plotdat = ifsf_plotrange(/auto,mapgd=map[igd],divinit=200d,$
                                  ncbdivmax=ncbdivmax)
         mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
            min=plotdat[0],max=plotdat[1])
         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos[*,nyabs*nx+i*nx],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
         cgtext,xran_kpc[0]-0.17*(xran_kpc[1]-xran_kpc[0]),$
                (yran_kpc[0]+yran_kpc[1])/2d,$
                emlab[i],align=0.5,orient=90d,chars=1.25
         cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,/ver,/right,charsize=0.6
;        v50
         map = nademvel[*,*,i]
;         igd = where(map ne bad)
;         ibd = where(map eq bad)
         igd_thiscomp = where(map ne bad)
         ibd_thiscomp = where(map eq bad)
         igd = cgsetintersection(igd_nadfitemweq,igd_thiscomp)
         ibd = cgsetunion(ibd_nadfitemweq,ibd_thiscomp)
         map[ibd] = bad
         plotdat = ifsf_plotrange(/auto,mapgd=map[igd],divinit=200d,$
                                  ncbdivmax=ncbdivmax)
         mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
            min=plotdat[0],max=plotdat[1])
         cgloadct,22,/brewer,/reverse
;         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos[*,1+nyabs*nx+i*nx],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
         cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,/ver,/right,charsize=0.6
;        v98
         map = nademv98[*,*,i]
;         igd = where(map ne bad)
;         ibd = where(map eq bad)
         igd_thiscomp = where(map ne bad)
         ibd_thiscomp = where(map eq bad)
         igd = cgsetintersection(igd_nadfitemweq,igd_thiscomp)
         ibd = cgsetunion(ibd_nadfitemweq,ibd_thiscomp)
         map[ibd] = bad
         plotdat = ifsf_plotrange(/auto,mapgd=map[igd],divinit=200d,$
                                  ncbdivmax=ncbdivmax)
         mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
            min=plotdat[0],max=plotdat[1])
         cgloadct,22,/brewer,/reverse
;         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos[*,2+nyabs*nx+i*nx],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
         cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,/ver,/right,charsize=0.6
;
      endfor

      cgps_close

   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; NaD abs. fit tau and C_f
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if tag_exist(initdat,'donad') then begin

      nx = 2
      if donadabsonecomp then ny=1
      if donadabsmulticomp then ny+=maxnadabsncomp_act

      cgps_open,initdat.mapdir+initdat.label+'NaDfit_taucf.eps',charsize=1,/encap,$
                /inches,xs=plotquantum*nx,ys=plotquantum*ny,/qui

      pos = cglayout([nx,ny],ixmar=[2,3],iymar=[2,2],oxmar=[2,0],oymar=[0,2],$
                     xgap=0,ygap=0,unit=!D.X_PX_CM/3.0)
      abslab = ['Abs (1 Comp)']
      if nyabs gt 1 then $
         for i=1,nyabs-1 do $
            abslab = [abslab,string('Abs (',i,' of ',nyabs-1,' Comp)',$
                                    format=('(A0,I0,A0,I0,A0)'))]

     for i=0,ny-1 do begin
;        tau
         map = nadabstau[*,*,i]
;         igd = where(map ne bad)
;         ibd = where(map eq bad)
         igd_thiscomp = where(map ne bad)
         ibd_thiscomp = where(map eq bad)
         igd = cgsetintersection(igd_nadfitabsweq,igd_thiscomp)
         ibd = cgsetunion(ibd_nadfitabsweq,ibd_thiscomp)
         map[ibd]=bad
         plotdat = ifsf_plotrange(/auto,mapgd=map[igd],divinit=0.5d,$
                                  ncbdivmax=ncbdivmax)
         mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
                         min=plotdat[0],max=plotdat[1])
         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos[*,0+i*nx],opos=truepos,$
                 noerase=i ne 0,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
         if i eq 0 then begin
            cgtext,(xran_kpc[0]+xran_kpc[1])/2d,$
                   yran_kpc[1]+0.1*(yran_kpc[1]-yran_kpc[0]),$
                  '$\tau$',chars=1.25,align=0.5
            ifsf_plotcompass,xarr_kpc,yarr_kpc
         endif
         cgtext,xran_kpc[0]-0.17*(xran_kpc[1]-xran_kpc[0]),$
                (yran_kpc[0]+yran_kpc[1])/2d,$
                abslab[i],align=0.5,orient=90d,chars=1.25
         cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         cbform = '(I0)'
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,/ver,/right,charsize=0.6
         
;        C_f
         map = nadabscf[*,*,i]
;         igd = where(map ne bad)
;         ibd = where(map eq bad)
         igd_thiscomp = where(map ne bad)
         ibd_thiscomp = where(map eq bad)
         igd = cgsetintersection(igd_nadfitabsweq,igd_thiscomp)
         ibd = cgsetunion(ibd_nadfitabsweq,ibd_thiscomp)
         map[ibd]=bad
         plotdat = ifsf_plotrange(/auto,mapgd=map[igd],divinit=0.2d,$
                                  ncbdivmax=ncbdivmax)
         mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
                         min=plotdat[0],max=plotdat[1])
         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos[*,1+i*nx],opos=truepos,$
                 /noerase,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
         if i eq 0 then $
            cgtext,(xran_kpc[0]+xran_kpc[1])/2d,$
                   yran_kpc[1]+0.1*(yran_kpc[1]-yran_kpc[0]),$
                  'C$\downf$',chars=1.25,align=0.5
         cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
         ticknames = string(dindgen(plotdat[3]+1)*plotdat[2]/double(plotdat[3]) - $
                            (plotdat[2] - plotdat[1]),format=cbform)
         cbform = '(D0.1)'
         cgcolorbar,position=cbpos,divisions=plotdat[3],$
                    ticknames=ticknames,/ver,/right,charsize=0.6
      endfor

      cgps_close

   endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Outflow maps
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if tag_exist(initmaps,'compof') AND $
    ~ tag_exist(initmaps,'noemlinfit') then begin

;     Quantities to plot
      if tag_exist(initmaps,'cvdf_oftags') then tags=initmaps.cvdf_oftags $
      else tags=['f','sig','v50','v84','v98']
      if tag_exist(initmaps,'cvdf_oftagtypes') then $
         tagtypes=initmaps.cvdf_oftagtypes $
      else tagtypes=['flux','vel','vel','vel','vel']
      if tag_exist(initmaps,'cvdf_oftitles') then $
         titles=initmaps.cvdf_oftitles $
      else titles = ['F$\downtot$','$\sigma$','v$\down50$',$
                     'v$\down84$','v$\down98$']

;     Size of plot grid
      if tag_exist(initmaps,'cvdf_ofnpx') then npx=initmaps.cvdf_ofnpx $
      else npx=3
      if tag_exist(initmaps,'cvdf_ofnpy') then npy=initmaps.cvdf_ofnpy $
      else npy=2
      
;     Loop through emission lines
      foreach line,outlines do begin

;        Get syntax of linelabel right; otherwise call to DEVICE chokes
         linelab=line
         ilb = strpos(linelab,'[')
         if ilb ne -1 then $
            linelab = strmid(linelab,0,ilb)+'\'+strmid(linelab,ilb)
         irb = strpos(linelab,']')
         if irb ne -1 then $
            linelab = strmid(linelab,0,irb)+'\'+strmid(linelab,irb)

         cgps_open,initdat.mapdir+initdat.label+linelab+'_of_c.eps',charsize=1,/encap,$
                   /inches,xs=plotquantum*npx,ys=plotquantum*npy*aspectrat,/qui
         pos = cglayout([npx,npy],ixmar=[3,3],iymar=[3,3],oxmar=[0,0],oymar=[0,0],$
                        xgap=0,ygap=0,unit=!D.X_PX_CM/3.0)

;        loop through plot types
         for j=0,n_elements(tags)-1 do begin

            iplot = j ; plot index
;           Get map and scale
            itag = where(strcmp(ofpars_tags,tags[j],/fold_case) eq 1,cttag)
            if cttag ne 0 then begin

            map = ofpars[line].(itag)
            ibd = where(map eq bad AND ~ finite(map),ctbd)
            inan = where(~finite(map),ctnan)
            igd = where(map ne bad AND map ne 0 AND finite(map),ctgd)

            title=titles[j]
            
            if ctgd gt 0 then begin
               
               if tagtypes[j] eq 'flux' then begin
 
                  cbform = '(D0.1)' ; colorbar syntax
               
                  if tag_exist(initmaps,'fluxfactor') then $
                     map[igd] *= initmaps.fluxfactor
                                 
                  zran=[0,1]
                  dzran = 1
                  ncbdiv = 5
                  zmax_flux = max(map[igd])

                  if hasrangefile then begin
                     ithisline = where(line eq rangeline AND $
                                       'of_'+tags[j] eq rangequant,ctthisline)
                     if ctthisline eq 1 then zmax_flux = rangelo[ithisline]
                  endif
                
                  map[igd] = map[igd]/zmax_flux[0]
                  title += ' ('+string(zmax_flux,format='(E0.2)')+')'
                  cgloadct,65,/reverse

               endif else begin
       
                  cbform = '(I0)' ; colorbar syntax
                  hasrange = 0
                  if hasrangefile then begin
                     ithisline = where(line eq rangeline AND $
                                       'of_'+tags[j] eq rangequant,ctthisline)
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

                  cgloadct,74,/reverse
 
               endelse
 
               if ctnan gt 0 then map[inan] = bad
               mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
                               min=zran[0],max=zran[1])

;              Plot image
               cgimage,mapscl,/keep,pos=pos[*,iplot],opos=truepos,$
                       noerase=iplot ne 0,missing_value=bad,missing_index=255,$
                       missing_color='white'
               cgplot,[0],xsty=5,ysty=5,position=truepos,$
                         /nodata,/noerase,title=title
               ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
;              Colorbar
               cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
               ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
                                  (dzran - zran[1]),format=cbform)
               cgcolorbar,position=cbpos,divisions=ncbdiv,$
                          ticknames=ticknames,/ver,/right,charsize=0.6

            endif
            endif
            
         endfor
 
         cgps_close

      endforeach

   endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Compute outflow properties
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   m_ha = 0d
   dmdt_ha = 0d
   p_ha = 0d
   dpdt_ha = 0d
   e_ha = 0d
   dedt_ha = 0d
   
   m_nad = 0d
   dmdt_nad = 0d
   p_nad = 0d
   dpdt_nad = 0d
   e_nad = 0d
   dedt_nad = 0d

   if tag_exist(initmaps,'compof') then begin

      ;;;;;;;;;;;
      ; IONIZED ;
      ;;;;;;;;;;;

;     Use CVDF to get outflow properties

;     Blueshifted
      map_fha = ofpars['Halpha'].f
      igd = where(map_fha ne bad,ctgd)
;     Convert from flux/arcsec^2 back to straight flux if requested
      if tag_exist(initmaps,'fluxunfactor') AND ctgd gt 0 then $
         map_fha[igd] *= initmaps.fluxunfactor
      map_vha = ofpars['Halpha'].v50*1d5
      map_sigha = ofpars['Halpha'].sig*1d5
                                ; velocities in cm/s

;     Redshifted
      mapr_fha = ofpars['Halpha'].fr
      igd = where(mapr_fha ne bad,ctgd)
;     Convert from flux/arcsec^2 back to straight flux if requested
      if tag_exist(initmaps,'fluxunfactor') AND ctgd gt 0 then $
         mapr_fha[igd] *= initmaps.fluxunfactor
      mapr_vha = ofpars['Halpha'].v50r*1d5
      mapr_sigha = ofpars['Halpha'].sigr*1d5
                                ; velocities in cm/s

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

;     Radius of wind in cm
      Rcm_ha = 3.08567802d18 * initmaps.Rha * 1d3
;     Radius of wind in spaxels
      Rpix_ha = (initmaps.Rha / kpc_per_as)/initdat.platescale
  
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
                       map_fha gt 0,ctgood_ha)
      igoodr_ha = where(map_r+0.5 lt Rpix_ha AND $
                        mapr_fha ne bad AND $
                        mapr_fha gt 0,ctgood_ha)

; blueshifted
      map_vrad_ha[igood_ha] = abs(map_vha[igood_ha]) / $
                                  map_costheta_ha[igood_ha]
; 1d-15 here is for Gemini default flux normalization; other 1d-7 is to go
; from ergs/s to W
      map_l_ha[igood_ha] = linelum(map_fha[igood_ha]*1d-22,ldist,/ergs)
      map_m_ha[igood_ha] = mumpsm * map_l_ha[igood_ha] / $
                           volemis / elecden
                                ; in units of msun
; The following formula comes from dividing dM/dt^avg_thin by M_thin
; (eq. 7 and 8 in RVS05b). Note that this is basically equivalent to
; dividing each pixel by its own dynamical time, tdyn ~ v/R.
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

; redshifted
      mapr_vrad_ha[igoodr_ha] = abs(mapr_vha[igoodr_ha]) / $
                                   map_costheta_ha[igoodr_ha]
      mapr_l_ha[igoodr_ha] = linelum(mapr_fha[igoodr_ha]*1d-22,ldist,/ergs)
      mapr_m_ha[igoodr_ha] = mumpsm * mapr_l_ha[igoodr_ha] / $
                             volemis / elecden
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


      m_ha = 2d*(total(map_m_ha)+total(mapr_m_ha))
      dmdt_ha = 2d*(total(map_dmdt_ha)+total(mapr_dmdt_ha))
      p_ha = 2d*(total(map_p_ha)+total(mapr_p_ha))
      dpdt_ha = 2d*(total(map_dpdt_ha)+total(mapr_dpdt_ha))
      e_ha = 2d*(total(map_e_ha)+total(mapr_e_ha))
      dedt_ha = 2d*(total(map_dedt_ha)+total(mapr_dedt_ha))
      omega_ha = 2d/4/!DPi*(total(map_domega_ha[igood_ha])+$
                 total(map_domega_ha[igoodr_ha]))

      ;;;;;;;;;;;
      ; NEUTRAL ;
      ;;;;;;;;;;;

;     Sum over outflowing components to get outflow properties

      if tag_exist(initdat,'donad') then begin

         ;;;;;;;;;;;;;;
         ; CRFW model ;
         ;;;;;;;;;;;;;;
  
         Rcm_nad = 3.08567802d18 * initmaps.Rnad * 1d3
            ; radius of wind in cm
         Rpix_nad = (initmaps.Rnad / kpc_per_as)/initdat.platescale
            ; radius of wind in spaxels
         map_dphi_sintheta_nad = 1/Rpix_nad
         map_dtheta_nad = asin((map_r+0.5)/Rpix_nad) - asin((map_r-0.5)/Rpix_nad)
         map_domega_nad = map_dphi_sintheta_nad*map_dtheta_nad
         map_costheta_nad = cos(asin(map_r/Rpix_nad))

;        Cycle through components
         for i=0,maxnadabsncomp_act do begin
            cnh =  nadabscnh[*,*,i]
            cf = nadabscf[*,*,i]
            vel = nadabsvel[*,*,i]*1d5
            sig = nadabssig[*,*,i]*1d5
            igdof = where(cnh ne bad,ctgd)
            if ctgd gt 0 then begin

               map_vrad_comp_nad = dblarr(dx,dy)
               map_sig_comp_nad = dblarr(dx,dy)
               map_m_comp_nad = dblarr(dx,dy)
               map_dmdt_comp_nad = dblarr(dx,dy)
               map_p_comp_nad = dblarr(dx,dy)
               map_dpdt_comp_nad = dblarr(dx,dy)
               map_e_comp_nad = dblarr(dx,dy)
               map_dedt_comp_nad = dblarr(dx,dy)

               map_vrad_comp_nad[igdof] = $
                  abs(vel[igdof]) / map_costheta_nad[igdof]
               map_sig_comp_nad[igdof] = sig[igdof]
               map_m_comp_nad[igdof] = $
                  mumpsm*Rcm_nad^2d*cnh[igdof]*map_domega_nad[igdof]*cf[igdof]
               map_dmdt_comp_nad[igdof] = $
                  mumpsm*Rcm_nad*cnh[igdof]*map_vrad_comp_nad[igdof]*speryr*$
                  map_domega_nad[igdof]*cf[igdof]
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


               m_nad += 2d*total(map_m_comp_nad)
               dmdt_nad += 2d*total(map_dmdt_comp_nad)
               p_nad += 2d*total(map_p_comp_nad)
               dpdt_nad += 2d*total(map_dpdt_comp_nad)
               e_nad += 2d*total(map_e_comp_nad)
               dedt_nad += 2d*total(map_dedt_comp_nad)

            endif
 
         endfor
      
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
         string(qso_fitpar[2]*2d,format='(D0.4)'),$
         ' spaxels'
      printf,lun_stats,'  FWHM = ',$
         string(qso_fitpar[2]*2d*initdat.platescale,format='(D0.4)'),$
         ' arcseconds'
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
         string(empsf_fitpar[2]*2d,format='(D0.4)'),$
         ' spaxels'
      printf,lun_stats,'  FWHM = ',$
         string(empsf_fitpar[2]*2d*initdat.platescale,format='(D0.4)'),$
         ' arcseconds'
      printf,lun_stats,'  index = ',$
         string(empsf_fitpar[7],format='(D0.4)')
   endif
   did_distance_intro = 0b
   if dohst then begin
      if tag_exist(initmaps.hst,'fithstpeak') then begin
         printf,lun_stats,'-----','------------','--------','--------','--------',$
                          '--------','--------','-----',format='(A-5,A12,5A8,A5)'
         printf,lun_stats,'Distances between assumed peak and actual (fitted) peak'
         printf,lun_stats,'-----','------------','--------','--------','--------',$
                          '--------','--------','-----',format='(A-5,A12,5A8,A5)'
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
            printf,lun_stats,'-----','------------','--------','--------','--------',$
                             '--------','--------','-----',format='(A-5,A12,5A8,A5)'
            printf,lun_stats,'Distances between assumed peak and actual (fitted) peak'
            printf,lun_stats,'-----','------------','--------','--------','--------',$
                                '--------','--------','-----',format='(A-5,A12,5A8,A5)'
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
            printf,lun_stats,'-----','------------','--------','--------','--------',$
                             '--------','--------','-----',format='(A-5,A12,5A8,A5)'
            printf,lun_stats,'Distances between assumed peak and actual (fitted) peak'
            printf,lun_stats,'-----','------------','--------','--------','--------',$
                             '--------','--------','-----',format='(A-5,A12,5A8,A5)'
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
 
   if tag_exist(initdat,'donad') then begin

   printf,lun_stats,'-----','------------','--------','--------','--------',$
          '--------','--------','-----',format='(A-5,A12,5A8,A5)'
   printf,lun_stats,'ABSORPTION'
   printf,lun_stats,'-----','------------','--------','--------','--------',$
          '--------','--------','-----',format='(A-5,A12,5A8,A5)'
   printf,lun_stats,'#Cmp','Quantity','Mean','Median','Min','Max','StdDev','#',$
          format='(A-5,A12,5A8,A5)'
   printf,lun_stats,'-----','------------','--------','--------','--------',$
          '--------','--------','-----',format='(A-5,A12,5A8,A5)'
   icomp = 0
   arr = cvd_nad_pars.v50
   igd = where(arr ne bad,ctgd)
   arr = arr[igd]
   printf,lun_stats,string(icomp,format='(I0)'),'Vel',mean(arr),median(arr),$
          min(arr),max(arr),stddev(arr),ctgd,$
          format='(A-5,A12,5D8.0,I5)'
   a_vel_stats = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]
   arr = cvd_nad_pars.sig
   igd = where(arr ne bad)
   arr = arr[igd] * 2d * sqrt(2d*alog(2d))
   printf,lun_stats,string(icomp,format='(I0)'),'FWHM',mean(arr),median(arr),$
          min(arr),max(arr),stddev(arr),ctgd,$
          format='(A-5,A12,5D8.0,I5)'
   a_fwhm_stats = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]
   arr = cvd_nad_pars.v98
   igd = where(arr ne bad)
   arr = arr[igd]
   printf,lun_stats,string(icomp,format='(I0)'),'V98',mean(arr),median(arr),$
          min(arr),max(arr),stddev(arr),ctgd,$
          format='(A-5,A12,5D8.0,I5)'
   a_v98_stats = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]

   if tag_exist(initmaps,'compof') then begin
      printf,lun_stats,$
             '-----','------------','--------','--------','--------',$
             '--------','--------','-----',format='(A-5,A12,5A8,A5)'
      printf,lun_stats,'CRFW model'
      printf,lun_stats,$
             '-----','------------','--------','--------','--------','--------',$
             '--------','-----',format='(A-5,A12,5A8,A5)'
;        printf,lun_stats,$
;           'Total solid angle: ',string(omega_nad,format='(D0.3)'),$
;           ' [fraction of 4Pi]'
;        printf,lun_stats,$
;           'Total solid angle times Cf: ',$
;           string(omega_timescf_nad,format='(D0.3)'),' [fraction of 4Pi]'
      printf,lun_stats,$
             'Total mass: ',string(m_nad,format='(E0.2)'),' M_sun'
      printf,lun_stats,$
             'Total mass outflow rate: ',string(dmdt_nad,format='(E0.2)'),$
             ' M_sun/yr'
      printf,lun_stats,$
             'Total momentum: ',string(m_nad,format='(E0.2)'),' dyne s'
      printf,lun_stats,$
             'Total momentum outflow rate * c: ',$
             string(dpdt_nad,format='(E0.2)'),' Lsun'
      printf,lun_stats,$
             'Total energy: ',string(e_nad,format='(E0.2)'),' erg'
      printf,lun_stats,$
             'Total energy outflow rate: ',$
             string(dedt_nad,format='(E0.2)'),' erg/s'
   endif

   endif else begin
      
      a_vel_stats = dblarr(6)+bad
      a_fwhm_stats = dblarr(6)+bad
      a_v98_stats = dblarr(6)+bad
      
   endelse

   if ~ tag_exist(initmaps,'ofparline') then ofparline='Halpha' $
   else ofparline=initmaps.ofparline

   if ~ tag_exist(initdat,'noemlinfit') then begin

   printf,lun_stats,'-----','------------','--------','--------','--------',$
                    '--------','--------','-----',format='(A-5,A12,5A8,A5)'
   printf,lun_stats,'EMISSION'
   printf,lun_stats,'-----','------------','--------','--------','--------',$
                    '--------','--------','-----',format='(A-5,A12,5A8,A5)'
   printf,lun_stats,'CVDF: total'
   printf,lun_stats,'-----','------------','--------','--------','--------',$
                    '--------','--------','-----',format='(A-5,A12,5A8,A5)'
   icomp = 0
   arr = linspecpars[ofparline].v50
   igd = where(arr ne bad,ctgd)
   arr = arr[igd]
   printf,lun_stats,string(icomp,format='(I0)'),'Vel',mean(arr),median(arr),$
          min(arr),max(arr),stddev(arr),ctgd,format='(A-5,A12,5D8.0,I5)'
   e_vel_stats_all = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]
   arr = linspecpars[ofparline].sig
   igd = where(arr ne bad,ctgd)
   arr = arr[igd] * 2d * sqrt(2d*alog(2d))
   printf,lun_stats,string(icomp,format='(I0)'),'FWHM',mean(arr),median(arr),$
          min(arr),max(arr),stddev(arr),ctgd,format='(A-5,A12,5D8.0,I5)'
   e_fwhm_stats_all = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]
   arr = linspecpars[ofparline].v98
   igd = where(arr ne bad,ctgd)
   arr = arr[igd]
   printf,lun_stats,string(icomp,format='(I0)'),'V98',mean(arr),median(arr),$
          min(arr),max(arr),stddev(arr),ctgd,format='(A-5,A12,5D8.0,I5)'
   e_v98_stats_all = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]

   if tag_exist(initmaps,'compof') then begin
      printf,lun_stats,'-----','------------','--------','--------','--------',$
                       '--------','--------','-----',format='(A-5,A12,5A8,A5)'
      printf,lun_stats,'CVDF: outflow only'
      printf,lun_stats,'-----','------------','--------','--------','--------',$
                       '--------','--------','-----',format='(A-5,A12,5A8,A5)'
      icomp = 0
      arr = ofpars[ofparline].v50
      igd = where(arr ne bad,ctgd)
      arr = arr[igd]
      printf,lun_stats,string(icomp,format='(I0)'),'Vel',mean(arr),median(arr),$
             min(arr),max(arr),stddev(arr),ctgd,format='(A-5,A12,5D8.0,I5)'
      e_vel_stats = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]
      arr = ofpars[ofparline].sig
      igd = where(arr ne bad,ctgd)
      arr = arr[igd] * 2d * sqrt(2d*alog(2d))
      printf,lun_stats,string(icomp,format='(I0)'),'FWHM',mean(arr),median(arr),$
             min(arr),max(arr),stddev(arr),ctgd,format='(A-5,A12,5D8.0,I5)'
      e_fwhm_stats = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]
      arr = ofpars[ofparline].v98
      igd = where(arr ne bad,ctgd)
      arr = arr[igd]
      printf,lun_stats,string(icomp,format='(I0)'),'V98',mean(arr),median(arr),$
             min(arr),max(arr),stddev(arr),ctgd,format='(A-5,A12,5D8.0,I5)'
      e_v98_stats = [mean(arr),median(arr),stddev(arr),min(arr),max(arr),ctgd]
  
      printf,lun_stats,$
             '-----','------------','--------','--------','--------',$
             '--------','--------','-----',format='(A-5,A12,5A8,A5)'
      printf,lun_stats,'CRFW model'
      printf,lun_stats,$
             '-----','------------','--------','--------','--------','--------',$
             '--------','-----',format='(A-5,A12,5A8,A5)'
      printf,lun_stats,$
             'Total mass: ',string(m_ha,format='(E0.2)'),' M_sun'
      printf,lun_stats,$
             'Total mass outflow rate: ',string(dmdt_ha,format='(E0.2)'),$
             ' M_sun/yr'
      printf,lun_stats,$
             'Total momentum: ',string(p_ha,format='(E0.2)'),' dyne s'
      printf,lun_stats,$
             'Total momentum outflow rate * c: ',$
             string(dpdt_ha,format='(E0.2)'),' Lsun'
      printf,lun_stats,$
             'Total energy: ',string(e_ha,format='(E0.2)'),' erg'
      printf,lun_stats,$
             'Total energy outflow rate: ',$
             string(dedt_ha,format='(E0.2)'),' erg/s'
   endif else begin
            
      e_vel_stats = dblarr(6)+bad
      e_fwhm_stats = dblarr(6)+bad
      e_v98_stats = dblarr(6)+bad
      e_vel_stats_all = dblarr(6)+bad
      e_fwhm_stats_all = dblarr(6)+bad
      e_v98_stats_all = dblarr(6)+bad
      
   endelse

   endif else begin

   e_vel_stats = dblarr(6)+bad
   e_fwhm_stats = dblarr(6)+bad
   e_v98_stats = dblarr(6)+bad
   e_vel_stats_all = dblarr(6)+bad
   e_fwhm_stats_all = dblarr(6)+bad
   e_v98_stats_all = dblarr(6)+bad

   endelse

   free_lun,lun_stats

   if tag_exist(initmaps,'compof') then nem = 0 else nem = -1
   windstr = {$
              nem: nem,$
;             Statistics
              a_vel_stats: a_vel_stats,$
              a_fwhm_stats: a_fwhm_stats,$
              a_v98_stats: a_v98_stats, $
              e_vel_stats_all: e_vel_stats_all,$
              e_fwhm_stats_all: e_fwhm_stats_all,$
              e_v98_stats_all: e_v98_stats_all,$
              e_vel_stats: e_vel_stats,$
              e_fwhm_stats: e_fwhm_stats,$
              e_v98_stats: e_v98_stats,$
;             CRFW results
              a_m:m_nad,$
              a_dmdt:dmdt_nad,$
              a_p:p_nad,$
              a_dpdt:dpdt_nad,$
              a_e:e_nad,$
              a_dedt:dedt_nad,$
              e_m:m_ha,$
              e_dmdt:dmdt_ha,$
              e_p:p_ha,$
              e_dpdt:dpdt_ha,$
              e_e:e_ha,$
              e_dedt:dedt_ha $
             }
   save,windstr,file=initdat.mapdir+initdat.label+'.xdr'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; OTHER PLOTS (GALAXY-SPECIFIC)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   plotinfo = {dx: dx,$
               dy: dy,$
               map_r: map_r,$
               map_rkpc_ifs: map_rkpc_ifs,$
               kpc_per_as: kpc_per_as,$
               xran_kpc: xran_kpc,$
               yran_kpc: yran_kpc,$
               xarr_kpc: xarr_kpc,$
               yarr_kpc: yarr_kpc,$
               carr: carr,$
               center_nuclei: center_nuclei,$
               center_nuclei_kpc_x: center_nuclei_kpc_x,$
               center_nuclei_kpc_y: center_nuclei_kpc_y}

   if tag_exist(initmaps,'fcn_oplots') then begin
      if tag_exist(initmaps,'tags_oplots') then begin
         tags = initmaps.tags_oplots
         args_oplots = create_struct(tags[0], scope_varfetch(tags[0]))
         for i=1,n_elements(tags)-1 do $
            args_oplots = $
               struct_addtags(args_oplots,$
                              create_struct(tags[i],scope_varfetch(tags[i])))
         call_procedure,initmaps.fcn_oplots,initdat,initmaps,plotinfo,$
                        _extra=args_oplots
      endif else begin
         call_procedure,initmaps.fcn_oplots,initdat,initmaps,plotinfo
      endelse
   endif

badinput:

end
