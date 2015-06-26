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
;    comprange: in, optional, type=byte
;      Adjust ranges separately by component, rather than simultaneously for all
;      components.
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
pro ifsf_makemaps,initproc,comprange=comprange

   fwhm2sig = 2d*sqrt(2d*alog(2d))
   plotquantum = 2.5 ; in inches
   bad = 1d99
   c_kms = 299792.458d
   ncbdivmax = 7
   maxnadabscomp = 3
   maxnademcomp = 3

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
   hasratpar=0
   if tag_exist(initmaps,'center_axes') then $
      center_axes = initmaps.center_axes
   if tag_exist(initmaps,'center_nuclei') then $
      center_nuclei = initmaps.center_nuclei
   if tag_exist(initmaps,'argslinratmaps') then begin
      argslinratmaps = initmaps.argslinratmaps
      hasratpar=1
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

      size_tmp = size(linmaps[outlines[0]])
      dx = size_tmp[1]
      dy = size_tmp[2]
      if center_axes[0] eq -1 then center_axes = [double(dx)/2d,double(dy)/2d]
      if center_nuclei[0] eq -1 then center_nuclei = center_axes
   endif
   
;  Restore continuum parameters
   if tag_exist(initdat,'startempfile') then begin
      restore,file=initdat.outdir+initdat.label+'.cont.xdr'
      size_tmp = size(contcube.cont_fit_stel_tot)
      dx = size_tmp[1]
      dy = size_tmp[2]
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
   dohstbl=0
   dohstrd=0
   dohstsm=0
   dohstcol=0
   dohstcolsm=0
   if tag_exist(initmaps,'hst') AND tag_exist(initmaps,'hstbl') then begin
      dohstbl=1
      hstbl = readfits(initmaps.hstbl.file,header,/silent,/ext)
      hst_sm_ifsfov = dblarr(4,2)
      hst_big_ifsfov = dblarr(4,2)
      if tag_exist(initmaps.hstbl,'platescale') then $
         hstpsbl = initmaps.hstbl.platescale $
      else hstpsbl = 0.05d
      bhst_sm = ifsf_hstsubim(hstbl,[initmaps.hst.subim_sm,$
                              initmaps.hst.subim_sm],$
                              [dx,dy],initdat.platescale,$
                              initdat.positionangle,center_nuclei,$
                              initmaps.hst.refcoords,$
                              initmaps.hstbl.scllim,$
                              sclargs=initmaps.hstbl.sclargs_sm,$
                              ifsbounds=hst_sm_ifsfov,hstps=hstpsbl)
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
      hst_sm_ifsfov = dblarr(4,2)
      hst_big_ifsfov = dblarr(4,2)
      if tag_exist(initmaps.hstrd,'platescale') then $
         hstpsrd = initmaps.hstrd.platescale $
      else hstpsrd = 0.05d
      rhst_sm = ifsf_hstsubim(hstrd,[initmaps.hst.subim_sm,$
                              initmaps.hst.subim_sm],$
                              [dx,dy],initdat.platescale,$
                              initdat.positionangle,center_nuclei,$
                              initmaps.hst.refcoords,$
                              initmaps.hstrd.scllim,$
                              sclargs=initmaps.hstrd.sclargs_sm,$
                              ifsbounds=hst_sm_ifsfov,hstps=hstpsrd)
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
                               initmaps.hst.refcoords,$
                               initmaps.hstrd.scllim,$
                               sclargs=initmaps.hstrd.sclargs_fov,$
                               /fov,hstps=hstpsrd)
      rhst_fov_ns = ifsf_hstsubim(hstrd,[0,0],[dx,dy],initdat.platescale,$
                                  initdat.positionangle,center_nuclei,$
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
   if tag_exist(initmaps,'hst') AND tag_exist(initmaps,'hstcol') then begin
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
      chst_sm = ifsf_hstsubim(hstcol,[initmaps.hst.subim_sm,$
                              initmaps.hst.subim_sm],$
                              [dx,dy],initdat.platescale,$
                              initdat.positionangle,center_nuclei,$
                              initmaps.hst.refcoords,$
                              initmaps.hstcol.scllim,$
                              sclargs=initmaps.hstcol.sclargs,hstps=hstpsbl)
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
      cshst_sm = ifsf_hstsubim(hstcolsm,[initmaps.hst.subim_sm,$
                               initmaps.hst.subim_sm],$
                               [dx,dy],initdat.platescale,$
                               initdat.positionangle,center_nuclei,$
                               initmaps.hst.refcoords,$
                               initmaps.hstcolsm.scllim,$
                               sclargs=initmaps.hstcolsm.sclargs,hstps=hstpsbl)
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
; Compute some more things
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  coordinates in kpc
   xran_kpc = double([-(center_axes[0]-0.5),dx-(center_axes[0]-0.5)]) $
              * initdat.platescale * kpc_per_as
   yran_kpc = double([-(center_axes[1]-0.5),dy-(center_axes[1]-0.5)]) $
              * initdat.platescale * kpc_per_as
   center_nuclei_kpc_x = (center_nuclei[0,*] - center_axes[0]) $
                         * initdat.platescale * kpc_per_as  
   center_nuclei_kpc_y = (center_nuclei[1,*] - center_axes[1]) $
                         * initdat.platescale * kpc_per_as

;  Radii in kpc
;  GMOS FOV
   map_x = rebin(dindgen(dx)+1,dx,dy)
   map_y = rebin(transpose(dindgen(dy)+1),dx,dy)
   map_r = sqrt((map_x - center_axes[0])^2d + (map_y - center_axes[1])^2d)
   map_rkpc_ifs = map_r * initdat.platescale * kpc_per_as
;  HST FOV
   if (dohstrd OR dohstbl) then begin
      if dohstbl then size_subim = size(bhst_fov) $
      else size_subim = size(rhst_fov)
      map_x_hst = rebin(dindgen(size_subim[1])+1,size_subim[1],size_subim[2])
      map_y_hst = rebin(transpose(dindgen(size_subim[2])+1),$
                        size_subim[1],size_subim[2])
      center_axes_hst = center_axes * double(size_subim[1]/dx) - 1
      map_r_hst = sqrt((map_x_hst - center_axes_hst[0])^2d + $
                       (map_y_hst - center_axes_hst[1])^2d)
      if dohstbl then begin
         if dohstrd then begin
            if initmaps.hstbl.platescale ne initmaps.hstrd.platescale then begin
               print,'WARNING: HST blue and red plate scales differ;'
               print,'         using blue platescale for radius calculations.'
            endif
         endif
         map_rkpc_hst = map_r_hst * initmaps.hstbl.platescale * kpc_per_as
      endif else begin
         map_rkpc_hst = map_r_hst * initmaps.hstrd.platescale * kpc_per_as      
      endelse
      if dohstbl then $
         if tag_exist(initmaps.hstbl,'nucoffset') then begin
            map_r_bhst = $
               sqrt((map_x_hst - $
                     (center_axes_hst[0]+initmaps.hstbl.nucoffset[0]))^2d + $
                    (map_y_hst - $
                     (center_axes_hst[1]+initmaps.hstbl.nucoffset[1]))^2d)
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
      

;  Line ratios
   if hasratpar then linrats = ifsf_lineratios(linmaps,linelist)

;  Cumulative velocity distribution functions
   if ~ tag_exist(initmaps,'noemlinfit') then begin
      linspecmaps = hash()
      linspecpars = hash()
      foreach line,outlines do begin
;        if requested, apply E(B-V) only to those components that are specified
;        and tied to Halpha or to the same line as Halpha
         if tag_exist(initmaps,'applyebv') then begin
            if (initdat.linetie)[line] eq 'Halpha' OR $
               (initdat.linetie)[line] eq (initdat.linetie)['Halpha'] then $
                  doebv = $
                     linrats['ebv'] * $
                     rebin(reform(initmaps.applyebv,1,1,initdat.maxncomp),$
                           dx,dy,initdat.maxncomp)
         endif else doebv=0
         linspecmaps[line] = $
            ifsf_cmplinspecmaps(linmaps[line,*,*,*,4],$
                                linmaps[line,*,*,*,2],$
                                linmaps[line,*,*,*,3],$
                                initdat.maxncomp,linelist[line],$
                                initdat.zsys_gas,ebv=doebv)

         linspecpars[line] = ifsf_cmplinspecpars(linspecmaps[line])
      endforeach         
      linspecpars_tags = tag_names(linspecpars[outlines[0]])
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
   xarr_kpc[0] = xarr0_norm * (xran_kpc[1]-xran_kpc[0]) + xran_kpc[0]
   xarr_kpc[1] = xarr_kpc[0] + rarr_norm * (xran_kpc[1]-xran_kpc[0])*sinangarr
   xarr_kpc[2] = xarr_kpc[0] - rarr_norm * (xran_kpc[1]-xran_kpc[0])*cosangarr
   xarr_kpc[3] = xarr_kpc[0] + (rarr_norm+rlaboff_norm) * $
                               (xran_kpc[1]-xran_kpc[0])*sinangarr
   xarr_kpc[4] = xarr_kpc[0] - (rarr_norm+rlaboff_norm) * $
                               (xran_kpc[1]-xran_kpc[0])*cosangarr
   yarr_kpc[0] = yarr0_norm * (yran_kpc[1]-yran_kpc[0]) + yran_kpc[0]
   yarr_kpc[1] = yarr_kpc[0] + rarr_norm * (yran_kpc[1]-yran_kpc[0])*cosangarr
   yarr_kpc[2] = yarr_kpc[0] + rarr_norm * (yran_kpc[1]-yran_kpc[0])*sinangarr
   yarr_kpc[3] = yarr_kpc[0] + (rarr_norm+rlaboff_norm) * $
                               (yran_kpc[1]-yran_kpc[0])*cosangarr
   yarr_kpc[4] = yarr_kpc[0] + (rarr_norm+rlaboff_norm) * $
                               (yran_kpc[1]-yran_kpc[0])*sinangarr

   minyarr_kpc = min(yarr_kpc)
   if minyarr_kpc lt yran_kpc[0] then yarr_kpc -= minyarr_kpc - yran_kpc[0]
   maxxarr_kpc = max(xarr_kpc)
   if maxxarr_kpc gt xran_kpc[1] then xarr_kpc -= maxxarr_kpc - xran_kpc[1]

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
   if dohstrd AND dohstbl then begin
      cap1 = textoidl('F435W+F814W')
      cap2 = textoidl('F435W+F814W')
      cap3 = textoidl('F435W+F814W')
      cap4 = textoidl('IFS cont.')
   endif else if dohstrd then begin
      cap1 = textoidl('F814W')
      cap2 = textoidl('F814W')
      cap3 = textoidl('F814W')
      cap4 = textoidl('IFS cont.')
   endif else begin
      cap1 = textoidl('F435W')
      cap2 = textoidl('F435W')
      cap3 = textoidl('F435W')
      cap4 = textoidl('IFS cont.')
   endelse
;  arrays for positions for zoom box
   posbox1x = dblarr(2)
   posbox1y = dblarr(2)
   posbox2x = dblarr(2)
   posbox2y = dblarr(2)

   cgps_open,initdat.mapdir+initdat.label+'cont.eps',charsize=1,/encap,$
             /inches,xs=plotquantum*npx,ys=plotquantum*npy*aspectrat,/qui
   pos = cglayout([npx,npy],ixmar=[2d,2d],iymar=[2d,2d],$
                  oxmar=[1,0],oymar=[0,0],xgap=0,ygap=0)

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
      cgtext,size_subim[1]*0.05,size_subim[2]*0.9,textoidl('21\times21 kpc'),$
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
         mapscl = rebin(mapscl,3,size_subim[1]*10,size_subim[2]*10,/sample)
      endif else begin
         mapscl = bytarr(size_subim[1],size_subim[2])
         if dohstrd then mapscl = rhst_fov
         if dohstbl then mapscl = bhst_fov
         mapscl = rebin(mapscl,size_subim[1]*10,size_subim[2]*10,/sample)
      endelse
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos[*,i],opos=truepos,$
              noerase=i ne 0,missing_value=bad,missing_index=255,$
              missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
             /nodata,/noerase,title=cap2      
      ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
      cgtext,-2,1.8,'IFS FOV',/data,color='white'
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
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase,title=cap3
         ifsf_plotaxesnuc,xran_kpc,yran_kpc,center_nuclei_kpc_x,center_nuclei_kpc_y
         cgtext,-2,1.8,'IFS FOV, conv.',/data,color='white'
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
      cap1 = textoidl('F435W-F814W')
      cap2 = textoidl('F435W-F814W')
      cap3 = textoidl('F435W-F814W')
      cap4 = textoidl('IFS cont.')
      cbform='(D0.1)'

      cgps_open,initdat.mapdir+initdat.label+'color.eps',charsize=1,/encap,$
                /inches,xs=plotquantum*npx,ys=plotquantum*npy,/qui
      pos = cglayout([npx,npy],ixmar=[2,2],iymar=[2,2],oxmar=[0,0],oymar=[0,0],$
                     xgap=0,ygap=0)

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
      cgtext,size_subim[1]*0.05,size_subim[2]*0.9,textoidl('21\times21 kpc'),$
             color='black'
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

   if ~ tag_exist(initmaps,'noemlinfit') then begin

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

         cgps_open,initdat.mapdir+initdat.label+linelab+'.eps',charsize=1,/encap,$
                   /inches,xs=plotquantum*npx,ys=plotquantum*npy,/qui
         pos = cglayout([npx,npy],ixmar=[3,3],iymar=[3,3],oxmar=[0,0],oymar=[0,0],$
                        xgap=0,ygap=0)

;        loop through plot types
         for j=0,2 do begin

;           Set up colorbar labeling
            if j eq 0 then cbform = '(D0.1)' else cbform = '(I0)'


;           Set up ranges for all components at once
            if ~ keyword_set(comprange) then begin
        
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
               igd = where(map ne bad AND fluxmap ne 0 AND finite(map),ctgd)

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
                  if keyword_set(comprange) AND ~hasrange then begin
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plots of individual emission lines, components summed
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if ~ tag_exist(initmaps,'noemlinfit') then begin

;     Quantities to plot
      vtags = ['sig','vpk','v50','v84','v98']
      ftags = ['ftot','fpk','fv50','fv84','fv98']

;     Size of plot grid
      npx = n_elements(vtags)
      npy = 2

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

         cgps_open,initdat.mapdir+initdat.label+linelab+'_c.eps',charsize=1,/encap,$
                   /inches,xs=plotquantum*npx,ys=plotquantum*npy,/qui
         pos = cglayout([npx,npy],ixmar=[3,3],iymar=[3,3],oxmar=[0,0],oymar=[0,0],$
                        xgap=0,ygap=0)

;        loop through plot types
         for j=0,npx-1 do begin

;           FLUXES

            cbform = '(D0.1)' ; colorbar syntax
            iplot = j ; plot index
;           Get map and scale
            itag = where(strcmp(linspecpars_tags,ftags[j],/fold_case) eq 1)
            map = linspecpars[line].(itag)
            ibd = where(map eq bad AND ~ finite(map),ctbd)
            inan = where(~finite(map),ctnan)
            igd = where(map ne bad AND map ne 0 AND finite(map),ctgd)

            if ctgd gt 0 then begin
               
               if tag_exist(initmaps,'fluxfactor') then $
                  map[igd] *= initmaps.fluxfactor
                                 
               zran=[0,1]
               dzran = 1
               ncbdiv = 5
               zmax_flux = max(map[igd])
               
               map[igd] = map[igd]/zmax_flux
               if ctnan gt 0 then map[inan] = bad
               mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
                               min=zran[0],max=zran[1])
   
;              Plot image
               cgloadct,65,/reverse
               title='flux'
               title += ' ('+string(zmax_flux,format='(E0.2)')+')'
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

;           VELOCITIES
                          
            cbform = '(I0)' ; colorbar syntax
            iplot = npx+j ; plot index

;           Get map and scale
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
   
;              Plot image
               cgloadct,74,/reverse
               title=vtags[j]
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
 
         endfor
 
         cgps_close

      endforeach

   endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plots of line ratios
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if hasratpar then begin

;     Compute line ratios
      linrats = ifsf_lineratios(linmaps,linelist)

      keys = argslinratmaps->keys()

;     Loop through plot files
      foreach plt,keys do begin

;        Size of plot grid
         arrsiz = size(argslinratmaps[plt])
         npx = arrsiz[1]
         if arrsiz[0] gt 1 then begin
            npy = arrsiz[2]
            nplots = npx*npy
         endif else begin
            npy = 1
            nplots = npx
         endelse

         cgps_open,initdat.mapdir+initdat.label+plt+'.eps',charsize=1,/encap,$
            /inches,xs=plotquantum*npx,ys=plotquantum*npy,/qui
         pos = cglayout([npx,npy],ixmar=[3,3],iymar=[3,3],oxmar=[0,0],oymar=[0,0],$
            xgap=0,ygap=0)
         cbform = '(D0.2)'

;        Loop through plot panels
         for i=0,nplots-1 do begin

            tmpstr = strsplit(argslinratmaps[plt,i],'_',/extract)
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
               if keyword_set(comprange) AND ~hasrange then begin
                  zran = [min(map[igd]),max(map[igd])]
                  dzran = zran[1] - zran[0]
                  ncbdiv = ifsf_cbdiv(zran,0.5,7)
               endif

               title='c'+string(vcomp,format='(I0)')+' '
               if ptype eq 'n2ha' then title+=textoidl('[NII]/H\alpha')
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

      cgps_open,initdat.mapdir+initdat.label+'NaDempweq.eps',charsize=1,/encap,$
         /inches,xs=plotquantum*2,ys=plotquantum*3,/qui

      pos = cglayout([2,3],ixmar=[2,2],iymar=[2,2],oxmar=[0,0],oymar=[0,0],$
         xgap=0,ygap=0)
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

      cgps_close

   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD EQUIVALENT WIDTH (FITS)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if tag_exist(initdat,'donad') then begin

      cgps_open,initdat.mapdir+initdat.label+'NaDfitweq.eps',charsize=1,/encap,$
         /inches,xs=plotquantum*2,ys=plotquantum*3,/qui

      pos = cglayout([2,3],ixmar=[3,3],iymar=[3,3],oxmar=[0,0],oymar=[0,0],$
         xgap=0,ygap=0)
      cbform = '(I0)'

;
;     ABSORPTION
;
      map = nadfit.weqabs[*,*,0]
      maperrlo = nadfit.weqabserr[*,*,0]
      maperrhi = nadfit.weqabserr[*,*,1]
      maperrlo_emp = nadcube.weq[*,*,1]
      maperrhi_emp = nadcube.weq[*,*,1]
      igd = where(map gt 0d AND map ne bad)
      ibd = where(map eq 0d OR map eq bad)

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
      map = nadfit.weqem[*,*,0]
      maperrlo = nadfit.weqemerr[*,*,0]
      maperrhi = nadfit.weqemerr[*,*,1]
      maperrlo_emp = nadcube.weq[*,*,3]
      maperrhi_emp = nadcube.weq[*,*,3]
      igd = where(abs(map) gt 0d AND map ne bad)
      ibd = where(map eq 0d OR map eq bad)

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
      igd = where(abs(map) gt 0d AND map ne bad)
      ibd = where(map eq 0d OR map eq bad)

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
                     xgap=0,ygap=0,aspect=1)
      map = nadabsncomp
      ibd = where(map eq bad)
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

      cgps_open,initdat.mapdir+initdat.label+'NaDempvel.eps',charsize=1,/encap,$
         /inches,xs=plotquantum*3,ys=plotquantum*2,/qui

      pos = cglayout([3,2],ixmar=[3,3],iymar=[3,3],oxmar=[0,0],oymar=[0,0],$
         xgap=0,ygap=0)
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
         xgap=0,ygap=0)
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
         igd = where(map ne bad)
         ibd = where(map eq bad)
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
         igd = where(map ne bad)
         ibd = where(map eq bad)
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
         igd = where(map ne bad)
         ibd = where(map eq bad)
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
         igd = where(map ne bad)
         ibd = where(map eq bad)
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
         igd = where(map ne bad)
         ibd = where(map eq bad)
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
         igd = where(map ne bad)
         ibd = where(map eq bad)
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
