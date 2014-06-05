; docformat = 'rst'
;
;+
;
; This procedure makes maps of various quantities.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    IDL save file (.xdr)
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
;      2014jun02, DSNR, updated to allow use without previous emission-line fit
;                       with IFSF
;      2014jun04, DSNR, updated to plot velocities using cumulative velocity
;                       distributions
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
pro ifsf_makemaps,initproc,comprange=comprange

   fwhm2sig = 2d*sqrt(2d*alog(2d))
   plotquantum = 2.5 ; in inches
   bad = 1d99
   c_kms = 299792.458d
   ncbdivmax = 7

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Load initialization parameters and line data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  Get fit initialization
   initmaps={dumy: 1}
   initdat=call_function(initproc,initmaps=initmaps)

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

;  Get NaD equivalent widths (empirical)
   if tag_exist(initdat,'donad') then begin
      restore,file=initdat.outdir+initdat.label+'.nadspec.xdr'
      if tag_exist(initmaps,'noemlinfit') then begin
         size_tmp = size(nadcube.weq)
         dx = size_tmp[1]
         dy = size_tmp[2]
         if center_axes[0] eq -1 then center_axes = [double(dx)/2d,double(dy)/2d]
         if center_nuclei[0] eq -1 then center_nuclei = center_axes
      endif
   endif
   
   if ~ tag_exist(initdat,'donad') AND $
      ~ tag_exist(initmaps,'noemlinfit') then begin
      print,'IFSF_MAKEMAPS: No emission line or absorption line data specified.'
      print,'               Aborting.'
      goto,badinput
   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Load and process continuum data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  Data cube
   if not tag_exist(initdat,'datext') then datext=1 else datext=initdat.datext
   if not tag_exist(initdat,'varext') then varext=2 else varext=initdat.varext
   if not tag_exist(initdat,'dqext') then dqext=3 else dqext=initdat.dqext
   ctcube = ifsf_readcube(initdat.infile,/quiet,oned=oned,$
                          datext=datext,varext=varext,dqext=dqext)

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
      bhst_sm = ifsf_hstsubim(hstbl,[initmaps.hst.subim_sm,$
                              initmaps.hst.subim_sm],$
                              [dx,dy],initdat.platescale,$
                              initdat.positionangle,center_nuclei,$
                              initmaps.hst.refcoords,$
                              initmaps.hstbl.scllim,$
                              sclargs=initmaps.hstbl.sclargs,$
                              ifsbounds=hst_sm_ifsfov)
      bhst_big = ifsf_hstsubim(hstbl,[initmaps.hst.subim_big,$
                               initmaps.hst.subim_big],$
                               [dx,dy],initdat.platescale,$
                               initdat.positionangle,center_nuclei,$
                               initmaps.hst.refcoords,$
                               initmaps.hstbl.scllim,$
                               sclargs=initmaps.hstbl.sclargs,$
                               ifsbounds=hst_big_ifsfov)
      bhst_fov = ifsf_hstsubim(hstbl,[0,0],[dx,dy],initdat.platescale,$
                               initdat.positionangle,center_nuclei,$
                               initmaps.hst.refcoords,$
                               initmaps.hstbl.scllim,$
                               sclargs=initmaps.hstbl.sclargs,$
                               /fov)
      if tag_exist(initmaps,'hstblsm') then begin
         dohstsm=1
         hstblsm = filter_image(hstbl,fwhm=initmaps.hst.smoothfwhm,/all)
         bhst_fov_sm = ifsf_hstsubim(hstblsm,[0,0],[dx,dy],$
                                     initdat.platescale,$
                                     initdat.positionangle,center_nuclei,$
                                     initmaps.hst.refcoords,$
                                     initmaps.hstblsm.scllim,$
                                     stretch=initmaps.hstblsm.stretch,$
                                     sclargs=initmaps.hstblsm.sclargs,$
                                     /fov)
         
      endif
   endif      
   if tag_exist(initmaps,'hst') AND tag_exist(initmaps,'hstrd') then begin
      dohstrd=1
      hstrd = readfits(initmaps.hstrd.file,header,/silent,/ext)
      hst_sm_ifsfov = dblarr(4,2)
      hst_big_ifsfov = dblarr(4,2)
      rhst_sm = ifsf_hstsubim(hstrd,[initmaps.hst.subim_sm,$
                              initmaps.hst.subim_sm],$
                              [dx,dy],initdat.platescale,$
                              initdat.positionangle,center_nuclei,$
                              initmaps.hst.refcoords,$
                              initmaps.hstrd.scllim,$
                              sclargs=initmaps.hstrd.sclargs,$
                              ifsbounds=hst_sm_ifsfov)
      rhst_big = ifsf_hstsubim(hstrd,[initmaps.hst.subim_big,$
                               initmaps.hst.subim_big],$
                               [dx,dy],initdat.platescale,$
                               initdat.positionangle,center_nuclei,$
                               initmaps.hst.refcoords,$
                               initmaps.hstrd.scllim,$
                               sclargs=initmaps.hstrd.sclargs,$
                               ifsbounds=hst_big_ifsfov)
      rhst_fov = ifsf_hstsubim(hstrd,[0,0],[dx,dy],initdat.platescale,$
                               initdat.positionangle,center_nuclei,$
                               initmaps.hst.refcoords,$
                               initmaps.hstrd.scllim,$
                               sclargs=initmaps.hstrd.sclargs,$
                               /fov)
      if tag_exist(initmaps,'hstrdsm') then begin
         dohstsm=1
         hstrdsm = filter_image(hstrd,fwhm=initmaps.hst.smoothfwhm,/all)
         rhst_fov_sm = ifsf_hstsubim(hstrdsm,[0,0],[dx,dy],$
                                     initdat.platescale,$
                                     initdat.positionangle,center_nuclei,$
                                     initmaps.hst.refcoords,$
                                     initmaps.hstrdsm.scllim,$
                                     stretch=initmaps.hstrdsm.stretch,$
                                     sclargs=initmaps.hstrdsm.sclargs,$
                                     /fov)
         
      endif
   endif      
   if tag_exist(initmaps,'hst') AND tag_exist(initmaps,'hstcol') then begin
      dohstcol=1
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
      hstcol = alog10(hstrd/hstbl)
      hstcol[ibdcol] = 1d99
;     Extract and scale
      chst_sm = ifsf_hstsubim(hstcol,[initmaps.hst.subim_sm,$
                              initmaps.hst.subim_sm],$
                              [dx,dy],initdat.platescale,$
                              initdat.positionangle,center_nuclei,$
                              initmaps.hst.refcoords,$
                              initmaps.hstcol.scllim,$
                              stretch=initmaps.hstcol.stretch,$
                              sclargs=initmaps.hstcol.sclargs)
      chst_big = ifsf_hstsubim(hstcol,[initmaps.hst.subim_big,$
                               initmaps.hst.subim_big],$
                               [dx,dy],initdat.platescale,$
                               initdat.positionangle,center_nuclei,$
                               initmaps.hst.refcoords,$
                               initmaps.hstcol.scllim,$
                               stretch=initmaps.hstcol.stretch,$
                               sclargs=initmaps.hstcol.sclargs)
      chst_fov = ifsf_hstsubim(hstcol,[0,0],[dx,dy],initdat.platescale,$
                               initdat.positionangle,center_nuclei,$
                               initmaps.hst.refcoords,$
                               initmaps.hstcol.scllim,$
                               stretch=initmaps.hstcol.stretch,$
                               sclargs=initmaps.hstcol.sclargs,$
                               /fov)
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
      hstcolsm = alog10(hstrdsm/hstblsm)
      hstcolsm[ibdcol] = 1d99
;     Extract and scale
      cshst_sm = ifsf_hstsubim(hstcolsm,[initmaps.hst.subim_sm,$
                               initmaps.hst.subim_sm],$
                               [dx,dy],initdat.platescale,$
                               initdat.positionangle,center_nuclei,$
                               initmaps.hst.refcoords,$
                               initmaps.hstcolsm.scllim,$
                               stretch=initmaps.hstcolsm.stretch,$
                               sclargs=initmaps.hstcolsm.sclargs)
      cshst_big = ifsf_hstsubim(hstcolsm,[initmaps.hst.subim_big,$
                                initmaps.hst.subim_big],$
                                [dx,dy],initdat.platescale,$
                                initdat.positionangle,center_nuclei,$
                                initmaps.hst.refcoords,$
                                initmaps.hstcolsm.scllim,$
                                stretch=initmaps.hstcolsm.stretch,$
                                sclargs=initmaps.hstcolsm.sclargs)
      cshst_fov = ifsf_hstsubim(hstcolsm,[0,0],[dx,dy],initdat.platescale,$
                                initdat.positionangle,center_nuclei,$
                                initmaps.hst.refcoords,$
                                initmaps.hstcolsm.scllim,$
                                stretch=initmaps.hstcolsm.stretch,$
                                sclargs=initmaps.hstcolsm.sclargs,$
                                /fov)
;     Extract unscaled color image and convert to same pixel scale as IFS data
      cshst_fov_ns = ifsf_hstsubim(hstcolsm,[0,0],[dx,dy],initdat.platescale,$
                                   initdat.positionangle,center_nuclei,$
                                   initmaps.hst.refcoords,[0,0],/noscl,/fov)
      cshst_fov_rb = congrid(cshst_fov_ns,dx,dy,/interp)
      size_fov_hst = size(cshst_fov)
      cshst_fov_rb /= size_fov_hst[1]*size_fov_hst[2]/(dx*dy)
   endif
   hstrd=0
   hstbl=0
   hstcol=0
   hstrdsm=0
   hstblsm=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Compute some things
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  Luminosity and angular size distances
   ldist = lumdist(initdat.zsys_gas,H0=73,Omega_m=0.27,Lambda0=0.73,/silent)
   kpc_per_as = ldist/(1+initdat.zsys_gas^2)*1000d/206265d

;  coordinates in kpc
   xran_kpc = double([-(center_axes[0]-0.5),dx-(center_axes[0]-0.5)]) $
              * initdat.platescale * kpc_per_as
   yran_kpc = double([-(center_axes[1]-0.5),dy-(center_axes[1]-0.5)]) $
              * initdat.platescale * kpc_per_as
   center_nuclei_kpc_x = (center_nuclei[0,*] - center_axes[0]) $
                         * initdat.platescale * kpc_per_as  
   center_nuclei_kpc_y = (center_nuclei[1,*] - center_axes[1]) $
                         * initdat.platescale * kpc_per_as

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
   cap1 = textoidl('F435W+F814W')
   cap2 = textoidl('F435W+F814W')
   cap3 = textoidl('IFS cont.')

   cgps_open,initdat.mapdir+initdat.label+'cont.eps',charsize=1,/encap,$
             /inches,xs=plotquantum*npx,ys=plotquantum*npy,/qui
   pos = cglayout([npx,npy],ixmar=[3,3],iymar=[3,3],oxmar=[0,0],oymar=[0,0],$
                   xgap=0,ygap=0)

   if (dohstrd OR dohstbl) then begin
      i = 0
;     HST continuum, large scale
      if dohstbl then size_subim = size(bhst_big) $
      else size_subim = size(rhst_big)
      mapscl = bytarr(3,size_subim[1],size_subim[2])
      if dohstrd then mapscl[0,*,*] = rhst_big
      if dohstbl then mapscl[2,*,*] = bhst_big
      if dohstrd AND dohstbl then $
         mapscl[1,*,*] = byte((double(rhst_big)+double(bhst_big))/2d)
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos[*,i],opos=truepos,$
              noerase=i ne 0,missing_value=bad,missing_index=255,$
              missing_color='white'
      cgplot,[0],xsty=5,ysty=5,xran=[0,size_subim[1]],$
             yran=[0,size_subim[2]],position=truepos,$
             /nodata,/noerase,title=cap1
      cgoplot,[hst_big_ifsfov[*,0],hst_big_ifsfov[0,0]],$
              [hst_big_ifsfov[*,1],hst_big_ifsfov[0,1]],color='Red'


      i = 1
;     HST continuum, IFS FOV
      if dohstbl then size_subim = size(bhst_fov) $
      else size_subim = size(rhst_fov)
      mapscl = bytarr(3,size_subim[1],size_subim[2])
      if dohstrd then mapscl[0,*,*] = rhst_fov
      if dohstbl then mapscl[2,*,*] = bhst_fov
      if dohstrd AND dohstbl then $
         mapscl[1,*,*] = byte((double(rhst_fov)+double(bhst_fov))/2d)
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos[*,i],opos=truepos,$
              noerase=i ne 0,missing_value=bad,missing_index=255,$
              missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
             /nodata,/noerase,title=cap2
      cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
      cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
      cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
      cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
      cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1

      i = 2
;     smoothed HST continuum, IFS FOV
      if dohstsm then begin
         mapscl = bytarr(3,size_subim[1],size_subim[2])
         if dohstrd then mapscl[0,*,*] = rhst_fov_sm
         if dohstbl then mapscl[2,*,*] = bhst_fov_sm
         if dohstrd AND dohstbl then $
            mapscl[1,*,*] = byte((double(rhst_fov_sm)+double(bhst_fov_sm))/2d)
         cgloadct,65,/reverse
         cgimage,mapscl,/keep,pos=pos[*,i],opos=truepos,$
                 noerase=i ne 0,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase,title=cap2
         cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
         cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
         cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
         cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
         cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
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
             /nodata,/noerase,title=cap3
      cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
      cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
      cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
      cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
      cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
   endif

   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Continuum color plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if dohstcol then begin
      npx = 2
      npy = 2
   endif else begin
      npx = 1
      npy = 1
   endelse
   cap1 = textoidl('F435W+F814W')
   cap2 = textoidl('F435W+F814W')
   cap3 = textoidl('IFS cont.')

   cgps_open,initdat.mapdir+initdat.label+'color.eps',charsize=1,/encap,$
             /inches,xs=plotquantum*npx,ys=plotquantum*npy,/qui
   pos = cglayout([npx,npy],ixmar=[3,3],iymar=[3,3],oxmar=[0,0],oymar=[0,0],$
                   xgap=0,ygap=0)

   if dohstcol then begin
      i = 0
      cgloadct,65
      cgimage,chst_big,/keep,pos=pos[*,i],opos=truepos,$
              noerase=i ne 0,missing_value=bad,missing_index=255,$
              missing_color='white'
      cgplot,[0],xsty=5,ysty=5,xran=[0,size_subim[1]],$
             yran=[0,size_subim[2]],position=truepos,$
             /nodata,/noerase,title=cap1
      cgoplot,[hst_big_ifsfov[*,0],hst_big_ifsfov[0,0]],$
              [hst_big_ifsfov[*,1],hst_big_ifsfov[0,1]],color='Red'


      i = 1
      cgloadct,65
      cgimage,chst_fov,/keep,pos=pos[*,i],opos=truepos,$
              noerase=i ne 0,missing_value=bad,missing_index=255,$
              missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
             /nodata,/noerase,title=cap2
      cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
      cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
      cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
      cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
      cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1

      i = 2
;     smoothed HST continuum, IFS FOV
      if dohstcolsm then begin
         cgloadct,65
         cgimage,cshst_fov,/keep,pos=pos[*,i],opos=truepos,$
                 noerase=i ne 0,missing_value=bad,missing_index=255,$
                 missing_color='white'
         cgplot,[0],xsty=5,ysty=5,position=truepos,$
                /nodata,/noerase,title=cap2
         cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
         cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
         cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
         cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
         cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
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
             /nodata,/noerase,title=cap3
      cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
      cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
      cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
      cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
      cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
   endif

   cgps_close

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
               ibd = where(map eq bad AND ~ finite(map),ctbd)
               inan = where(~finite(map),ctnan)
               igd = where(map ne bad AND map ne 0 AND finite(map),ctgd)

               if ctgd gt 0 then begin
            
                  if j eq 1 then begin
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
;                 Plot axes in kpc
                  cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
                  cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
                  cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
                  cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
;                 Plot nuclei
                  cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
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
;              Plot axes in kpc
               cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
               cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
               cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
               cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
;              Plot nuclei
               cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
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
;              Plot axes in kpc
               cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
               cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
               cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
               cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
;              Plot nuclei
               cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
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
               cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
               cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
               cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
               cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
               cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
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
                     /iso
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

      nadabsweq_snrthresh=4d
      nademweq_snrthresh=4d

      cgps_open,initdat.mapdir+initdat.label+'NaDweq.eps',charsize=1,/encap,$
         /inches,xs=plotquantum*2,ys=plotquantum*3,/qui

      pos = cglayout([2,3],ixmar=[3,3],iymar=[3,3],oxmar=[0,0],oymar=[0,0],$
         xgap=0,ygap=0)
      cbform = '(I0)'

;
;     ABSORPTION
;
      map = nadcube.weq[*,*,0]
      igd = where(map ge nadabsweq_snrthresh*nadcube.weq[*,*,1] AND $
                  map gt 0d AND map ne bad)
      ibd = where(map lt nadabsweq_snrthresh*nadcube.weq[*,*,1] OR $
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

;     replace bad points with 0
      map[ibd] = 0d

      mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
                      min=zran[0],max=zran[1])
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos[*,0],opos=truepos,$
         noerase=i ne 0,missing_value=bad,missing_index=255,$
         missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
         /nodata,/noerase,title='W$\down eq$(NaD abs, $\angstrom$)'
      cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
      cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
      cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
      cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
      cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
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
         /nodata,/noerase,title='W/$\delta$W(NaD abs)'
         cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
      cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
      cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
      cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
      cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6

;
;     EMISSION
;
      map = nadcube.weq[*,*,2]
      igd = where(abs(map) ge nademweq_snrthresh*nadcube.weq[*,*,3] AND $
                  abs(map) gt 0d AND map ne bad)
      ibd = where(abs(map) lt nademweq_snrthresh*nadcube.weq[*,*,3] OR $
                  map eq 0d AND map ne bad)

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

;     replace bad points with 0
      map[ibd] = 0d

      mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
         min=zran[0],max=zran[1])
      cgloadct,65
      cgimage,mapscl,/keep,pos=pos[*,2],opos=truepos,$
         noerase=i ne 0,missing_value=bad,missing_index=255,$
         missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
         /nodata,/noerase,title='W$\down eq$(NaD em, $\angstrom$)'
         cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
      cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
      cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
      cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
      cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
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
         /nodata,/noerase,title='W/$\delta$W(NaD em)'
      cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
      cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
      cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
      cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
      cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6



;     Flux
      cbform = '(D0.2)'
      map = nadcube.emflux[*,*,0]
      igd = where(abs(nadcube.weq[*,*,2]) ge $
                  nademweq_snrthresh*nadcube.weq[*,*,3] AND $
                  abs(map) gt 0d AND map ne bad)
      ibd = where(abs(nadcube.weq[*,*,2]) lt $
                  nademweq_snrthresh*nadcube.weq[*,*,3] OR $
                  map eq 0d AND map ne bad)

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
         zmax_flux = zran[1]
         zran=[0,1]
         dzran = 1
         ncbdiv = 5
      endif
      
;     replace bad points with 0
      map[ibd] = 0d
;     normalize good points
      map[igd] /= zmax_flux

      mapscl = bytscl(rebin(map,dx*20,dy*20,/sample),$
         min=zran[0],max=zran[1])
      cgloadct,65,/reverse
      cgimage,mapscl,/keep,pos=pos[*,4],opos=truepos,$
         noerase=i ne 0,missing_value=bad,missing_index=255,$
         missing_color='white'
      cgplot,[0],xsty=5,ysty=5,position=truepos,$
         /nodata,/noerase,title='Flux (NaD em)'
      cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
      cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
      cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
      cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
      cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6


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
      igd = where(nadcube.weq[*,*,0] ge nadabsweq_snrthresh*nadcube.weq[*,*,1] $
                  AND map gt 0d AND map ne bad)
      ibd = where(nadcube.weq[*,*,0] lt nadabsweq_snrthresh*nadcube.weq[*,*,1] $
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
         cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
      cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
      cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
      cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
      cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6

      map = nadcube.vel[*,*,1]
      igd = where(nadcube.weq[*,*,0] ge nadabsweq_snrthresh*nadcube.weq[*,*,1] $
                  AND map ne bad)
      ibd = where(nadcube.weq[*,*,0] lt nadabsweq_snrthresh*nadcube.weq[*,*,1] $
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
         cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
      cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
      cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
      cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
      cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6


      map = nadcube.vel[*,*,2]
      igd = where(nadcube.weq[*,*,0] ge nadabsweq_snrthresh*nadcube.weq[*,*,1] $
                  AND map ne bad)
      ibd = where(nadcube.weq[*,*,0] lt nadabsweq_snrthresh*nadcube.weq[*,*,1] $
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
         cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
      cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
      cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
      cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
      cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
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
                  nademweq_snrthresh*nadcube.weq[*,*,3] $
                  AND map gt 0d AND map ne bad)
      ibd = where(abs(nadcube.weq[*,*,2]) lt $
                  nademweq_snrthresh*nadcube.weq[*,*,3] $
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
         cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
      cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
      cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
      cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
      cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6

      map = nadcube.vel[*,*,4]
      igd = where(abs(nadcube.weq[*,*,2]) ge $
                  nademweq_snrthresh*nadcube.weq[*,*,3] $
                  AND map ne bad)
      ibd = where(abs(nadcube.weq[*,*,2]) lt $
                  nademweq_snrthresh*nadcube.weq[*,*,3] $
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
         cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
      cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
      cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
      cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
      cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6

      map = nadcube.vel[*,*,5]
      igd = where(abs(nadcube.weq[*,*,2]) ge $
                  nademweq_snrthresh*nadcube.weq[*,*,3] $
                  AND map ne bad)
      ibd = where(abs(nadcube.weq[*,*,2]) lt $
                  nademweq_snrthresh*nadcube.weq[*,*,3] $
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
         cgaxis,xaxis=0,xran=xran_kpc,/xsty,/save
      cgaxis,xaxis=1,xran=xran_kpc,xtickn=replicate(' ',60),/xsty
      cgaxis,yaxis=0,yran=yran_kpc,/ysty,/save
      cgaxis,yaxis=1,yran=yran_kpc,ytickn=replicate(' ',60),/ysty
      cgoplot,center_nuclei_kpc_x,center_nuclei_kpc_y,psym=1
      cbpos=[truepos[2],truepos[1],truepos[2]+0.01,truepos[3]]
      ticknames = string(dindgen(ncbdiv+1)*dzran/double(ncbdiv) - $
         (dzran - zran[1]),format=cbform)
      cgcolorbar,position=cbpos,divisions=ncbdiv,$
         ticknames=ticknames,/ver,/right,charsize=0.6

      cgps_close

   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NaD: Color vs. equivalent width and emission-line flux
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if tag_exist(initdat,'donad') AND $
      tag_exist(initmaps,'hst') AND $
      tag_exist(initmaps,'hstcolsm') then begin

      cgps_open,initdat.mapdir+initdat.label+'NaDweq_v_col.eps',charsize=0.75,$
         /encap,/inches,xs=plotquantum*2,ys=plotquantum*1,/qui
      pos = cglayout([2,1],ixmar=[3,3],iymar=[3,3],oxmar=[0,0],oymar=[8,4],$
         xgap=0)

      cmap = cshst_fov_rb
      
;     ABSORPTION
      map = nadcube.weq[*,*,0]
      igd = where(nadcube.weq[*,*,0] ge nadabsweq_snrthresh*nadcube.weq[*,*,1] $
                  AND map gt 0d AND map ne bad)
      ibd = where(nadcube.weq[*,*,0] lt nadabsweq_snrthresh*nadcube.weq[*,*,1] $
                  OR map eq 0d OR map eq bad)
      map[ibd] = bad
      hasrange = 0
      if hasrangefile then begin
         ithisline = where(rangeline eq 'NaDabs' AND $
            rangequant eq 'weq',ctthisline)
         if ctthisline eq 1 then begin
            yran = [rangelo[ithisline],rangehi[ithisline]]
            hasrange = 1
         endif
      endif
      if ~hasrange then yran = [min(map[igd]),max(map[igd])]

      xran = [min(cmap[igd]),max(cmap[igd])]
      cgplot,cmap,map,/xsty,/ysty,xran=xran,yran=yran,psym=16,pos=pos[*,0],$
         xtit='log(red/blue)',ytit='W$\down eq$(NaD abs, $\angstrom$)'

;     EMISSION
      map = nadcube.emflux[*,*,0]
      igd = where(abs(nadcube.weq[*,*,2]) ge $
                  nademweq_snrthresh*nadcube.weq[*,*,3] AND $
                  abs(map) gt 0d AND map ne bad)
      ibd = where(abs(nadcube.weq[*,*,2]) lt $
                  nademweq_snrthresh*nadcube.weq[*,*,3] OR $
                  map eq 0d AND map ne bad)
      map[ibd] = bad
      hasrange = 0
      if hasrangefile then begin
         ithisline = where(rangeline eq 'NaDabs' AND $
            rangequant eq 'weq',ctthisline)
         if ctthisline eq 1 then begin
            yran = [rangelo[ithisline],rangehi[ithisline]]
            hasrange = 1
         endif
      endif
      if ~hasrange then yran = [min(map[igd]),max(map[igd])]

      xran = [min(cmap[igd]),max(cmap[igd])]
      cgplot,cmap,map,/xsty,/ysty,xran=xran,yran=yran,psym=16,pos=pos[*,1],$
         /noerase,xtit='log(red/blue)',ytit='Flux (NaD em)'

      cgps_close

   endif

badinput:

end
