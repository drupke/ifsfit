; docformat = 'rst'
;
;+
; 
; Extract subimages from HST images, with an option to rotate the subimage
; and extract the IFS FOV. Byte scale the data using CGIMGSCL. 
; Default stretch is ASINH; select different stretch using SCLARGS
; keyword.
;
; In this routine there are three coordinates of interest:
;   1. the reference point, rref
;   2. the IFS field center, ric
;   3. the center of the HST subimage, rsc
;   
; There are several evolving coordinate systems:
;   0. the IFS image (_i)
;   1. the original HST image (_h)
;   2. the subimage (_hs)
;   3. the recentered subimage (_hsr)
;   4. the rotated, recentered subimage (_hsrr)
;   5. the re-recentered, rotated, recentered subimage (_hsrrr)
;   6. the shifted, re-re ... (_hsrrrs)
;
; The trick is to keep track of the coordinates of interest in whatever
; coordinate system using vector addition and rotation.
; 
; Because most coordinate systems are referenced to the lower-left pixxel, 
; zero-offset coordinates are used to make vector addition and rotation 
; straightforward.
;
; It is assumed that HST pixel scale is an integer factor of the IFS spaxel
; scale.
;
; The center of a spaxel is an integer, and a spaxel edge is a half-integer.
;
; Coordinate system rotation is opposite to what might seem obvious because the
; HST image data are rotated opposite to the IFS PA, which then makes the 
; IFS in a row-column orientation w.r.t. the HST image.
; 
; The procedure requires IDLUTILS because of the SSHIFTROTATE routine.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;   Image array.
;   
; :Params:
;    image: in, required, type=dblarr(2)
;      HST image, as read by READFITS.
;    subimsize: in, required, type=dblarr(2)
;      Size of output subimage, in arcseconds. This parameter is ignored if
;      the FOV keyword is set.
;    ifsdims: in, required, type=dblarr(2)
;      Dimensions of IFS FOV in spaxels.
;    ifsps: in, required, type=double
;      Plate scale of IFS spaxels.
;    ifspa: in, required, type=dblarr(2,2)
;      Position angle of IFS field (east of north, in degrees).
;    ifsrefcoords: in, required, type=dblarr(2)
;      IFS coordinates of reference point for mutual centering, in single-offset
;      spaxel coordinates. Spaxel centers are integers, and edges are half-integers.
;      This location and HSTREFCOORDS in principle represent exactly the 
;      same physical location.
;    hstrefcoords: in, required, type=dblarr(2)
;      HST coordinates of reference point for mutual centering, in single-offset
;      pixel coordinates. Pixel centers are integers, and edges are half-integers.
;      This location and IFSREFCOORDS in principle represent exactly the 
;      same physical location.
;    scllim: in, required, type=dblarr(2)
;      Intensity limits for byte scale routine.
;    
; :Keywords:
;    hstps: in, optional, type=double, default=0.05
;      Plate scale of input image, in arcseconds.
;    ifsbounds: out, optional, type=dblarr(4,2)
;      Returns the coordinates of the corners of the IFS FOV in the
;      HST subimage coordinate system.
;    fov: in, optional, type=byte
;      Select the entire IFS FOV as the subimage. If selected, HST data
;      will be rotated to align its pixels with the IFS spaxels.
;    noscl: in, optional, type=byte
;      Skip byte scaling.
;    sclargs: in, optional, type=structure
;      Optional arguments to byte scaling function.
;    
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
;      2014may20, DSNR, created
;      2014dec15, DSNR, removed STRETCH keyword; use SCLARGS to change stretch
;      2015jul08, DSNR, re-written to get centering and rotation correct
;                       replaced SSHIFTROTATE rotation with ROT; don't 
;                       understand exact nature of rotation
;      2015jul27, DSNR, bug fix: rotation matrix was working CW, not CCW
;      2015sep22, DSNR, approximate fix for IFS bounding box
;      2016may13, DSNR, routine re-write because of issues with IFSBOUNDS
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
function ifsf_hstsubim,image,subimsize,ifsdims,ifsps,ifspa,ifsrefcoords,$
                       hstrefcoords,scllim,hstps=hstps,ifsbounds=ifsbounds,$
                       fov=fov,sclargs=sclargs,noscl=noscl
                       
;  HST platescale, in arcseconds
   if ~ keyword_set(hstps) then hstps = 0.05d

;  Rotation matrix to transform from one coordinate system to another, in 
;  degrees T counterclockwise.
;  Changing sin(T) to -sin(T) rotates the other direction.
;        |  cos(T) -sin(T) | (x1) = (x2)
;        |  sin(T)  cos(T) | (y1)   (y2)
   ifspa_rad = double(ifspa)*!DPi/180d
   sinifspa = sin(ifspa_rad)
   cosifspa = cos(ifspa_rad)
;  This matrix rotates counter-clockwise
   rotmatccw = [[cosifspa,sinifspa],[-sinifspa,cosifspa]]
;  This matrix rotates clockwise
   rotmatcw = [[cosifspa,-sinifspa],[sinifspa,cosifspa]]

;  Half dimensions of IFS field of view in IFS pixels
   dxifs_half = double(ifsdims[0])/2d
   dyifs_half = double(ifsdims[1])/2d
;  Half dimensions of IFS field of view in HST pixels
   dxifs_half_hst = dxifs_half*ifsps/hstps
   dyifs_half_hst = dyifs_half*ifsps/hstps

;1. Define coordinates

;  Convert reference coordinates from 1-offset to 0-offset
   rref_h = hstrefcoords - 1d
   rref_i = ifsrefcoords - 1d

;  The IFS field center in IFS coordinates
   ric_i = [double(ifsdims[0])/2d,double(ifsdims[1])/2d]-0.5d
;  The IFS field center in HST coordinates
   ric_h = ifsps/hstps * rotmatccw # (ric_i - rref_i) + rref_h

;2. Define and extract subimage

;  Case for creating subimage that's the size of the IFS FOV.
   if keyword_set(fov) then begin
      buffer = max(ifsdims)*ifsps/hstps
      subimsize = double(ifsdims)*ifsps
   endif else begin
      buffer=0
   endelse
   
;  Add buffer*2 pixels to subimage
   subimsize_pix=round(subimsize/hstps)+buffer*2
;  Make sure subimsize is odd
   if not fix(subimsize_pix[0]) then subimsize_pix[0]++
   if not fix(subimsize_pix[1]) then subimsize_pix[1]++

;  Round ric_h to nearest integer to determine subimage center
   rsc_h = double(round(ric_h))

;  Half of subimage size in HST pixels. This will be a half-integer.
   subimsize_pix_half=double(subimsize_pix)/2d

;  Coordinates of subimage corners, in coordinate system of 
;  original HST image, units of HST pixels, zero-offset coordinates.
   subimcrd_h = dblarr(4)
   subimcrd_h[0] = rsc_h[0] - (subimsize_pix_half[0] -0.5d)
   subimcrd_h[1] = rsc_h[0] + (subimsize_pix_half[0] -0.5d)
   subimcrd_h[2] = rsc_h[1] - (subimsize_pix_half[1] -0.5d)
   subimcrd_h[3] = rsc_h[1] + (subimsize_pix_half[1] -0.5d)
;  to make sure result is an integer and not v. close to an integer
   subimcrd = round(subimcrd_h)

;  Create subimage and byte scale
   hst_subim = temporary(image[subimcrd[0]:subimcrd[1],$
                               subimcrd[2]:subimcrd[3]])
   if ~ keyword_set(noscl) then $                       
      if keyword_set(sclargs) then $
         hst_subim = call_function('cgimgscl',hst_subim,minval=scllim[0],$
                                   max=scllim[1],_extra=sclargs) $
      else hst_subim = call_function('cgimgscl',hst_subim,minval=scllim[0],$
                                     max=scllim[1],stretch=5)

;  IC and SC in new coordinate system
   ric_hs = ric_h - double([subimcrd[0],subimcrd[2]])
   rsc_hs = rsc_h - double([subimcrd[0],subimcrd[2]])

;  If the HST image is being rotated / clipped to match the IFS field of view
   if keyword_set(fov) then begin

;3. Switch to coordinate system where rsc = [0,0]

      ric_hsr = ric_hs - rsc_hs
;
;4. Rotate the HST subimage so that it matches the IFS orientation.
;
;     (Documentation of ROT doesn't specify that rotation center
;     is in zero-offset coordinates but I tested it.) Rotation can't be 
;     around a non-integer pixel center. The pivot keyword is unnecessary
;     since we've defined rsc_hs to be the center of the subimage.
      hst_subim = rot(hst_subim,ifspa,1d,rsc_hs[0],rsc_hs[1],cubic=-0.5,$
                      /pivot)

;     Rotated coordinates
      ric_hsrr = rotmatcw # ric_hsr

;5. Switch to coordinate system where the lower left corner of the subimage
;   is [0,0]

      ric_hsrrr = rsc_hs + ric_hsrr

;6. Shift so IC is on a pixel center or edge. If, for a given dimension, the 
;   number of HST pixels is even, then the IC must lie on the edge of two pixels.
;   If it's odd, it must lie in the center of a pixel.
;
;   Recall that odd integers are true and even integers are false.
;   
;   A positive shift in SSHIFTROTATE means that the underlying image data 
;   moves to lower coordinate values.
;   
;   SSHIFTROTATE uses sinc interpolation
      shift = dblarr(2)
      ifsdims_hstpix = ifsdims * ifsps/hstps
;     If there is an odd number of pixels, we want to move IC to a pixel center.
;     E.g., if IC = 1.1, then it rounds to 1.0.
;     If there is an even number of pixels, we want to move IC to a pixel
;     edge, and we'll shift by an extra 0.5 pixels.
      shift = -(ric_hsrrr - double(round(ric_hsrrr)))
      if not round(ifsdims_hstpix[0]) then shift[0]-=0.5d
      if not round(ifsdims_hstpix[1]) then shift[1]-=0.5d
      hst_subim = sshiftrotate(hst_subim,0d,xshift=shift[0],yshift=shift[1])
      ric_hsrrrs = ric_hsrrr + shift
         
;7. Extract IFS field from subimage.

;     Coordinates of IFS field
      subimcrd_ifs = dblarr(4)
      subimcrd_ifs[0] = ric_hsrrrs[0] - (dxifs_half_hst-0.5d)
      subimcrd_ifs[1] = ric_hsrrrs[0] + (dxifs_half_hst-0.5d)
      subimcrd_ifs[2] = ric_hsrrrs[1] - (dyifs_half_hst-0.5d)
      subimcrd_ifs[3] = ric_hsrrrs[1] + (dyifs_half_hst-0.5d)
;     to make sure result is an integer and not v. close to an integer
      subimcrd_ifs = round(subimcrd_ifs)

      hst_subim = hst_subim[subimcrd_ifs[0]:subimcrd_ifs[1],$
                            subimcrd_ifs[2]:subimcrd_ifs[3]]

   endif

;  IFS boundaries in HST pixel coordinates, relative to the subimage.
   if keyword_set(ifsbounds) then begin
      ifsbounds[0,*] = $
         ifsps/hstps*rotmatccw # ([0d,0d]-0.5d - rref_i) + rref_h - $
         double([subimcrd[0],subimcrd[2]])
       ifsbounds[1,*] = $
         ifsps/hstps*rotmatccw # ([double(ifsdims[0]),0d]-0.5d - rref_i) + rref_h - $
         double([subimcrd[0],subimcrd[2]])
       ifsbounds[2,*] = $
         ifsps/hstps*rotmatccw # (double(ifsdims)-0.5d - rref_i) + rref_h - $
         double([subimcrd[0],subimcrd[2]])
       ifsbounds[3,*] = $
         ifsps/hstps*rotmatccw # ([0d,double(ifsdims[1])]-0.5d - rref_i) + rref_h - $
         double([subimcrd[0],subimcrd[2]])
   endif

   return,hst_subim
   
end
