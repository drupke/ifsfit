; docformat = 'rst'
;
;+
;
; Extract subimages from HST images, with an option to extract the
; full IFS FOV and rotate the HST data so that the HST pixels are
; oriented the same way as the IFS spaxels. Byte scale the data using
; CGIMGSCL. Default stretch is ASINH; select different stretch using SCLARGS
; keyword.
;
; The procedure for matching HST and IFS subimages is to use a reference
; coordinate (specified as IFSREFCOORDS and HSTREFCOORDS). The exact center of 
; the IFS field is used to link these two systems. It assumes that HST pixel
; scale is an integer factor of the IFS spaxel scale. The key is that the IFS
; reference coordinate and field center remain constant during this procedure;
; the only thing that changes is the HST subimage position and rotation.
; 
; The procedure is as follows:
;
; 1. Extract HST subimage that is the size of the IFS FOV (in pixels) times 
;    some factor, plus a buffer on each side. The coordinates of the true IFS
;    field center are tracked in the subimage, assuming that the IFS and HST 
;    reference coordinates line up and given the position of the IFS reference
;    coordinate w.r.t. the field center. Because the HST subimage has to be
;    sampled discretely, the true IFS field center in the subimage doesn't 
;    necessarily lie at a pixel center or edge.
;    
; 2. Shift HST subimage a fractional number of pixels (using SSHIFTROTATE) so 
;    that the center of the HST subimage is truly the IFS field center. The 
;    IFS and HST reference coordinates should now line up.
;    
; 3. Rotate HST subimage to match orientation of IFS field. The rotation point
;    is the (integer) pixel nearest the IFS field center. The two reference
;    points no longer line up, because the HST subimage is rotated.
;    
; 4. Shift HST subimage again to re-match reference coordinates. Compute 
;    shift by rotating HST reference coordinate around rotation point.
; 
; NOTE: This procedure requires IDLUTILS because of the SSHIFTROTATE
; routine.
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
;                       replaced SSHIFTROTATE rotation with ROT; don't understand
;                       exact nature of rotation
;      2015jul27, DSNR, bug fix: rotation matrix was working CW, not CCW
;      2015sep22, DSNR, approximate fix for IFS bounding box
;    
; :Copyright:
;    Copyright (C) 2014-2015 David S. N. Rupke
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

;  Rotation matrix to transform from sky coordinate system to instrument 
;  coordinate system. Changing sin(T) to -sin(T) rotates clockwise back to 
;  sky coordinate system.
;        |  cos(T) sin(T) | (x_sky) = (x_instr)
;        | -sin(T) cos(T) | (y_sky)   (y_instr)
;  where T is the PA (E of N) of the instrument FOV.
   ifspa_rad = double(ifspa)*!DPi/180d
   sinifspa = sin(ifspa_rad)
   cosifspa = cos(ifspa_rad)
;  This matrix rotates counter-clockwise
   rotmatcc = [[cosifspa,sinifspa],[-sinifspa,cosifspa]]
;  This matrix rotates clockwise
   rotmatcw = [[cosifspa,-sinifspa],[sinifspa,cosifspa]]

;  The decimal X and Y offsets of the input IFS reference point from the IFS
;  field center, in IFS spaxels. The center of a spaxel is an integer, and a 
;  spaxel edge is a half-integer. Spaxel coordinates are single-offset.
   ifsoffset = [double(ifsdims[0])/2d,double(ifsdims[1])/2d] + 0.5d - $
               [double(ifsrefcoords[0]),double(ifsrefcoords[1])]
;  Decimal coordinates of IFS field center in the original HST image, in HST 
;  pixel coords, if the image were not rotated to match the IFS orientation. 
;  The center of a pixel is an integer, and a pixel edge is a half-integer. 
;  Pixel coordinates are single-offset.
   hst_center = double(hstrefcoords) + ifsoffset*ifsps/hstps
; 
;  HST subimage coordinates
;
;  Case for creating subimage that's the size of the IFS FOV.
   if keyword_set(fov) then begin
       buffer_ifs = max(ifsdims)
       buffer = buffer_ifs*ifsps/hstps
       subimsize = double(ifsdims)*ifsps
       rotang = ifspa
   endif else begin
;     Size chosen to deal with edge effects when shifting and resampling with
;     SSHIFTROTATE
      buffer=10
      rotang=0d
   endelse
   
;  Add buffer*2 pixels to subimage
   subimsize_pix=round(subimsize/hstps)+buffer*2
   
;  Half of subimage size in HST pixels. This will be an integer or half-
;  integer.
   subimsize_pix_half=double(subimsize_pix)/2d
;  Decimal coordinates of subimage corners, in coordinate system of 
;  original HST image, units of HST pixels, in zero-offset coordinates
;  so that IDL arrays can be addressed.
   subimcrd_dec = dblarr(4)
   subimcrd_dec[0] = hst_center[0]-1d - (subimsize_pix_half[0] -0.5d)
   subimcrd_dec[1] = hst_center[0]-1d + (subimsize_pix_half[0] -0.5d)
   subimcrd_dec[2] = hst_center[1]-1d - (subimsize_pix_half[1] -0.5d)
   subimcrd_dec[3] = hst_center[1]-1d + (subimsize_pix_half[1] -0.5d)
   subimcrd = round(subimcrd_dec)

;  Redefine the decimal coordinates of the IFS field center in the coordinate 
;  system of the subimage. The center is still expressed in single-offset 
;  coordinates, with pixel centers being integers. This expression also
;  accounts for the fact that the IFS field center may not be at the center
;  of an HST pixel.
   hst_center = [hst_center[0] - double(subimcrd[0]),$
                 hst_center[1] - double(subimcrd[2])]

;  Create subimage and byte scale
   hst_subim = temporary(image[subimcrd[0]:subimcrd[1],$
                               subimcrd[2]:subimcrd[3]])
   if ~ keyword_set(noscl) then $                       
      if keyword_set(sclargs) then $
         hst_subim = call_function('cgimgscl',hst_subim,minval=scllim[0],$
                                   max=scllim[1],_extra=sclargs) $
      else hst_subim = call_function('cgimgscl',hst_subim,minval=scllim[0],$
                                     max=scllim[1],stretch=5)

;
;  Shift to match reference coordinates.
;
;  Subimage center in single-offset coordinates.
   subim_cent = double(subimsize_pix) / 2d + 0.5d
;  Difference between subimage center and IFS field center
   subim_centdiff = subim_cent - hst_center
;  Shift subimage so that the center of the HST subimage is the IFS field center.
   hst_subim = sshiftrotate(hst_subim,0d,$
                            xshift=subim_centdiff[0],$
                            yshift=subim_centdiff[1])
;
;  Rotate the HST subimage so that it matches the IFS orientation.
;
;  Rotation is about the pixel center nearest the IFS field center, in zero-
;  offset coordinates. (Documentation of ROT doesn't specify that rotation center
;  is in zero-offset coordinates but I tested it.) Rotation can't be around a 
;  non-integer pixel center.
   rotcent = round(subim_cent)-1
   hst_subim = rot(hst_subim,-(rotang),1d,rotcent[0],rotcent[1],cubic=-0.5,$
                   /pivot)
;
;  Shift again to re-match reference coordinates, since rotation did not
;  occur around reference coordinate.
;  
;  Difference between IFS field center (subimage center) and rotation center, 
;  in HST pixels
   dxifs = subim_cent - (double(rotcent)+1d)
;  Transform x and y coordinates of IFS reference point from coordinate system
;  w.r.t. field center (subimage center) to coordinate system 
;  w.r.t. rotation center, in HST pixels
   ifsoffset_hst_rotcent = -(ifsoffset)*ifsps/hstps + dxifs
;  Rotate reference coordinate counterclockwise by PA
   ifsoffset_hst_rotcent_rot = rotmatcc # ifsoffset_hst_rotcent
;  Calculate difference between old and new position of reference coordinate
   refdiff = ifsoffset_hst_rotcent - ifsoffset_hst_rotcent_rot
;  Shift again so that reference coordinates line up.
   hst_subim = sshiftrotate(hst_subim,0d,$
                            xshift=refdiff[0],$
                            yshift=refdiff[1])
;  Trim subimage.
   hst_subim = hst_subim[buffer:subimsize_pix[0]-1-buffer,$
                         buffer:subimsize_pix[1]-1-buffer]

;  IFS boundaries in HST pixel coordinates, relative to the subimage.
   if keyword_set(ifsbounds) then begin
      dxifs_as = double(ifsdims[0])*ifsps/2d
      dyifs_as = double(ifsdims[1])*ifsps/2d
      dxifs = fix(dxifs_as/hstps)
      dyifs = fix(dyifs_as/hstps)
      difs_1 = [-dxifs,-dyifs]
      difs_2 = [ dxifs,-dyifs]
      difs_3 = [ dxifs, dyifs]
      difs_4 = [-dxifs, dyifs]
      ifsbounds[0,*] = (hst_center-buffer) + rotmatcw # difs_1 ; - [subimcrd[0],subimcrd[2]]
      ifsbounds[1,*] = (hst_center-buffer) + rotmatcw # difs_2 ; - [subimcrd[0],subimcrd[2]]
      ifsbounds[2,*] = (hst_center-buffer) + rotmatcw # difs_3 ; - [subimcrd[0],subimcrd[2]]
      ifsbounds[3,*] = (hst_center-buffer) + rotmatcw # difs_4 ; - [subimcrd[0],subimcrd[2]]
   endif

   return,hst_subim
   
end
