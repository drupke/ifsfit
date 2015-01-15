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
;      Size of subimage, in arcseconds.
;    ifsdims: in, required, type=dblarr(2)
;      Dimensions of IFS FOV in spaxels.
;    ifsps: in, required, type=double
;      Plate scale of IFS spaxels.
;    ifspa: in, required, type=dblarr(2,2)
;      Position angle of IFS field (east of north, in degrees).
;    ifsrefcoords: in, required, type=dblarr(2)
;      IFS coordinates of reference point for mutual centering.
;    hstrefcoords: in, required, type=dblarr(2)
;      HST coordinates of reference point for mutual centering.
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
function ifsf_hstsubim,image,subimsize,ifsdims,ifsps,ifspa,ifsrefcoords,$
                       hstrefcoords,scllim,hstps=hstps,ifsbounds=ifsbounds,$
                       fov=fov,sclargs=sclargs,noscl=noscl

;  HST platescale, in arcseconds
   if ~ keyword_set(hstps) then hstps = 0.05d

;  Rotation matrix to transform from sky coordinate system to instrument 
;  coordinate system:
;        |  cos(T) sin(T) | (x_sky) = (x_instr)
;        | -sin(T) cos(T) | (y_sky)   (y_instr)
;  where T is the PA (E of N) of the instrument FOV.
   ifspa_rad = ifspa*!DPi/180d
   sinifspa = sin(ifspa_rad)
   cosifspa = cos(ifspa_rad)
   rotmat = [[cosifspa,sinifspa],[-sinifspa,cosifspa]]

;  Center of IFS FOV in HST pixel coords. The offset specifies the X and Y 
;  offset of the reference point (e.g., the galaxy nucleus) from the IFS field 
;  center.
   ifsoffset = [double(ifsdims[0])/2d,double(ifsdims[1])/2d] - $
               [double(ifsrefcoords[0]),double(ifsrefcoords[1])]
   hst_center = hstrefcoords + (rotmat # ifsoffset)*$
                ifsps/hstps
 
; 
;  HST subimage coordinates
;
;  Case for creating subimage that's the size of the IFS FOV.
   if keyword_set(fov) then begin
;     The first subimage is a dummy one that's big enough to handle any rotation.
      twice_maxifsdims = dblarr(2)+2d*double(max(ifsdims))
      subimsizedumy = twice_maxifsdims*ifsps
;     Subimage size in HST pixels
      subimsize_pix=round(subimsizedumy/hstps)
;     Make sure subimage size in pixels is odd
      for j=0,1 do if (not subimsize_pix[j]) then subimsize_pix[j]+=1
;     Half of subimage size (minus a half pixel; i.e., rounded down)
      subimsize_pix_half=floor(subimsize_pix/2)
;     Find corners of subimage.
      subimcrddumy = intarr(4) ; coordinates of subimage corners
      subimcrddumy[0] = hst_center[0] - subimsize_pix_half[0]
      subimcrddumy[1] = hst_center[0] + subimsize_pix_half[0]
      subimcrddumy[2] = hst_center[1] - subimsize_pix_half[1]
      subimcrddumy[3] = hst_center[1] + subimsize_pix_half[1]
;     Subimage size for IFS FOV
      subimsize = double(ifsdims)*ifsps
;     For this case we redefine the center w.r.t. the dummy subimage.
      hst_center = [fix(double(subimcrddumy[1]-subimcrddumy[0])/2d +1d),$
                    fix(double(subimcrddumy[3]-subimcrddumy[2])/2d +1d)]
   endif
   
;  Subimage size in HST pixels
   subimsize_pix=round(subimsize/hstps)
;  Make sure subimage size in pixels is odd
   for j=0,1 do if (not subimsize_pix[j]) then subimsize_pix[j]+=1
;  Half of subimage size (minus a half pixel; i.e., rounded down)
   subimsize_pix_half=floor(subimsize_pix/2)
;  Find corners of subimage.
   subimcrd = intarr(4) ; coordinates of subimage corners
   subimcrd[0] = hst_center[0] - subimsize_pix_half[0]
   subimcrd[1] = hst_center[0] + subimsize_pix_half[0]
   subimcrd[2] = hst_center[1] - subimsize_pix_half[1]
   subimcrd[3] = hst_center[1] + subimsize_pix_half[1]

;  Create subimage and byte scale
   if keyword_set(fov) then use_subimcrd = subimcrddumy $
   else use_subimcrd = subimcrd
   hst_subim = temporary(image[use_subimcrd[0]:use_subimcrd[1],$
                               use_subimcrd[2]:use_subimcrd[3]])
   if ~ keyword_set(noscl) then $                       
      if keyword_set(sclargs) then $
         hst_subim = call_function('cgimgscl',hst_subim,minval=scllim[0],$
                                   max=scllim[1],_extra=sclargs) $
      else hst_subim = call_function('cgimgscl',hst_subim,minval=scllim[0],$
                                     max=scllim[1],stretch=5)

;  rotate and trim image
   if keyword_set(fov) then begin
      hst_subim = sshiftrotate(hst_subim,-(ifspa),xcen=hst_center[0],$
                               ycen=hst_center[1])
      hst_subim = hst_subim[subimcrd[0]:subimcrd[1],$
                            subimcrd[2]:subimcrd[3]]
   endif

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
      ifsbounds[0,*] = hst_center + rotmat # difs_1 - [subimcrd[0],subimcrd[2]]
      ifsbounds[1,*] = hst_center + rotmat # difs_2 - [subimcrd[0],subimcrd[2]]
      ifsbounds[2,*] = hst_center + rotmat # difs_3 - [subimcrd[0],subimcrd[2]]
      ifsbounds[3,*] = hst_center + rotmat # difs_4 - [subimcrd[0],subimcrd[2]]
   endif

   return,hst_subim
   
end
