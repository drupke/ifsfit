; docformat = 'rst'
;
;+
;
; Take an arbitrarily oriented cross section of an IFS map. From
; D. Fanning, of course!
;
; 
; :Categories:
;    IFSF
;
; :Returns:
;
; :Params:
;    map: in, required, type=dblarr(dx,dy)
;      IFS map.
;    center: in, required, type=dblarr(2)
;      Center for cross section, in zero-offset coordinates.
;    length: in, required, type=double
;      Length of cross section, in spaxels.
;    angle: in, required, type=double
;      Angle for cross section, CCW from +y axis (i.e., E of N), in degrees.
; 
; :Keywords:
;    bad: in, optional, type=double, default=1d99
;      Value for bad values in map.
;    nearest: in, optional, type=byte
;      Default is binlinear interpolation. Select this for nearest-neighbor.
;    ends: out, optional, type=dblarr(2,2)
;      2x2 array with spaxel endpoints for slit
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
;      2013may14, DSNR, created
;      2015nov24, DSNR, documented
;      2015dec08, DSNR, finished re-writing; added bad data masking and radial
;                       computation
;      2016jul14, DSNR, renamed to match convention
;    
; :Copyright:
;    Copyright (C) 2015--2016 David S. N. Rupke
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
function ifsf_xsec,map,center,length,angle,nearest=nearest,bad=bad,$
                   ends=ends

   if not keyword_set(bad) then bad = 1d99

;  Create bad data mask. 0 = good, 1 = bad
   mapflg = map*0d
   ibd = where(map eq bad,ctbad)
   if ctbad gt 0 then mapflg[ibd]=1d

   halflength=length/2d
   sinangle = sin(angle*!DPi/180d)
   cosangle = cos(angle*!DPi/180d)
   xends = center[0]+[halflength*sinangle,$
                      -halflength*sinangle]
   yends = center[1]+[-halflength*cosangle,$
                       halflength*cosangle]

;  Suppose your two endpoints were given by (x1,y1) and (x2,y2). To
;  calculate an image profile along the line, you must know the image
;  values at discrete points along the line. There are, of course, an
;  infinity of points along the line from which you can calculate the
;  image value, but in practice I like to use the larger of the
;  endpoint differences. My IDL code looks like this:
   nPoints = ABS(xends[1]-xends[0]+1) > ABS(yends[1]-yends[0]+1)

;  Next, you must construct two vectors containing the X and Y
;  locations, respectively, of these points. This is easily done in
;  IDL, like this:
   xloc = xends[0] + (xends[1] - xends[0]) * dindgen(nPoints) / (nPoints - 1)
   yloc = yends[0] + (yends[1] - yends[0]) * dindgen(nPoints) / (nPoints - 1)

;  Finally, the profile values are calculated by interpolating the
;  image values at these locations along the line with the Interpolate
;  command. (These will be bilinearly interpolated values, by default.
;  If you prefer nearest neighbor interpolation, you can use the Round
;  function in IDL.)
   if keyword_set(nearest) then begin
      xsec = Interpolate(map, Round(xloc), Round(yloc))
;     Masking out bad values by interpolating bad data mask.
      if ctbad gt 0 then begin
         xsecflg = Interpolate(mapflg, Round(xloc), Round(yloc))
         ibad = where(xsecflg gt 0)
         xsec[ibad]=bad
      endif
   endif else begin
      xsec = Interpolate(map, xloc, yloc)
      if ctbad gt 0 then begin
         xsecflg = Interpolate(mapflg, xloc, yloc)
         ibad = where(xsecflg gt 0)
         xsec[ibad]=bad
      endif
   endelse

   ends=[[xends],[yends]]

;  Radial location of xsec, in spaxels
   xloc_c = xloc - center[0]
   yloc_c = yloc - center[1]
   ipos = where(xloc_c ge 0)
   ineg = where(xloc_c lt 0)
   rloc = dblarr(n_elements(xloc_c))
   rloc[ipos] = sqrt(xloc_c[ipos]^2d + yloc_c[ipos]^2d)
   rloc[ineg] = -sqrt(xloc_c[ineg]^2d + yloc_c[ineg]^2d)

   return,[[rloc],[xsec]]

end
