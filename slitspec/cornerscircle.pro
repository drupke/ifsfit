;+
; NAME:
;     cornerscircle
;
; PURPOSE:
;     get corners of circle of diameter d center xcen,ycen
;
; EXPLANATION:
;
; CALLING SEQUENCE
;     
;
; INPUTS:
;     d - diameter of circle
;     xcen - center x
;     ycen - center y
;
; OUTPUT:
;     corners of rectangle cw from upper left corner, (e.g. when at 0 angle)
;     array "corners"
; METHOD:
;
; REVISION HISTORY:
;     10Nov26  Jeff Rich created
;-

function cornerscircle,d,xcen,ycen

d = double(d)
r = d/2d
xcen = double(xcen)
ycen = double(ycen)

;number of corners for circle
npts = 6250

corners = dblarr(npts,2)

for i=0,npts-1 do begin

    ang = ((360d/npts) * i) * (!pi/180d)

    x = r * cos(ang) + xcen
    y = r * sin(ang) + ycen

    corners[i,0] = x
    corners[i,1] = y

endfor

return, corners

END
