;+
; NAME:
;     cornersbox
;
; PURPOSE:
;     get corners of rectangle of length l width w angle cw center xcen,ycen
;
; EXPLANATION:
;
; CALLING SEQUENCE
;     
;
; INPUTS:
;     l - length of rectangle (clock hand)
;     w - width of rectangle (perp)
;     ang-angle of clock hand (cw from up)
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

function cornersbox,l,w,ang,xcen,ycen

l = double(l)
w = double(w)
ang = double(ang)
xcen = double(xcen)
ycen = double(ycen)

radang = (ang/180d)*!pi

phi1 = (!pi/2) - radang + atan(w/l)
phi2 = (!pi/2) - radang - atan(w/l)

r = sqrt((l/2d)^2 + (w/2d)^2)

x1 = r*cos(phi1)
y1 = r*sin(phi1)
x2 = r*cos(phi2)
y2 = r*sin(phi2)

;print,x1,y1
;print,x2,y2

corners = dblarr(4,2)

corners[0,0] = xcen+x1
corners[1,0] = xcen+x2
corners[2,0] = xcen-x1
corners[3,0] = xcen-x2

corners[0,1] = ycen+y1
corners[1,1] = ycen+y2
corners[2,1] = ycen-y1
corners[3,1] = ycen-y2

;for i=0,3 do print,corners[i,0],corners[i,1]


;plot,[corners[0,0],corners[1,0]],[corners[0,1],corners[1,1]],xran=[-17,17],yran=[-10,10],/xsty,/ysty
;oplot,[corners[1,0],corners[2,0]],[corners[1,1],corners[2,1]]
;oplot,[corners[2,0],corners[3,0]],[corners[2,1],corners[3,1]]
;oplot,[corners[3,0],corners[0,0]],[corners[3,1],corners[0,1]]

;stop

return, corners

END
