;+
; NAME:
;     slitcens
;
; PURPOSE:
;    Part of script for generating pseudo-images from ifu cubes
;    Given a cube, slit element length, PA, xcen and ycen  
;    Generates an array of x & y values for centers of 
;    slit elements
;
; EXPLANATION:
;
; CALLING SEQUENCE
;     
;
; INPUTS:
;     cube - WiFeS data cube 
;     d - slit element length
;     ang-slit angle clockwise from up (w of n if PA=0)
;     xcen
;     ycen
;
; OUTPUT:
;     array of slit center locations
;
; METHOD:
;
; REVISION HISTORY:
;     10Nov25  Jeff Rich created
;-

function slitcens,cube,d,ang,xcen,ycen

radang = (ang/180d)*!pi

xmax=n_elements(cube[*,0,0])
ymax=n_elements(cube[0,*,0])

if ang ne 90 AND ang ne 270 then begin

    yo = sqrt((d^2)/(tan(radang)^2 + 1))
    xo = yo*tan(radang)

endif else begin

    yo = 0
    xo = d

endelse

;get num slit elements to left and right of center

if ang ge 0 AND ang lt 90 then begin
nrx = abs((xmax-xcen)/xo)
nry = abs((ymax-ycen)/yo)
nlx = abs(xcen/xo)
nly = abs(ycen/yo)
endif
if ang ge 90 AND ang lt 180 then begin
nlx = abs((xmax-xcen)/xo)
nry = abs((ymax-ycen)/yo)
nrx = abs(xcen/xo)
nly = abs(ycen/yo)
endif
if ang ge 180 AND ang lt 270 then begin
nrx = abs((xmax-xcen)/xo)
nry = abs((ymax-ycen)/yo)
nlx = abs(xcen/xo)
nly = abs(ycen/yo)
endif
if ang ge 270 AND ang lt 360 then begin
nlx = abs((xmax-xcen)/xo)
nry = abs((ymax-ycen)/yo)
nrx = abs(xcen/xo)
nly = abs(ycen/yo)
endif


if ~finite(nlx) then nlx = 1
if ~finite(nly) then nly = 1
if ~finite(nrx) then nrx = 1
if ~finite(nry) then nry = 1


if nlx lt nly then nl = nlx else nl = nly
if nrx lt nry then nr = nrx else nr = nry

;print,nrx,nry,nr
;print, nlx,nly,nl


nr = round(nr)
nl = round(nl)

n = nr+nl

cens = dblarr(n+1,2)


i=0

for x=n-nr,n do begin
cens[x,0] = xcen + (i*xo)
cens[x,1] = ycen + (i*yo)
i+=1
endfor

i=0

for x=n-nr,0,-1 do begin
cens[x,0] = xcen - (i*xo)
cens[x,1] = ycen - (i*yo)
i+=1
endfor

;stop


return, cens

END
