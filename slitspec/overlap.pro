;+
; NAME:
;     overlap
;
; PURPOSE:
;     generate list of overlapping grid spaxels & percent overlap
;     given a particular box with corners 
; EXPLANATION:
;
; CALLING SEQUENCE
;     
;
; INPUTS:
;     ;cube - cube for overlap;assuming 1 arcsec square grid
;     corners - corners of slit formated from cornersbox.pro
;     
;
; OUTPUT:
;     list of spaxels & cooresponding overlap percentage/weight
;     list[0,*] = [x,y,percent]
; METHOD:
;
; REVISION HISTORY:
;     10Nov26  Jeff Rich created
;-

function overlap,corners

;get mins/maxes from corners
xmin = min(corners[*,0])
xmax = max(corners[*,0])
ymin = min(corners[*,1])
ymax = max(corners[*,1])
;print,xmin,ymin,xmax,ymax
;if xmin lt 0 then xmin = floor(xmin) else xmin = ceil(xmin)
;if xmax lt 0 then xmax = floor(xmax) else xmax = ceil(xmax)
;if ymin lt 0 then ymin = floor(ymin) else ymin = ceil(ymin)
;if ymax lt 0 then ymax = floor(ymax) else ymax = ceil(ymax)
xmin = floor(xmin)
ymin = floor(ymin)
xmax = ceil(xmax)
ymax = ceil(ymax)
;ordered arrays of corners points so I can ignore orientation
xord = corners[sort(corners[*,0]),0]
yord = corners[sort(corners[*,0]),0]


i=0

;get length of list array
for x=xmin,xmax do begin

    for y=ymin,ymax do begin
        xs = corners[*,0]
        ys = corners[*,1]
        
        polyclip,x,y,xs,ys

        if n_elements(xs) ne 1 then i+=1

    endfor
endfor

;print,i
;create list
if i ne 0 then list = dblarr(i,3) else begin 
    list = dblarr(1,3)
    list[0,*] = -1
endelse

i=0

for x=xmin,xmax do begin

    for y=ymin,ymax do begin


        
        xs = corners[*,0]
        ys = corners[*,1]
        
        polyclip,x,y,xs,ys

        if n_elements(xs) ne 1 then begin
            
            list[i,0] = x
            list[i,1] = y
            list[i,2] = double(poly_area(xs,ys))
                                    
            i+=1

 ;           stop

        endif
        
    endfor

endfor

return, list

END
