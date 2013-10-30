
;cens syntax slitcens(cube,arcsec,angle,center,center

xmax = n_elements(red[*,0,0])
ymax = n_elements(red[0,*,0])

;plot box size of array
;define sides
lvx = [0,0] & lvy=[0,ymax]
rvx = [xmax,xmax] & rvy =[0,ymax]
bhx = [0,xmax] & bhy=[0,0]
uhx = [0,xmax] & uhy =[ymax,ymax]


plot,lvx,lvy,xran=[-10,xmax+10],yran=[-10,ymax+10],/xsty,/ysty
oplot,rvx,rvy
oplot,bhx,bhy
oplot,uhx,uhy

corners = cornersbox(5d,5d,ang,0,0)

plot,[corners[0,0],corners[1,0]],[corners[0,1],corners[1,1]],xran=[-17,17],yran=[-10,10],/xsty,/ysty,/nodata

for i = 0,36 do begin

ang = i*10d

corners = cornersbox(0.75d,0.75d,ang,0,0)

oplot,[corners[0,0],corners[1,0]],[corners[0,1],corners[1,1]]
oplot,[corners[1,0],corners[2,0]],[corners[1,1],corners[2,1]]
oplot,[corners[2,0],corners[3,0]],[corners[2,1],corners[3,1]]
oplot,[corners[3,0],corners[0,0]],[corners[3,1],corners[0,1]]


endfor



for i = 0,12 do begin

ang = i*30d

corners = cornersbox(1.25d,1.25d,ang,0,0)

oplot,[corners[0,0],corners[1,0]],[corners[0,1],corners[1,1]]
oplot,[corners[1,0],corners[2,0]],[corners[1,1],corners[2,1]]
oplot,[corners[2,0],corners[3,0]],[corners[2,1],corners[3,1]]
oplot,[corners[3,0],corners[0,0]],[corners[3,1],corners[0,1]]


endfor


for i = 0,10 do begin

ang = i*5d

corners = cornersbox(4d,2d,ang,0,0)

oplot,[corners[0,0],corners[1,0]],[corners[0,1],corners[1,1]]
oplot,[corners[1,0],corners[2,0]],[corners[1,1],corners[2,1]]
oplot,[corners[2,0],corners[3,0]],[corners[2,1],corners[3,1]]
oplot,[corners[3,0],corners[0,0]],[corners[3,1],corners[0,1]]


endfor

for i = 0,10 do begin

ang = i*5d

corners = cornersbox(2d,4d,ang,0,0)

oplot,[corners[0,0],corners[1,0]],[corners[0,1],corners[1,1]]
oplot,[corners[1,0],corners[2,0]],[corners[1,1],corners[2,1]]
oplot,[corners[2,0],corners[3,0]],[corners[2,1],corners[3,1]]
oplot,[corners[3,0],corners[0,0]],[corners[3,1],corners[0,1]]


endfor


for i = 0,10 do begin

ang = i*18d

corners = cornersbox(6d,4d,ang,0,0)

oplot,[corners[0,0],corners[1,0]],[corners[0,1],corners[1,1]]
oplot,[corners[1,0],corners[2,0]],[corners[1,1],corners[2,1]]
oplot,[corners[2,0],corners[3,0]],[corners[2,1],corners[3,1]]
oplot,[corners[3,0],corners[0,0]],[corners[3,1],corners[0,1]]


endfor


for i = 0,10 do begin

ang = i*18d

corners = cornersbox(10d,10d,ang,0,0)

oplot,[corners[0,0],corners[1,0]],[corners[0,1],corners[1,1]]
oplot,[corners[1,0],corners[2,0]],[corners[1,1],corners[2,1]]
oplot,[corners[2,0],corners[3,0]],[corners[2,1],corners[3,1]]
oplot,[corners[3,0],corners[0,0]],[corners[3,1],corners[0,1]]

endfor

for i = 0,18 do begin

ang = i*10d

corners = cornersbox(7d,7d,ang,0,0)

oplot,[corners[0,0],corners[1,0]],[corners[0,1],corners[1,1]]
oplot,[corners[1,0],corners[2,0]],[corners[1,1],corners[2,1]]
oplot,[corners[2,0],corners[3,0]],[corners[2,1],corners[3,1]]
oplot,[corners[3,0],corners[0,0]],[corners[3,1],corners[0,1]]

endfor

for i = 0,30 do begin

ang = i*9d

corners = cornersbox(14d,14d,ang,0,0)

oplot,[corners[0,0],corners[1,0]],[corners[0,1],corners[1,1]]
oplot,[corners[1,0],corners[2,0]],[corners[1,1],corners[2,1]]
oplot,[corners[2,0],corners[3,0]],[corners[2,1],corners[3,1]]
oplot,[corners[3,0],corners[0,0]],[corners[3,1],corners[0,1]]

endfor



END
