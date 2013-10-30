if n_elements(red) eq 0 then red = wifes_readcube('/Users/jrich/wifes/cubes/IRASF10257-4339_R7_res_t_x2.fits',wave_red,err_red)
if n_elements(blue) eq 0 then blue = wifes_readcube('/Users/jrich/wifes/cubes/IRASF10257-4339_B3_res_t_x2.fits',wave_blue,err_blue)

xmax = n_elements(red[*,0,0])
ymax = n_elements(red[0,*,0])
ang=337d
xcen=20d
ycen=20d
len=2d
wid=3d

nwvred =n_elements(red[0,0,*]) 
nwvblue=n_elements(blue[0,0,*]) 

cens = slitcens(red,len,ang,xcen,ycen)

;initial plot
plot,[0,0],[0,ymax],xran=[-2,xmax+2],yran=[-2,ymax+2],/xsty,/ysty
oplot,[0,xmax],[ymax,ymax]
oplot,[xmax,xmax],[ymax,0]
oplot,[0,xmax],[0,0]

num = n_elements(cens[*,0])

redslit = dblarr(num,nwvred)
blueslit= dblarr(num,nwvblue)

for n=0,num-1 do begin

;corner of box for slit n
corners = cornersbox(len,wid,ang,cens[n,0],cens[n,1])

oplot,[corners[0,0],corners[1,0]],[corners[0,1],corners[1,1]]
oplot,[corners[1,0],corners[2,0]],[corners[1,1],corners[2,1]]
oplot,[corners[2,0],corners[3,0]],[corners[2,1],corners[3,1]]
oplot,[corners[3,0],corners[0,0]],[corners[3,1],corners[0,1]]

list = overlap(corners)

;make slit element
;red slit
for p=0,nwvred-1 do begin
    for q=0,n_elements(list[*,0])-1 do begin
        ;pixel list coords/frac
        xl = list[q,0]
        yl = list[q,1]
        fl = list[q,2]
        if xl ge 0 AND xl lt xmax AND yl ge 0 AND yl lt ymax then redslit[n,p] += red[list[q,0],list[q,1],p]*list[q,2]
    endfor
endfor
;blue slit
for p=0,nwvblue-1 do begin
    for q=0,n_elements(list[*,0])-1 do begin
        ;pixel list coords/frac
        xl = list[q,0]
        yl = list[q,1]
        fl = list[q,2]
        if xl ge 0 AND xl lt xmax AND yl ge 0 AND yl lt ymax then blueslit[n,p] += blue[list[q,0],list[q,1],p]*list[q,2]
    endfor
endfor
          

;if i eq 0 then listall = list else listall = [listall,list]

for i=0,n_elements(list[*,0])-1 do begin
    
    xi = list[i,0] & xf = list[i,0]+1
    yi = list[i,1] & yf = list[i,1]+1
    oplot,[xi,xi],[yi,yf],color=fsc_color('red')
    oplot,[xi,xf],[yf,yf],color=fsc_color('red')
    oplot,[xf,xf],[yf,yi],color=fsc_color('red')
    oplot,[xi,xf],[yi,yi],color=fsc_color('red')
endfor
    
endfor

for n=0,num-1 do begin

;corner of box for slit n
corners = cornersbox(len,wid,ang,cens[n,0],cens[n,1])

list = overlap(corners)

for i=0, n_elements(list[*,0])-1 do begin    
    if list[i,2] eq 1 then begin
    xi = list[i,0] & xf = list[i,0]+1
    yi = list[i,1] & yf = list[i,1]+1
    oplot,[xi,xi],[yi,yf],color=fsc_color('blue')
    oplot,[xi,xf],[yf,yf],color=fsc_color('blue')
    oplot,[xf,xf],[yf,yi],color=fsc_color('blue')
    oplot,[xi,xf],[yi,yi],color=fsc_color('blue')
    endif
    
endfor
endfor

END
