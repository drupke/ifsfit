name = 'IRASF10257-4339'

cuber = '/Users/jrich/wifes/cubes/'+name+'_R7_res_t_x2.fits'
cubeb = '/Users/jrich/wifes/cubes/'+name+'_B3_res_t_x2.fits'



if n_elements(red) eq 0 then red = wifes_readcube(cuber,wave_red,err_red)
if n_elements(blue) eq 0 then blue = wifes_readcube(cubeb,wave_blue,err_blue)

imor = mrdfits(cuber,0,hdror)
imr  = mrdfits(cuber,1,hdrimr)
errr = mrdfits(cuber,2,hdrerrr)
imob = mrdfits(cubeb,0,hdrob)
imb  = mrdfits(cubeb,1,hdrimb)
errb = mrdfits(cubeb,2,hdrerrb)

xmax = n_elements(red[*,0,0])
ymax = n_elements(red[0,*,0])
ang=90d
xcen=25.5
ycen=20.5
len=12.8d
wid=12.8d

plotname= name+'_'+strcompress(long(ang),/remove_all)+'_'+strcompress(long(xcen),/remove_all)+'_'+strcompress(long(ycen),/remove_all)+'_slit.ps'


nwvred =n_elements(red[0,0,*]) 
nwvblue=n_elements(blue[0,0,*]) 

cens = slitcens(red,len,ang,xcen,ycen)

;initial plot
set_plot,'ps'
device,filename=plotname,/color

!p.thick=2

plot,[0,0],[0,ymax],xran=[-2,xmax+2],yran=[-2,ymax+2],/xsty,/ysty
oplot,[0,xmax],[ymax,ymax]
oplot,[xmax,xmax],[ymax,0]
oplot,[0,xmax],[0,0]

num = n_elements(cens[*,0])

redslit = dblarr(num,nwvred)
blueslit= dblarr(num,nwvblue)
err_redslit = dblarr(num,nwvred)
err_blueslit= dblarr(num,nwvblue)

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
        if xl ge 0 AND xl lt xmax AND yl ge 0 AND yl lt ymax then begin 
            redslit[n,p] += red[list[q,0],list[q,1],p]*list[q,2]
            err_redslit[n,p] += err_red[list[q,0],list[q,1],p]*list[q,2]
        endif
        endfor
endfor
;blue slit
for p=0,nwvblue-1 do begin
    for q=0,n_elements(list[*,0])-1 do begin
        ;pixel list coords/frac
        xl = list[q,0]
        yl = list[q,1]
        fl = list[q,2]
        if xl ge 0 AND xl lt xmax AND yl ge 0 AND yl lt ymax then begin
            blueslit[n,p] += blue[list[q,0],list[q,1],p]*list[q,2]
            err_blueslit[n,p] += err_blue[list[q,0],list[q,1],p]*list[q,2]
        endif
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

device,/close
set_plot,'x'


;change header info

dlamr = sxpar(hdrimr,'cd3_3')
dlamb = sxpar(hdrimb,'cd3_3')
lamor = sxpar(hdrimr,'CRVAL3')
lamob = sxpar(hdrimb,'CRVAL3')
pixor = sxpar(hdrimr,'CRPIX3')
pixob = sxpar(hdrimb,'CRPIX3')

sxaddpar,hdrimr,'NAXIS1',num
sxaddpar,hdrerrr,'NAXIS1',num
sxaddpar,hdrimb,'NAXIS1',num
sxaddpar,hdrerrb,'NAXIS1',num

sxaddpar,hdrimr,'CD1_1',len
sxaddpar,hdrerrr,'CD1_1',len
sxaddpar,hdrimb,'CD1_1',len
sxaddpar,hdrerrb,'CD1_1',len

sxaddpar,hdrimr,'CD2_2',dlamr
sxaddpar,hdrerrr,'CD2_2',dlamr
sxaddpar,hdrimb,'CD2_2',dlamb
sxaddpar,hdrerrb,'CD2_2',dlamb

sxdelpar,hdrimr,'CD3_3'
sxdelpar,hdrerrr,'CD3_3'
sxdelpar,hdrimb,'CD3_3'
sxdelpar,hdrerrb,'CD3_3'

sxaddpar,hdrimr,'CRVAL2',lamor
sxaddpar,hdrerrr,'CRVAL2',lamor
sxaddpar,hdrimb,'CRVAL2',lamob
sxaddpar,hdrerrb,'CRVAL2',lamob

sxdelpar,hdrimr,'CRVAL3'
sxdelpar,hdrerrr,'CRVAL3'
sxdelpar,hdrimb,'CRVAL3'
sxdelpar,hdrerrb,'CRVAL3'

sxaddpar,hdrimr,'CRPIX2',pixor
sxaddpar,hdrerrr,'CRPIX2',pixor
sxaddpar,hdrimb,'CRPIX2',pixob
sxaddpar,hdrerrb,'CRPIX2',pixob

sxdelpar,hdrimr,'CRPIX3'
sxdelpar,hdrerrr,'CRPIX3'
sxdelpar,hdrimb,'CRPIX3'
sxdelpar,hdrerrb,'CRPIX3'

sxdelpar,hdrimr,'NAXIS3'
sxdelpar,hdrerrr,'NAXIS3'
sxdelpar,hdrimb,'NAXIS3'
sxdelpar,hdrerrb,'NAXIS3'

sxaddpar,hdrimr,'CDELT2',dlamr
sxaddpar,hdrerrr,'CDELT2',dlamr
sxaddpar,hdrimb,'CDELT2',dlamb
sxaddpar,hdrerrb,'CDELT2',dlamb

sxdelpar,hdrimr,'CDELT3'
sxdelpar,hdrerrr,'CDELT3'
sxdelpar,hdrimb,'CDELT3'
sxdelpar,hdrerrb,'CDELT3'

redname = name+'_'+strcompress(long(ang),/remove_all)+'_'+strcompress(long(xcen),/remove_all)+'_'+strcompress(long(ycen),/remove_all)+'_r.fits'
bluename= name+'_'+strcompress(long(ang),/remove_all)+'_'+strcompress(long(xcen),/remove_all)+'_'+strcompress(long(ycen),/remove_all)+'_b.fits'


mwrfits,[0,0,0],redname,hdror
mwrfits,redslit,redname,hdrimr
mwrfits,err_redslit,redname,hdrerrr

mwrfits,[0,0,0],bluename,hdrob
mwrfits,blueslit,bluename,hdrimb
mwrfits,err_blueslit,bluename,hdrerrb

!p.thick=1

END
