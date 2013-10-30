;read in data cube

infile = '/Users/jrich/wifes/cubes/IRASF10257-4339'

red = wifes_readcube(infile+'_R7_res_t_x2.fits',wave_red,red_err)
blue = wifes_readcube(infile+'_B3_res_t_x2.fits',wave_blue,blue_err)

;center
;brightest point in median red contin
xcen = 26
ycen = 21

;make a slit?
;slit center, x & y
slitXCEN = xcen+0.5
slitYCEN = ycen+0.5
;angle of slit degrees cc (e) of north
slitPA = 95
;slitwidth in arcsec
slitWID= 2
;slit length in arcsec
;slitLEN=
;For now make slit length of image

;assume slit is built with 1" by wid " elements from center to edge?





END
