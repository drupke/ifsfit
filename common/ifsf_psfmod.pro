pro ifsf_psfmod

   dx = 79
   dy = 95
   ; center in single-offset coordinates
   center_nuclei = [43d,44d]
   platescale = 0.2914d
   modfac = 11d
   nmodx = dx*modfac
   nmody = dy*modfac
   nmodlam = 101d
   xmodran = [3450d,5600d]
   ; model wavelength range
   lammod = dindgen(nmodlam)/nmodlam*(xmodran[1]-xmodran[0])+xmodran[0]
   ; model x and y values, single-offset coordinates, in units of original pixels
   xmod = (dindgen(nmodx)-modfac/2d + 0.5d)/modfac + 1d
   ymod = (dindgen(nmody)-modfac/2d + 0.5d)/modfac + 1d
   ;map_xmod = rebin(xmod,nmodx,nmody)
   ;map_ymod = rebin(transpose(ymod),nmodx,nmody)
   ; model PSF fluxes, and rebinned values
   psfmod = dblarr(nmodx,nmody,nmodlam)
   psfrb = dblarr(dx,dy,nmodlam)
   pkflxratmod = dblarr(nmodlam)
   pkflxrat = dblarr(nmodlam)
   ; Peak value for Gaussian
   refpeak = 1d
   ; Sigma in model pixels; first factor is seeing in " at reflam
   refsig_pix = 1d/platescale*modfac/(2d*sqrt(2d*alog(2d)))
   ; Reference wavelength where sigma = seeing
   reflam = 5250d
   ; get indices of center
   xval = dindgen(dx)+1d
   yval = dindgen(dy)+1d
   ixcent = where(abs(xval - center_nuclei[0]) lt 1d/2d)
   iycent = where(abs(yval - center_nuclei[1]) lt 1d/2d)
   ; indices to xmod and ymod of center; use this logic to get nearest pixel
   ixcentmod = where(abs(xmod - center_nuclei[0]) lt 1d/modfac/2d)
   iycentmod = where(abs(ymod - center_nuclei[1]) lt 1d/modfac/2d)
   ; cycle through wavelengths
   for i=0, nmodlam-1 do begin
      ; 1d sigma in model pixels
      sig_pix = refsig_pix*(lammod[i]/reflam)^(-0.2d)
      ; express centroid in single-offset coordinates, so add 1 pixel
      psfmod[*, *, i] = $
         psf_gaussian([[refpeak, refpeak], [ixcentmod, iycentmod], $
         [sig_pix, sig_pix]], NPIXEL=[nmodx, nmody], NDIMENSION=2, $
         /DOUBLE)
      ; total flux / peak flux for model
      pkflxratmod[i] = total(psfmod[*,*,i]) / refpeak
      ; rebin model fluxes back to dx, dy grid; assume rebin is averaging and 
      ; conserve total flux
      psfrb[*,*,i] = rebin(psfmod[*,*,i],dx,dy)*modfac^2d
      ; total flux / flux in peak spaxel
      pkflxrat[i] = total(psfmod[*,*,i])/psfrb[ixcent,iycent,i]
      ; check that our indices are hitting the right spots; these differences
      ; should be 0
      ;print,psfmod[ixcentmod,iycentmod,i]-max(psfmod[*,*,i]),$
      ;   psfrb[ixcent,iycent,i]-max(psfrb[*,*,i])
   endfor
   
   ireflam = value_locate(lammod,reflam)
   pkflxratmod /= pkflxratmod[ireflam]
   pkflxrat /= pkflxrat[ireflam]

   set_plot,'x'
   pkflxrat_anl = (lammod/reflam)^(-0.4)
   ;pkflxrat_anl /= pkflxrat_anl[ireflam]
   cgplot,lammod,pkflxrat_anl,xran=xmodran,yran=[0.9,1.06],$
      color='red',psym='open circle',/xsty,/ysty
   cgoplot,lammod,pkflxratmod,thick=2
   cgoplot,lammod,pkflxrat,color='blue',thick=4
   ;al_legend,['
   ;return
end