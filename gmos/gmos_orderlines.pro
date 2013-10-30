pro gmos_orderlines,struct,zarr

  ppoff = struct.param[0]
  ncomp = struct.param[1]
  ppoff0 = ppoff - 2*(ncomp-1)

  linepars = sepfitpars(struct.param,struct.perror)
  nlines = n_elements(linepars.flux[*,0])

  haind = where(struct.linelabel eq 'Halpha')
  hafluxes = dblarr(ncomp)
  zfitarr = dblarr(ncomp)
  for i=0,ncomp-1 do begin
     hafluxes[i] = linepars.fluxpk[haind,i]
     if i eq 0 then zfitarr[0] = struct.param[2] $
     else zfitarr[i] = struct.param[ppoff0+2*i-1]
  endfor

  sortha = reverse(sort(hafluxes))

  zarr = zarr[sortha]
  zfitarr = zfitarr[sortha]
  for i=0,ncomp-1 do begin
     foff = ppoff + i*nlines*3 + indgen(nlines)*3
     woff = foff + 1
     soff = foff + 2
     struct.param[foff] = linepars.fluxpk[*,sortha[i]]
     struct.perror[foff] = linepars.fluxpkerr[*,sortha[i]]
     struct.param[woff] = linepars.wave[*,sortha[i]]
     struct.param[soff] = linepars.sigma[*,sortha[i]]
     if i eq 0 then struct.param[2] = zfitarr[0] $
     else struct.param[ppoff0+2*i-1] = zfitarr[i]
  endfor

end
