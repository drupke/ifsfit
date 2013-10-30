function getbalmfluxes,labels,fluxes,errors=errors
;
; History
;   09may26  DSNR  created
;   09jun07  DSNR  modified to include errors
;

;blines = ['Halpha','Hbeta','Hgamma','Hdelta','Hepsilon',$
;          'H8','H9','H10','H11','H12']
;nlines = 10
  blines = ['Halpha','Hbeta','Hgamma','Hdelta','Hepsilon']
  nlines = 5

  fluxes_info = size(fluxes)
  if (fluxes_info[0] gt 1) then ncomp = fluxes_info[2] else ncomp=1

  if ~ keyword_set(errors) then bfluxes = {value: dblarr(nlines,ncomp) } $
  else bfluxes = {value: dblarr(nlines,ncomp), error: dblarr(nlines,ncomp)}

  for i=0,nlines-1 do begin
     ind = where(labels eq blines[i],ct)
     if ct gt 0 then begin
        for j=0,ncomp-1 do begin
           bfluxes.value[i,j] = fluxes[ind,j]
           if keyword_set(errors) then bfluxes.error[i,j] = errors[ind,j]
        endfor
     endif else begin
        bfluxes.value[i,*] = 0d
        if keyword_set(errors) then bfluxes.error[i,*] = 0d
     endelse
  endfor

  return,bfluxes

end
