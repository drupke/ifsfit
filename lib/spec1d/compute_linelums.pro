function compute_linelums, linefit, distance, distance_err, $
  select_lines=select_lines, photerror=photerror
; jm03may11uofa
; compute emission-line luminosities

    nspec = n_elements(linefit)
    
    if n_elements(photerror) eq 0L then photerror = 0.0 ; [%]
    
    lsun = 3.845D33             ; erg/s
    Mpc2cm = 3.086D24           ; conversion from Mpc to cm

    nselect = n_elements(select_lines)
    linelums = create_struct(select_lines[0]+'_lum',-999.0,select_lines[0]+'_lum_err',-999.0)
    for k = 1L, nselect-1L do linelums = create_struct(linelums,select_lines[k]+'_lum',-999.0,$
      select_lines[k]+'_lum_err',-999.0)
    linelums = replicate(linelums,nspec)
    
;   if distance lt -900.0 then return, linelums else dist = distance*Mpc2cm
;   if distance_err lt -900.0 then dist_err = 0.0 else dist_err = distance_err*Mpc2cm

    for j = 0L, nselect-1L do begin

       mini = struct_trimtags(linefit,select=select_lines[j])

       array = mini.(0)
       flux = reform(array[0,*])
       flux_err = reform(array[1,*])
;      flux_err = sqrt( (mini.(0))[1]^2.0 + (photerror*flux/100.0)^2.0 )

       good = where((distance gt -900.0) and (flux_err gt 0.0),ngood)
       if ngood ne 0L then begin

          dist = distance[good]*Mpc2cm
          nodisterr = where(distance_err[good] lt -900.0,nnodisterr)
          if nnodisterr ne 0L then distance_err[good] = 0.0
          dist_err = distance_err[good]*Mpc2cm

          flux = flux[good]
          flux_err = sqrt( (flux_err[good])^2.0 + (photerror*flux/100.0)^2.0 )
          
          lum = flux*4.0*!dpi*dist*dist                                                ; [erg/s]
          lum_err = 4.0*!dpi*sqrt((2*flux*dist*dist_err)^2.0+(dist*dist*flux_err)^2.0) ; [erg/s]

; take the logarithm and fill the structure       
       
          linelums[good].(2*j+1) = lum_err / lum / alog(10.0) ; luminosity error
          linelums[good].(2*j) = alog10(lum/lsun)             ; luminosity

       endif
          
    endfor

return, linelums
end    
