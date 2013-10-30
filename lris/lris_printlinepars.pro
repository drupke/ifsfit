pro lris_printlinepars,lines,fluxes,fluxerrors,outfile,speclabel,$
                       append=append,whichlines=whichlines,$
                       init=init,nobcor=nobcor,zeroflux=zeroflux
;
; History
;  09mayXX  DSNR  created
;  09jun07  DSNR  rewritten

  if ~ keyword_set(append) then append=0
  if ~ keyword_set(whichlines) then whichlines = lines
  if keyword_set(nobcor) then br=0 else br=1

  nlines = n_elements(whichlines)

  if keyword_set(zeroflux) then begin
     ncomp = 1
     br = 0
  endif else begin
     fluxes_info = size(fluxes)
     if (fluxes_info[0] gt 1) then ncomp = fluxes_info[2] else ncomp=1
  endelse

; Open file for output
  if append then openu,unit,outfile,/get_lun,append=append $
  else openw,unit,outfile,/get_lun
  
  for j=0,ncomp-1 do begin

     if keyword_set(zeroflux) then begin
        out = ''
        for i=0,nlines-1 do $
           out += string(0d,0d,format='(E12.4,E12.4)')
     endif else begin

        out = ''
        for i=0,nlines-1 do begin
           if keyword_set(init) then $
              out += string(whichlines[i],'1sigma',format='(A12,A12)') $
           else begin
              ind = where(lines eq whichlines[i])
              out += string(fluxes[ind,j],fluxerrors[ind,j],$
                            format='(E12.4,E12.4)')
           endelse
        endfor

     endelse
        
     if keyword_set(init) then $
        printf,unit,'# ID','comp','br?',out,format='(A-15,A5,A4,A0)' $
     else $
        printf,unit,speclabel,j+1,br,out,format='(A-15,I5,I3,A0)'

  endfor

  free_lun,unit

end
