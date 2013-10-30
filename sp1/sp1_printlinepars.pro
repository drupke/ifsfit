pro sp1_printlinepars,lines,fluxes,fluxerrors,outfile,$
                      append=append,whichlines=whichlines,init=init

;
; History
;   10nov09  DSNR  copied from GMOS_PRINTLINEPARS
;

  if ~ keyword_set(append) then append=0
  if ~ keyword_set(whichlines) then whichlines = lines

  nlines = n_elements(whichlines)

  fluxes_info = size(fluxes)
  if (fluxes_info[0] gt 1) then ncomp = fluxes_info[2] else ncomp=1

; Open file for output
  if (append) then openu,unit,outfile,/get_lun,append=append $
  else openw,unit,outfile,/get_lun

  for j=0,ncomp-1 do begin

     out = ''
     for i=0,nlines-1 do begin
        if keyword_set(init) then $
           out += string(whichlines[i],'1sigma',$
                         format='(A12,A12)') $
        else begin
           ind = where(lines eq whichlines[i])
           out += string(fluxes[ind,j],fluxerrors[ind,j],$
                         format='(E12.4,E12.4)')
        endelse
     endfor
        

     if (keyword_set(init)) then $
        printf,unit,'#Cmp',out,format='(A-4,A0)' $
     else $
        printf,unit,j+1,out,format='(I4,A0)'

  endfor

  free_lun,unit

end
