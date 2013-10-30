pro nifs_printlinepars,lines,fluxes,fluxerrors,outfile,col,row,$
                       append=append,whichlines=whichlines,init=init,$
                       doubleline=doubleline
;
; History
;  13mar05  DSNR  created
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
           ind = where(lines eq whichlines[i],ctlines)
           if ctlines gt 0 then begin
;             If line is specified to be fit, but out of the fit
;             range, it can return a positive flux and a zero
;             errors. Ignore these lines.
              if fluxerrors[ind,j] eq 0 then $
                 out += string(0,0,format='(E12.4,E12.4)') $
              else $
                 out += string(fluxes[ind,j],fluxerrors[ind,j],$
                               format='(E12.4,E12.4)')
;             If a line is specified to be output but NOT to be fit,
;             output 0 flux and flux errors.
           endif else begin
              out += string(0,0,format='(E12.4,E12.4)')
           endelse
        endelse
     endfor

     if (keyword_set(init)) then $
        printf,unit,'#Col','Row','Cmp',out,format='(A-4,2A4,A0)' $
     else begin
        cmp=j+1
        if doubleline eq 'paa' AND j eq 2 then cmp=j+2
        printf,unit,col,row,cmp,out,format='(3I4,A0)'
     endelse

  endfor

  free_lun,unit

end
