pro lrisbr_printlinepars,lines,$
                         fluxesB,fluxerrorsB,$
                         fluxesR,fluxerrorsR,$
                         outfile,speclabel,append=append,$
                         whichlines=whichlines,init=init,$
                         nobcor=nobcor
;+
; History
;  10mar19  DSNR  created
;-

  if ~ keyword_set(append) then append=0
  if ~ keyword_set(whichlines) then whichlines = lines
  if keyword_set(nobcor) then br=0 else br=1

  nlines = n_elements(whichlines)

  fluxesBinfo = size(fluxesB)
  fluxesRinfo = size(fluxesR)
  if (fluxesBinfo[0] gt 1) then ncompB = fluxesBinfo[2] else ncompB=1
  if (fluxesBinfo[1] eq 1) then ncompB = 0
  if (fluxesRinfo[0] gt 1) then ncompR = fluxesRinfo[2] else ncompR=1
  if (fluxesRinfo[1] eq 1) then ncompR = 0

; Open file for output
  if (append) then openu,unit,outfile,/get_lun,append=append $
  else openw,unit,outfile,/get_lun

  ncomp = max([ncompB,ncompR])

  for j=0,ncomp-1 do begin

     out = ''
     for i=0,nlines-1 do begin
        if keyword_set(init) then $
           out += string(whichlines[i],'1sigma',format='(A12,A12)') $
        else begin
           zero=1
           ind = where(lines eq whichlines[i])
           if ncompB ge j+1 then begin
              if fluxesB[ind,j] gt 0 then begin
                 out += string(fluxesB[ind,j],fluxerrorsB[ind,j],$
                               format='(E12.4,E12.4)')
                 zero=0
              endif
           endif
           if ncompR ge j+1 then begin
              if fluxesR[ind,j] gt 0 then begin
                 out += string(fluxesR[ind,j],fluxerrorsR[ind,j],$
                               format='(E12.4,E12.4)')
                 zero=0
              endif
           endif
           if zero then $
              out += string(0,0,format='(E12.4,E12.4)')
        endelse
     endfor
     
     if keyword_set(init) then $
        printf,unit,'# ID','comp','br?',out,format='(A-15,A5,A4,A0)' $
     else $
        printf,unit,speclabel,j+1,br,out,format='(A-15,I5,I3,A0)'

  endfor

  if ncomp eq 0 then begin

     out = ''
     for i=0,nlines-1 do begin
        out += string(0,0,format='(E12.4,E12.4)')
     endfor
     printf,unit,speclabel,1,br,out,format='(A-15,I5,I3,A0)'
     
  endif

  free_lun,unit
     
end
