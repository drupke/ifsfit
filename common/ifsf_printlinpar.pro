; docformat = 'rst'
;
;+
;
; Write emission line fluxes and flux errors to a text
; file. Optionally include Na D equivalent width.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    None.
;
; :Params:
;    lines: in, required, type=dblarr(Nlines)
;      Names of lines.
;    fluxes: in, required, type=dblarr(Nlines\,Ncomp)
;      Fluxes of lines.
;    fluxerrors: in, required, type=dblarr(Nlines\,Ncomp)
;      Flux errors.
;    outfile: in, required, type=string
;      Full path and name of output file.
;    col: in, required, type=string
;      IFS column of spectrum.
;    row: in, required, type=string
;      IFS row of spectrum.
;
; :Keywords:
;    append: in, optional, type=byte, default=0
;      Append to previously existing file?
;    init: in, optional, type=byte, default=0
;      Initialize a new file? This writes column headers but no data.
;    nadweq: in, optional, type=dblarr(2), default=0
;      Write Na D equivalent width and error as last two columns? If
;      so, set this to these values. Assumes one absorption component
;      and attaches to first emission-line component.
; 
; :Author:
;    David S. N. Rupke::
;      Rhodes College
;      Department of Physics
;      2000 N. Parkway
;      Memphis, TN 38104
;      drupke@gmail.com
;
; :History:
;    ChangeHistory::
;      2009may, DSNR, created
;      2009jun07, DSNR, rewritten
;      2013nov21, DSNR, documented, renamed, added license and copyright 
;    
; :Copyright:
;    Copyright (C) 2013 David S. N. Rupke
;
;    This program is free software: you can redistribute it and/or
;    modify it under the terms of the GNU General Public License as
;    published by the Free Software Foundation, either version 3 of
;    the License or any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;    General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program.  If not, see
;    http://www.gnu.org/licenses/.
;
;-
pro ifsf_printlinpar,lines,fluxes,fluxerrors,outfile,col,row,$
                     append=append,init=init,nadweq=nadweq

  if ~ keyword_set(append) then append=0

  nlines = n_elements(lines)

  fluxes_info = size(fluxes)
  if (fluxes_info[0] gt 1) then ncomp = fluxes_info[2] else ncomp=1

; Open file for output
  if (append) then openu,unit,outfile,/get_lun,append=append $
  else openw,unit,outfile,/get_lun

  for j=0,ncomp-1 do begin

     out = ''
     for i=0,nlines-1 do begin
        if keyword_set(init) then $
           out += string(lines[i],'1sigma',$
                         format='(A12,A12)') $
        else begin
           ind = where(lines eq lines[i])
           out += string(fluxes[ind,j],fluxerrors[ind,j],$
                         format='(E12.4,E12.4)')
        endelse
     endfor
     if keyword_set(nadweq) then begin
        if j eq 0 then begin
           if keyword_set(init) then $
              out += string('NaI D','1sigma',$
                            format='(A12,A12)') $
           else begin
              out += string(nadweq[0],nadweq[1],$
                            format='(E12.4,E12.4)')
           endelse
        endif else $
           out += string(-1,-1,format='(2I12)')
     endif

     if keyword_set(init) then $
        printf,unit,'#Col','Row','Cmp',out,format='(A-4,2A4,A0)' $
     else $
        printf,unit,col,row,j+1,out,format='(3I4,A0)'

  endfor

  free_lun,unit

end
