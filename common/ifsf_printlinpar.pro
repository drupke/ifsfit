; docformat = 'rst'
;
;+
;
; Write emission line parameters to a text file.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    None.
;
; :Params:
;    outlines: in, required, type=strarr(nlines)
;      Names of lines for which to print line fluxes.
;    lun: in/out, required, type=integer
;      Logical unit for output file. If outfile is set, the lun is output to 
;      this file. If outfile is not set, the lun to which to write the data is
;      contained in this parameter.
;    col: in, optional, type=integer
;      IFS column of spectrum. Only used if outfile is not set.
;    row: in, optional, type=integer
;      IFS row of spectrum. Only used if outfile is not set.
;    maxncomp: in, optional, type=integer
;      Maximum number of components among all emission lines being fit. Only 
;      used if outfile is not set.
;    linepars: in, optional, type=structure
;      Output of IFSF_SEPFITPARS, containing line parameters. Only used if 
;      outfile is not set.
;      
; :Keywords:
;    outfile: in, optional, type=string
;      Full path and name of output file. If set, initializes file by opening
;      it, creating a lun, writing column titles, and returning the lun. If not
;      set, input lun is opened and written to.
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
;      2014jan13, DSNR, re-written to use hashes, to open file only once, and to
;                       print central wavelength and sigma as well as flux
;      2014jan16, DSNR, bugfixes
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
pro ifsf_printlinpar,outlines,lun,col,row,maxncomp,linepars,outfile=outfile

;  Initialize file
   if keyword_set(outfile) then begin

      openw,lun,outfile,/get_lun
      linestr = ''
      collabstr = ''
      foreach line,outlines do begin
        linestr += string(line,format='(A-48)')
        collabstr += string('Flux','Flux Error','Wave(A)','Sigma(km/s)',$
                            format='(4A12)')
      endforeach
      printf,lun,'#',linestr,format='(A-12,A0)'
      printf,lun,'#Col','Row','Cmp',collabstr,format='(A-4,2A4,A0)'

;  Print fluxes
   endif else begin
    
;     Cycle through components
      for i=0,maxncomp-1 do begin
         linestr=''
         foreach line,outlines do begin
            linestr += string(linepars.flux[line,i],$
                              linepars.fluxerr[line,i],$
                              linepars.wave[line,i],$
                              linepars.sigma[line,i],$
                              format='(E12.4,E12.4,D12.2,D12.2)')
         endforeach
         printf,lun,col,row,i+1,linestr,format='(3I4,A0)'
      endfor

   endelse

end
