; docformat = 'rst'
;
;+
;
; Write Na D absorption line properties to a text file.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    None.
;
; :Params:
;    lun: in, required, type=integer
;      Logical unit for output file. If outfile is set, the lun is output to
;      this file. If outfile is not set, the lun to which to write the data is
;      contained in this parameter as an output.
;    col: in, optional, type=integer
;      IFS column of spectrum. Only used if outfile is not set.
;    row: in, optional, type=integer
;      IFS row of spectrum. Only used if outfile is not set.
;    param: in, optional, type=structure
;      Output of NaD fit, containing line parameters. Only used if 
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
;      2010jul22, DSNR, created
;      2013nov21, DSNR, documented, renamed, added license and copyright
;      2014may13, DSNR, re-written to match logic of IFSF_PRINTLINPAR,
;                       and for new fitting method
;    
; :Copyright:
;    Copyright (C) 2013-2014 David S. N. Rupke
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
pro ifsf_printnadpar,lun,col,row,param,outfile=outfile

;  Initialize file
   if keyword_set(outfile) then begin

      openw,lun,outfile,/get_lun
      linestr = string('#','HeI 5876 emission','NaD absorption','NaD emission',$
                       format='(A-12,A-32,A-42,A-40)')
      collabstr = string('#Col','Row','Cmp',$
                         'Wave','Sigma','Flux',$
                         'C_f','tau','Wave','sigma',$
                         'Wave','Sigma','Flux','D2/D1',$
                         format='(A-4,2A4,2A12,2A8,A10,4A12,2A8)')
      printf,lun,linestr,collabstr

;  Print line parameters
   endif else begin

      nhei=param[0]
      nnadabs=param[1]
      nnadem=param[2]
      nmax = max(param[0:2])


      for i=0,nmax-1 do begin
         outstr=''
         ilo = 3
         if nhei ge i+1 then $
            outstr += string(param[ilo+i*3:ilo+(i+1)*3-1],$
                              format='(D12.2,D12.2,D8.4)') $
         else outstr += string(0,0,0,format='(D12.2,D12.2,D8.4)')
         ilo = 3+nhei*3
         if nnadabs ge i+1 then $
            outstr += string(param[ilo+i*4:ilo+(i+1)*4-1],$
                              format='(D8.4,D10.4,D12.2,D12.2)') $
         else outstr += string(0,0,0,0,format='(D8.4,D10.4,D12.2,D12.2)')
         ilo = 3+nhei*3+nnadabs*4
         if nnadem ge i+1 then $
            outstr += string(param[ilo+i*4:ilo+(i+1)*4-1],$
                              format='(D12.2,D12.2,D8.4,D8.4)') $
         else outstr += string(0,0,0,0,format='(D12.2,D12.2,D8.4,D8.4)')

         printf,lun,col,row,i+1,outstr,format='(3I4,A0)'
     
      endfor
  
   endelse
  
end
