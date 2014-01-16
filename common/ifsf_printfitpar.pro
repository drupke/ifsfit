; docformat = 'rst'
;
;+
;
; Write fit parameters to a text file.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    None.
;
; :Params:
;    lun: in/out, required, type=integer
;      Logical unit for output file. If outfile is set, the lun is output to
;      this file. If outfile is not set, the lun to which to write the data is
;      contained in this parameter.
;    col: in, optional, type=integer
;      IFS column of spectrum. Only used if outfile is not set.
;    row: in, optional, type=integer
;      IFS row of spectrum. Only used if outfile is not set.
;    struct: in, optional, type=structure
;      Structure output from IFSF. Only used if outfile is not set.
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
;      2013jan13, DSNR, copied from IFSF_PRINTLINPAR and modified
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
pro ifsf_printfitpar,lun,col,row,struct,outfile=outfile

;  Initialize file
   if keyword_set(outfile) then begin

      openw,lun,outfile,/get_lun
      printf,lun,'#Col','Row','Cmp','Rchi2','Niter',$
             format='(A-4,2A4,2A6)'

;  Print data
   endif else begin
    
     printf,lun,col,row,struct.redchisq,struct.niter,$
            format='(I4,I4,I4,D6.2,I6)'


   endelse

end
