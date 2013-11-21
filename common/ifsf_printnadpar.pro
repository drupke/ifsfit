; docformat = 'rst'
;
;+
;
; Write Na D absorption line properties to a text file. This cycles
; through all spectra and writes the file all at once.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    None.
;
; :Params:
;    gal: in, required, type=string
;      Galaxy label in file naming.
;    cols: in, required, type=dblarr(Nspec)
;      List of column #s for each spectrum.
;    row: in, required, type=dblarr(Nspec)
;      List of row #s for each spectrum.
;    fitdir: in, required, type=string
;      Directory where normalized spectra and fit parameter files are
;      located.
;    outfile: in, required, type=string
;      Name of output text file. (Path specified by fitdir.)
;    z: in, required, type=double
;      Galaxy redshift.
;
; :Keywords:
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
pro ifsf_printnadpar,gal,cols,rows,fitdir,outfile,z

  nad1_rest = 5895.92d
  nad2_rest = 5889.95d

  nspec = n_elements(rows)

  openw,lun,fitdir+outfile+'.dat',/get_lun
  printf,lun,'#Col','Row','Cmp','z','FWHM','tau','C_f','W_eq',$
         format='(A-4,2A4,5A10)'

  for i=0,nspec-1 do begin

;    Syntax for file naming:
;      Emission-line subtracted spectrum: galaxy_col_row_nad_spec.dat
;      Best fit parameters: galaxy_col_row_nad_parout.dat
     spec = string(cols[i],'_',rows[i],format='(I04,A0,I04)')
     spectot= fitdir+gal+'_'+spec
     specin = spectot+'_nad_spec.dat'
     parout = spectot+'_nad_parout.dat'
;    Get best fit parameters
     ifsf_readnadpar,parout,abspars,empars,opars
;    Compute equivalent width
     dumy = ifsf_cmpnadfull(specin,parout,z,weq=weq)
     printf,lun,cols[i],rows[i],1,abspars[2,0]/nad1_rest-1d,$
            abspars[3,0]*2d*sqrt(alog(2d)),abspars[1,0],abspars[0,0],$
            weq,format='(3I4,D10.6,D10.2,D10.4,D10.4,D10.4)'
     
  endfor
  
  free_lun,lun

end
