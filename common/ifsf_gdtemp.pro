; docformat = 'rst'
;
;+
;
; Load Gonzalez-Delgado et al. 2005 population synthesis models into a
; structure of wavelength and flux arrays. Population synthesis models
; from Gonzalez-Delgado et al. 2005, MNRAS, 357, 945.  See
; http://www.iaa.csic.es/~rosa/research/synthesis/HRES/ESPS-HRES.html
; for more details.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    IDL save file. The file consists of a structure (named
;    'template') containing three tags: lambda, flux, and ages. Lambda
;    is an n-element array, while flux is and nxm array of the format
;    [n_wavelengths, n_ages]. Ages is the array of age values.
;
; :Params:
;    infile: in, required, type=string
;      Filename and path of input model table.
;    outfile: in, required, type=string
;      Filename and path of output file.
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
;      2009aug06, DSNR, created
;      2014feb17, DSNR, complete rewrite for ease of use/customization;
;                       added detailed documentation
;    
; :Copyright:
;    Copyright (C) 2014 David S. N. Rupke
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
pro ifsf_gdtemp,infile,outfile

  openr,lun,infile,/get_lun
  header = strarr(43)
  readf,lun,header
  ages = strarr(1)
  readf,lun,ages
  header = strarr(1)
  readf,lun,header
  data = dblarr(139,13321)
  readf,lun,data
  free_lun,lun

  ages = strsplit(ages,/extract)
  ages = double(ages[3:48])

  ilum = dindgen(46)*3+1
  lambda = data[0,*]
  flux = data[ilum,*]
  flux = transpose(flux)
  
  for i=0,n_elements(flux[0,*])-1 do $
     flux[*,i] /= max(flux[*,i])
  
  template = {lambda: lambda, flux: flux, ages: ages}

  save,template,filename=outfile
  
end
