; docformat = 'rst'
;
;+
;
; Load Valdes et al. 2004 (ApJS, 152, 221) stellar spectra into a
; structure of wavelength and flux arrays.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    IDL save file. The file consists of a structure (named
;    'template') containing three tags: lambda, flux, and temps. Lambda
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
;      2018apr14, DSNR, copied from IFSF_GDTEMP
;    
; :Copyright:
;    Copyright (C) 2018 David S. N. Rupke
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
pro ifsf_s7smoothvaldestemp,infile,outfile,z,stichwave=stitchwave

   bad = 1d99
   if ~ keyword_set(stitchwave) then stitchwave=5600d
   disp = 0.4d

;  get template
   restore,infile
   sizetemp = size(template.flux)
;  redshift template
   zlambda = template.labmda * (1d + z)
;  blue parameters
   iblue = value_locate(zlambda,stitchwave)
   sigblue = ((1.35d/2.35d)^2d)
   bluevscale = disp/mean([template.lambda[0],template.lambda[iblue]])*299792d
   redvscale = disp/mean([template.lambda[iblue+1],template.lambda[sizetemp[1]-1]])*299792d
   for i=0,sizetemp[2]-1 do begin
;     rebin blue side
      log_rebin,[zlambda[0],zlambda[iblue]],$
                 template.flux[0:iblue,i],bluelog,blueloglam,$
                 velscale=bluevscale
;     rebin red side
      log_rebin,[zlambda[iblue+1],zlambda[sizetemp[1]-1]],$
                 template.flux[iblue+1:sizetemp[1]-1,i],redlog,redloglam,$
                 velscale=redvscale
      bluelog_sm = $
         gauss_smooth(bluelog,
                      
                      bluesmooth
                      interpol(continuum_log,gdlambda_log,ALOG(gdlambda))
   endfor

   template = {lambda: wave, flux: flux, obj: obj, $
               teff: teff, logg: logg, feh: feh}

end
