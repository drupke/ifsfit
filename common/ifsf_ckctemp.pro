; docformat = 'rst'
;
;+
;
; Load CKC models into a structure of wavelength and flux arrays.
;
; From Christy:
; "The files are unpublished “CKC” models (Conroy, Kurucz, and someone else who 
; starts with a C) that Charlie gave me in 2014. They are fully theoretical. 
; The native resolution of the models is quite high (R~10,000). These have 
; been convolved to match MaGEs instrumental resolution of 30 km/s and 
; re-sampled in log-lambda wavelength bins (since that’s what PPXF wants). 
; They are also normalized near 5500 A. There are 43 ages here, which is 
; about as many as PPXF can handle. However, I often get away with using a 
; much sparser age grid (as few as 9 ages.) If you’d prefer the sparse age 
; grid, different wavelength coverage, or spectral resolution let me know. 
; I tend to custom repackage the original models for each task. (I can send 
; you the originals too if you’d like, though they are more unwieldy.)"
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
;      2020apr14, DSNR, copied from IFSF_VALDESTMP
;    
; :Copyright:
;    Copyright (C) 2020 David S. N. Rupke
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
pro ifsf_ckctemp,infile,outfile,wavelo=wavelo,wavehi=wavehi

   bad = 1d99

   if ~ keyword_set(wavelo) then wavelo = 3465.0d
   if ~ keyword_set(wavehi) then wavehi = 7000.0d
   
   s = mrdfits(infile, 1)

   template = {lambda: s.wave, flux: s.flux, $
               ages: s.age, sigma: s.sigma, dv: s.dv}

   save,template,filename=outfile
  
end
