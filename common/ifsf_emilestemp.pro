; docformat = 'rst'
;
;+
;
; Load E-MILES population synthesis models into a
; structure of fluxlength and flux arrays. Population synthesis models
; from Vazdekis et al. 2016.  See
; http://research.iac.es/proyecto/miles/pages/spectral-energy-distributions-seds/e-miles.php
; for more details.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    IDL save file. The file consists of a structure (named
;    'template') containing three tags: lambda, flux, and ages. Lambda
;    is an n-element array, while flux is and nxm array of the format
;    [n_fluxlengths, n_ages]. Ages is the array of age values.
;
; :Params:
;    indir: in, required, type=string
;      Path of input model FITS files.
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
;      2019nov13, DSNR, created
;    
; :Copyright:
;    Copyright (C) 2019 David S. N. Rupke
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
pro ifsf_emilestemp,indir,outfile,wavelo,wavehi,zran=zran,tran=tran

   bad = 1d99

;  get filenames
   fitsfiles = file_search(indir+'*.fits')
  
;  process each file
   firstfile = 1b
   foreach file,fitsfiles do begin
      nouse = 0b
      spec = ifsr_readspec(file,/nodcflag,/noctype)
      filepluspath = strsplit(file,'/',/extract)
      filenopath = filepluspath[n_elements(filepluspath)-1]
      zpos = strpos(filenopath,'Z')
      Zsign = filenopath.substring(zpos+1,zpos+1)
      Z = double(filenopath.substring(zpos+2,zpos+5))
      if Zsign eq 'm' then Z*=-1d
      tpos = strpos(filenopath,'T')
;     convert from Gyr to Myr
      T = double(filenopath.substring(tpos+1,tpos+7)) * 1d9
      thiswaveall = spec[*,0]
      ilo = value_locate(thiswaveall,wavelo)
      ihi = value_locate(thiswaveall,wavehi)
      if keyword_set(tran) then begin
         if T lt tran[0] OR T gt tran[1] then nouse=1b
      endif else if keyword_set(zran) then begin
         if Z lt zran[0] OR Z gt zran[1] then nouse=1b
      endif
      if not nouse then begin
         if firstfile then begin
            waveall = thiswaveall[ilo:ihi]
            fluxall = spec[ilo:ihi,1]
            zall = Z
            tall = T
            firstfile = 0b
         endif else begin
;            waveall = [[waveall],[thiswaveall[ilo:ihi]]]
            fluxall = [[fluxall],[spec[ilo:ihi,1]]]
            zall = [zall,Z]
            tall = [tall,T]
         endelse
      endif
   endforeach
   
   template = {lambda: waveall, flux: fluxall, ages: tall, logzh: zall}
   save,template,filename=outfile
      
end
