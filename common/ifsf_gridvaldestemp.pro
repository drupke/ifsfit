; docformat = 'rst'
;
;+
;
; Grid Valdes et al. 2004 (ApJS, 152, 221) stellar spectra 
; to select a representative subset. See Figure 1 in Valdes+04.
;
; Input file can first be created with IFSF_VALDESTEMP.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    IDL save file. The file consists of a structure (named
;    'template') containing several tags.
;
; :Params:
;    infile: in, required, type=string
;      Filename and path of input templates.
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
;      2018apr27, DSNR, created
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
pro ifsf_gridvaldestemp,infile,outfile

   restore,infile

;  Loop through g, z plane
   dlogg = 0.25d
   dfeh = 0.2d
   logg = 0d
   firsttemp = 1b
   while logg lt 5d do begin
      feh = -1d
      while feh lt 0.4d do begin
         ind = where(template.logg ge logg AND template.logg lt logg+dlogg AND $
                     template.feh ge feh AND template.feh lt feh+dfeh,ct)    
         if ct gt 0 then begin
            if firsttemp then begin
               newflux = template.flux[*,ind[0]]
               newobj = template.obj[ind[0]]
               newteff = template.teff[ind[0]]
               newlogg = template.logg[ind[0]]
               newfeh = template.feh[ind[0]]
               firsttemp = 0b
            endif else begin
               newflux = [[newflux],[template.flux[*,ind[0]]]]
               newobj = [newobj,template.obj[ind[0]]]
               newteff = [newteff,template.teff[ind[0]]]
               newlogg = [newlogg,template.logg[ind[0]]]
               newfeh = [newfeh,template.feh[ind[0]]]
            endelse
         endif
         feh+=dfeh
      endwhile
      logg+=dlogg
   endwhile

;  Loop through teff, z plane, upper left quadrant
   dlteff = 0.1d
   dfeh = 0.2d
   lteff = 3.9d
   logteff = alog10(template.teff)
   while lteff lt 4.5d do begin
      feh = -1d
      while feh lt 0.4d do begin
         ind = where(logteff ge lteff AND logteff lt lteff+dlteff AND $
                     template.feh ge feh AND template.feh lt feh+dfeh,ct)
         if ct gt 0 then begin
            isin = 0b
            for i=0,ct-1 do begin
               indtmp = where(newobj eq template.obj[ind[i]],cttmp)
               if cttmp gt 0 then isin = 1b
            endfor
            if isin eq 0b then begin
               newflux = [[newflux],[template.flux[*,ind[0]]]]
               newobj = [newobj,template.obj[ind[0]]]
               newteff = [newteff,template.teff[ind[0]]]
               newlogg = [newlogg,template.logg[ind[0]]]
               newfeh = [newfeh,template.feh[ind[0]]]
            endif
         endif
         feh+=dfeh
      endwhile
      lteff+=dlteff
   endwhile


   template = {lambda: template.lambda, flux: newflux, obj: newobj, $
               teff: newteff, logg: newlogg, feh: newfeh}
   save,template,filename=outfile

end