; docformat = 'rst'
;
;+
;
; Grid Valdes et al. 2004 (ApJS, 152, 221) stellar spectra 
; to select a representative subset. See Figure 1 in Valdes+04.
;
; Default is to pick the "zeroth" element of each bin. Randomize this by
; selecting the RANDOM keyword.
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
;    random: in, optional, type=boolean
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
pro ifsf_gridvaldestemp,infile,outfile,random=random

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
            if not keyword_set(random) then iuse = 0 $
            else begin
               ran = randomu(seed)
               iuse = fix((ct-1)*ran)
            endelse
            if firsttemp then begin
               newflux = template.flux[*,ind[iuse]]
               newobj = template.obj[ind[iuse]]
               newteff = template.teff[ind[iuse]]
               newlogg = template.logg[ind[iuse]]
               newfeh = template.feh[ind[iuse]]
               firsttemp = 0b
            endif else begin
               newflux = [[newflux],[template.flux[*,ind[iuse]]]]
               newobj = [newobj,template.obj[ind[iuse]]]
               newteff = [newteff,template.teff[ind[iuse]]]
               newlogg = [newlogg,template.logg[ind[iuse]]]
               newfeh = [newfeh,template.feh[ind[iuse]]]
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
               if not keyword_set(random) then iuse = 0 $
               else begin
                  ran = randomu(seed)
                  iuse = fix((ct-1)*ran)
               endelse
               newflux = [[newflux],[template.flux[*,ind[iuse]]]]
               newobj = [newobj,template.obj[ind[iuse]]]
               newteff = [newteff,template.teff[ind[iuse]]]
               newlogg = [newlogg,template.logg[ind[iuse]]]
               newfeh = [newfeh,template.feh[ind[iuse]]]
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