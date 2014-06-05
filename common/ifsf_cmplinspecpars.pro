; docformat = 'rst'
;
;+
;
; Compute velocities for emission lines based on cumulative velocity 
; distribution, as described in Zakamska & Green 2013 (or 2014?).
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Structure of velocities and fluxes at those velocities, with each tag
;    being an (Ncols x Nrows) array. See end of file for list of tags.
;
; :Params:
;    linmap: in, required, type=structure
;      Output line spectra from IFSF_CMPLINSPECMAPS.
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
;      2014jun03, DSNR, created
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
function ifsf_cmplinspecpars,linmap

   bad = 1d99
   c = 299792.458

   mapsize = size(linmap.flux)
   nvel = n_elements(linmap.vel)

   vpk = dblarr(mapsize[1],mapsize[2])+bad
   v98 = dblarr(mapsize[1],mapsize[2])+bad
   v84 = dblarr(mapsize[1],mapsize[2])+bad
   v50 = dblarr(mapsize[1],mapsize[2])+bad
   sig = dblarr(mapsize[1],mapsize[2])+bad
   fpk = dblarr(mapsize[1],mapsize[2])
   fv98 = dblarr(mapsize[1],mapsize[2])
   fv84 = dblarr(mapsize[1],mapsize[2])
   fv50 = dblarr(mapsize[1],mapsize[2])
   for i=0,mapsize[1]-1 do begin
      for j=0,mapsize[2]-1 do begin
         fpk[i,j] = max(linmap.flux[i,j,*],k)
         if k ne 0 then vpk[i,j] = linmap.vel[k]
         iv98 = value_locate(linmap.cumfluxnorm[i,j,*],0.02)
         iv84 = value_locate(linmap.cumfluxnorm[i,j,*],0.1)
         iv50 = value_locate(linmap.cumfluxnorm[i,j,*],0.5)
         if iv98 ne nvel-1 then begin
            v98[i,j] = linmap.vel[iv98]
            fv98[i,j] = linmap.flux[i,j,iv98]
         endif
         if iv84 ne nvel-1 then begin
            v84[i,j] = linmap.vel[iv84]
            fv84[i,j] = linmap.flux[i,j,iv84]
         endif
         if iv50 ne nvel-1 then begin
            v50[i,j] = linmap.vel[iv50]
            fv50[i,j] = linmap.flux[i,j,iv50]
            if iv84 ne nvel-1 then $
               sig[i,j] = linmap.vel[iv50] - linmap.vel[iv84]
         endif
      endfor
   endfor
   ftot = total(linmap.flux,3)   

   return,{vpk:vpk,v98:v98,v84:v84,v50:v50,sig:sig,$
           fpk:fpk,fv98:fv98,fv84:fv84,fv50:fv50,ftot:ftot}

end
