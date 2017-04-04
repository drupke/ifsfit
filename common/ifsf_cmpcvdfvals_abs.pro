; docformat = 'rst'
;
;+
;
; Compute velocities for emission lines based on cumulative velocity 
; distribution, as described in Zakamska & Greene 2013 (or 2014?).
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Hash absVEL[...
;
; :Params:
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
;      2016sep28, DSNR, created; copied from IFSF_CMPLINSPECPARS
;    
; :Copyright:
;    Copyright (C) 2016 David S. N. Rupke
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
function ifsf_cmpcvdfvals_abs,abscvdf

   bad = 1d99
   c = 299792.458

   mapsize = size(abscvdf['tau'])
   nvel = n_elements(abscvdf['vel'])

   absvel = hash()
   absvel['vpk'] = dblarr(mapsize[1],mapsize[2])+bad
   absvel['vsig'] = dblarr(mapsize[1],mapsize[2])+bad
   absvel['v%02'] = dblarr(mapsize[1],mapsize[2])+bad
   absvel['v%16'] = dblarr(mapsize[1],mapsize[2])+bad
   absvel['v%50'] = dblarr(mapsize[1],mapsize[2])+bad
   absvel['v%84'] = dblarr(mapsize[1],mapsize[2])+bad
   absvel['v%98'] = dblarr(mapsize[1],mapsize[2])+bad
   velkeys = absvel.keys()
   
   for i=0,mapsize[1]-1 do begin
      for j=0,mapsize[2]-1 do begin
         maxtau = max(abscvdf['tau',i,j,*],k)
         if k ne 0 then absvel['vpk',i,j] = abscvdf['vel',k]
         iv98 = value_locate(abscvdf['cumtaunorm',i,j,*],0.02)
         iv84 = value_locate(abscvdf['cumtaunorm',i,j,*],0.16)
         iv50 = value_locate(abscvdf['cumtaunorm',i,j,*],0.5)
         iv16 = value_locate(abscvdf['cumtaunorm',i,j,*],0.84)
         iv02 = value_locate(abscvdf['cumtaunorm',i,j,*],0.98)
         if iv98 ne -1 AND iv98 ne nvel-1 then $
            absvel['v%98',i,j] = abscvdf['vel',iv98]
         if iv84 ne -1 AND iv84 ne nvel-1 then $
            absvel['v%84',i,j] = abscvdf['vel',iv84]
         if iv50 ne -1 AND iv50 ne nvel-1 then $
            absvel['v%50',i,j] = abscvdf['vel',iv50]
         if iv16 ne -1 AND iv16 ne nvel-1 then $
            absvel['v%16',i,j] = abscvdf['vel',iv16]
         if iv02 ne -1 AND iv02 ne nvel-1 then $
            absvel['v%02',i,j] = abscvdf['vel',iv02]
         if iv50 ne -1 AND iv84 ne -1 AND iv16 ne -1 AND $
            iv50 ne nvel-1 AND iv84 ne nvel-1 AND iv16 ne nvel-1 then $
            absvel['vsig',i,j] = $
               ((abscvdf['vel',iv50] - abscvdf['vel',iv84])+$
               (abscvdf['vel',iv16] - abscvdf['vel',iv50]))/2d
      endfor
   endfor

   return,absvel

end
