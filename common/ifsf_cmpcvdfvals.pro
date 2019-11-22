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
;    Hash EMLVEL[...
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
function ifsf_cmpcvdfvals,emlcvdf,emlflx,emlflxerr,fluxvels=fluxvels,$
                          sigcut=sigcut

   bad = 1d99
   c = 299792.458

   outlines = emlflx['ftot']->keys()
   mapsize = size(emlflx['ftot',outlines[0]])
   nvel = n_elements(emlcvdf['vel'])

   emlvel = hash()
   emlvel['vpk'] = hash()
   emlvel['vsig'] = hash()
   emlvel['v%02'] = hash()
   emlvel['v%16'] = hash()
   emlvel['v%50'] = hash()
   emlvel['v%84'] = hash()
   emlvel['v%98'] = hash()
   velkeys = emlvel.keys()
   
   emlflx['ftotpk'] = hash()
   emlflxerr['ftotpk'] = hash()
   flxkeys = ['ftotpk']
   emlflx['fv%50'] = hash()
   emlflxerr['fv%50'] = hash()
   flxkeys = [flxkeys,'fv%50']
   if keyword_set(fluxvels) then begin
      foreach vel,fluxvels do begin
         fkey = 'fv'+string(vel,format='(I0)')
         emlflx[fkey] = hash()
         emlflxerr[fkey] = hash()
         flxkeys = [flxkeys,fkey]
      endforeach
   endif     
       
   foreach line,outlines do begin
      foreach vkey,velkeys do $
         emlvel[vkey,line] = dblarr(mapsize[1],mapsize[2])+bad
      foreach fkey,flxkeys do begin
         emlflx[fkey,line] = dblarr(mapsize[1],mapsize[2])
         emlflxerr[fkey,line] = dblarr(mapsize[1],mapsize[2])
      endforeach

      for i=0,mapsize[1]-1 do begin
         for j=0,mapsize[2]-1 do begin
            maxflx = max(emlcvdf['flux',line,i,j,*],k)
            if k ne 0 then begin
               emlflx['ftotpk',line,i,j] = maxflx
               emlflxerr['ftotpk',line,i,j] = emlcvdf['fluxerr',line,i,j,k]
               emlvel['vpk',line,i,j] = emlcvdf['vel',k]
            endif
            iv98 = value_locate(emlcvdf['cumfluxnorm',line,i,j,*],0.02d)
            iv84 = value_locate(emlcvdf['cumfluxnorm',line,i,j,*],0.16d)
            iv50 = value_locate(emlcvdf['cumfluxnorm',line,i,j,*],0.5d)
            iv16 = value_locate(emlcvdf['cumfluxnorm',line,i,j,*],0.84d)
            iv02 = value_locate(emlcvdf['cumfluxnorm',line,i,j,*],0.98d)
            if iv98 ne -1 AND iv98 ne nvel-1 then $
               emlvel['v%98',line,i,j] = emlcvdf['vel',iv98]
            if iv84 ne -1 AND iv84 ne nvel-1 then $
               emlvel['v%84',line,i,j] = emlcvdf['vel',iv84]
            if iv16 ne -1 AND iv16 ne nvel-1 then $
               emlvel['v%16',line,i,j] = emlcvdf['vel',iv16]
            if iv02 ne -1 AND iv02 ne nvel-1 then $
               emlvel['v%02',line,i,j] = emlcvdf['vel',iv02]
            if iv50 ne -1 AND iv50 ne nvel-1 then begin
               emlvel['v%50',line,i,j] = emlcvdf['vel',iv50]
               emlflx['fv%50',line,i,j] = emlcvdf['flux',line,i,j,iv50]
               emlflxerr['fv%50',line,i,j] = emlcvdf['fluxerr',line,i,j,iv50]
               if keyword_set(sigcut) AND emlflx['fv%50',line,i,j] gt 0 then $
                  if emlflx['fv%50',line,i,j]/$
                     emlflxerr['fv%50',line,i,j] lt sigcut then begin
                     emlflx['fv%50',line,i,j] = 0d
                     emlflxerr['fv%50',line,i,j] = 0d
                  endif
            endif
            if iv50 ne -1 AND iv84 ne -1 AND iv16 ne -1 AND $
               iv50 ne nvel-1 AND iv84 ne nvel-1 AND iv16 ne nvel-1 then $
                  emlvel['vsig',line,i,j] = $
                     ((emlcvdf['vel',iv50] - emlcvdf['vel',iv84])+$
                      (emlcvdf['vel',iv16] - emlcvdf['vel',iv50]))/2d
            if keyword_set(fluxvels) then $
               foreach vel,fluxvels do begin
                  ivel = value_locate(emlcvdf['vel'],vel)
                  fkey = 'fv'+string(vel,format='(I0)')
                  emlflx[fkey,line,i,j] = emlcvdf['flux',line,i,j,ivel]
                  emlflxerr[fkey,line,i,j] = emlcvdf['fluxerr',line,i,j,ivel]
                  if keyword_set(sigcut) AND emlflx[fkey,line,i,j] gt 0 then $
                     if emlflx[fkey,line,i,j]/$
                        emlflxerr[fkey,line,i,j] lt sigcut then begin
                        emlflx[fkey,line,i,j] = 0d
                        emlflxerr[fkey,line,i,j] = 0d
                     endif
               endforeach
         endfor
      endfor
   endforeach

   return,emlvel

end
