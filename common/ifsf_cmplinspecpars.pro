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
;    Structure of velocities and fluxes at those velocities, with each tag
;    being an (Ncols x Nrows) array. See end of file for list of tags.
;
; :Params:
;    linmap: in, required, type=structure
;      Output line spectra from IFSF_CMPLINSPECMAPS.
;
; :Keywords:
;    linpararr: in, optional, type=dblarr(Ncols,Nrows,5)
;      Option to output parameter fluxes in array form, for input to line ratio
;      routine.
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
;      2015may14, DSNR, added option to output in array form as well
;      2015sep05, DSNR, added calculation of outflow properties
;    
; :Copyright:
;    Copyright (C) 2014-2015 David S. N. Rupke
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
function ifsf_cmplinspecpars,linmap,linpararr=linpararr,ofpars=ofpars,$
                             ofthresh=ofthresh,sigthresh=sigthresh,$
                             ofignore=ofignore,diskrot=diskrot,$
                             diffthresh=diffthresh

   bad = 1d99
   c = 299792.458
   
   if ~ keyword_set(ofthresh) then ofthresh=150d
   if ~ keyword_set(sigthresh) then sigthresh=150d
   if ~ keyword_set(diffthresh) then diffthresh=50d

   mapsize = size(linmap.flux)
   nvel = n_elements(linmap.vel)

   ofignore_map = bytarr(mapsize[1],mapsize[2])
   if keyword_set(ofignore) then begin
      size_tmp = size(ofignore)
      for i=0,size_tmp[2]-1 do begin
         for j=ofignore[0,i],ofignore[2,i] do begin
            for k=ofignore[1,i],ofignore[3,i] do begin
               ofignore_map[j,k]=1b
            endfor
         endfor
      endfor
   endif

   vpk = dblarr(mapsize[1],mapsize[2])+bad
   v98 = dblarr(mapsize[1],mapsize[2])+bad
   v84 = dblarr(mapsize[1],mapsize[2])+bad
   v50 = dblarr(mapsize[1],mapsize[2])+bad
   v16 = dblarr(mapsize[1],mapsize[2])+bad
   v02 = dblarr(mapsize[1],mapsize[2])+bad
   sig = dblarr(mapsize[1],mapsize[2])+bad
   fpk = dblarr(mapsize[1],mapsize[2])
   fv98 = dblarr(mapsize[1],mapsize[2])
   fv84 = dblarr(mapsize[1],mapsize[2])
   fv50 = dblarr(mapsize[1],mapsize[2])
   fv16 = dblarr(mapsize[1],mapsize[2])
   fv02 = dblarr(mapsize[1],mapsize[2])
   if keyword_set(ofpars) then begin
      v98of = dblarr(mapsize[1],mapsize[2])+bad
      v84of = dblarr(mapsize[1],mapsize[2])+bad
      v50of = dblarr(mapsize[1],mapsize[2])+bad
      sigof = dblarr(mapsize[1],mapsize[2])+bad
      fof = dblarr(mapsize[1],mapsize[2])+bad  
      v98ofr = dblarr(mapsize[1],mapsize[2])+bad
      v84ofr = dblarr(mapsize[1],mapsize[2])+bad
      v50ofr = dblarr(mapsize[1],mapsize[2])+bad
      sigofr = dblarr(mapsize[1],mapsize[2])+bad
      fofr = dblarr(mapsize[1],mapsize[2])+bad  
      v98ofdiff = dblarr(mapsize[1],mapsize[2])+bad
      v84ofdiff = dblarr(mapsize[1],mapsize[2])+bad
      v50ofdiff = dblarr(mapsize[1],mapsize[2])+bad
      sigofdiff = dblarr(mapsize[1],mapsize[2])+bad
      fofdiff = dblarr(mapsize[1],mapsize[2])+bad  
   endif
   for i=0,mapsize[1]-1 do begin
      for j=0,mapsize[2]-1 do begin
         fpk[i,j] = max(linmap.flux[i,j,*],k)
         if k ne 0 then vpk[i,j] = linmap.vel[k]
         if ~ keyword_set(diskrot) then diskrot_spax=vpk[i,j] $
         else diskrot_spax=diskrot[i,j]      
         iv98 = value_locate(linmap.cumfluxnorm[i,j,*],0.02)
         iv84 = value_locate(linmap.cumfluxnorm[i,j,*],0.16)
         iv50 = value_locate(linmap.cumfluxnorm[i,j,*],0.5)
         iv16 = value_locate(linmap.cumfluxnorm[i,j,*],0.84)
         iv02 = value_locate(linmap.cumfluxnorm[i,j,*],0.98)
         if iv98 ne -1 then begin
            v98[i,j] = linmap.vel[iv98]
            fv98[i,j] = linmap.flux[i,j,iv98]
         endif
         if iv84 ne -1 then begin
            v84[i,j] = linmap.vel[iv84]
            fv84[i,j] = linmap.flux[i,j,iv84]
         endif
         if iv50 ne -1 then begin
            v50[i,j] = linmap.vel[iv50]
            fv50[i,j] = linmap.flux[i,j,iv50]
            if iv84 ne -1 AND iv16 ne nvel-1 then $
               sig[i,j] = ((linmap.vel[iv50] - linmap.vel[iv84])+$
                           (linmap.vel[iv16] - linmap.vel[iv50]))/2d
         endif
         if iv16 ne -1 then begin
            v16[i,j] = linmap.vel[iv16]
            fv16[i,j] = linmap.flux[i,j,iv16]
         endif
         if iv02 ne -1 then begin
            v02[i,j] = linmap.vel[iv02]
            fv02[i,j] = linmap.flux[i,j,iv02]
         endif
         if keyword_set(ofpars) AND $
            sig[i,j] gt sigthresh AND $
            not ofignore_map[i,j] then begin
            if k ne 0 then begin
;              lower and upper limits to outflow velocity
               ul_vof = diskrot_spax - ofthresh
               ll_vof = diskrot_spax + ofthresh
               i_ul_vof = value_locate(linmap.vel,ul_vof)
               i_ll_vof = value_locate(linmap.vel,ll_vof)
               ul_cumfluxnorm = linmap.cumfluxnorm[i,j,i_ul_vof]
               ll_cumfluxnorm = linmap.cumfluxnorm[i,j,i_ll_vof]
;              renormalize v50/v84/v98 w.r.t. outflow portion of CVDF
               iv16of = value_locate(linmap.cumfluxnorm[i,j,*],$
                                     0.84*ul_cumfluxnorm)
               iv50of = value_locate(linmap.cumfluxnorm[i,j,*],$
                                     0.5*ul_cumfluxnorm)
               iv84of = value_locate(linmap.cumfluxnorm[i,j,*],$
                                     0.16*ul_cumfluxnorm)
               iv98of = value_locate(linmap.cumfluxnorm[i,j,*],$
                                     0.02*ul_cumfluxnorm)
               iv16ofr = $
                  value_locate(linmap.cumfluxnorm[i,j,*],$
                               ll_cumfluxnorm+0.16*(1-ll_cumfluxnorm))
               iv50ofr = $
                  value_locate(linmap.cumfluxnorm[i,j,*],$
                               ll_cumfluxnorm+0.5*(1-ll_cumfluxnorm))
               iv84ofr = $
                  value_locate(linmap.cumfluxnorm[i,j,*],$
                               ll_cumfluxnorm+0.84*(1-ll_cumfluxnorm))
               iv98ofr = $
                  value_locate(linmap.cumfluxnorm[i,j,*],$
                               ll_cumfluxnorm+0.98*(1-ll_cumfluxnorm))
               if iv98of ne -1 then v98of[i,j] = linmap.vel[iv98of]-diskrot_spax
               if iv84of ne -1 then v84of[i,j] = linmap.vel[iv84of]-diskrot_spax
               if iv50of ne -1 then begin
                  v50of[i,j] = linmap.vel[iv50of]-diskrot_spax
                  if iv84of ne -1 AND iv16of ne -1 then $
                     sigof[i,j] = ((linmap.vel[iv50of]-linmap.vel[iv84of]) + $
                                   (linmap.vel[iv16of]-linmap.vel[iv50of]))/2d
               endif
               if iv98ofr ne nvel-1 then $
                  v98ofr[i,j] = linmap.vel[iv98ofr]-diskrot_spax
               if iv84ofr ne nvel-1 then $
                  v84ofr[i,j] = linmap.vel[iv84ofr]-diskrot_spax
               if iv50ofr ne nvel-1 then begin
                  v50ofr[i,j] = linmap.vel[iv50ofr]-diskrot_spax
                  if iv84ofr ne nvel-1 AND iv16ofr ne nvel-1 then $
                     sigofr[i,j] = ((linmap.vel[iv84ofr]-linmap.vel[iv50ofr]) + $
                                    (linmap.vel[iv50ofr]-linmap.vel[iv16ofr]))/2d
               endif
;              total emission line flux in outflow; integrate using 
;              trapezoidal rule
               fof[i,j] = total((linmap.flux[i,j,0:i_ul_vof-1] + $
                                linmap.flux[i,j,1:i_ul_vof])/2d * $
                                (linmap.wave[1:i_ul_vof] - $
                                 linmap.wave[0:i_ul_vof-1]))
               fofr[i,j] = total((linmap.flux[i,j,i_ll_vof:nvel-1] + $
                                 linmap.flux[i,j,i_ll_vof:nvel-1])/2d * $
                                 (linmap.wave[i_ll_vof:nvel-1] - $
                                  linmap.wave[i_ll_vof-1:nvel-2]))
;              blue/red difference
               if iv50of ne -1 AND iv50ofr ne nvel-1 then $
                  v50ofdiff[i,j] = abs(v50ofr[i,j])-abs(v50of[i,j])
               if iv84of ne -1 AND iv84ofr ne nvel-1 then $
                  v84ofdiff[i,j] = abs(v84ofr[i,j])-abs(v84of[i,j])
               if iv98of ne -1 AND iv98ofr ne nvel-1 then $
                  v98ofdiff[i,j] = abs(v98ofr[i,j])-abs(v98of[i,j])
               sigofdiff[i,j] = sigofr[i,j]-sigof[i,j]
               fofdiff[i,j] = fofr[i,j]-fof[i,j]
;              reject outflow if v50diff is too small
               if abs(v50ofdiff[i,j]) lt diffthresh then begin
                  v50of[i,j]=bad
                  v84of[i,j]=bad
                  v98of[i,j]=bad
                  sigof[i,j]=bad
                  fof[i,j]=bad
                  v50ofr[i,j]=bad
                  v84ofr[i,j]=bad
                  v98ofr[i,j]=bad
                  sigofr[i,j]=bad
                  fofr[i,j]=bad
                  v50ofdiff[i,j]=bad
                  v84ofdiff[i,j]=bad
                  v98ofdiff[i,j]=bad
                  sigofdiff[i,j]=bad
                  fofdiff[i,j]=bad
               endif
            endif
         endif
      endfor
   endfor
;  Sum total fluxes over components
   ftot = total(linmap.flux,3)   

   if keyword_set(linpararr) then $
      linpararr = [[[ftot]],[[fpk]],[[fv50]],[[fv84]],[[fv98]]]

   if keyword_set(ofpars) then $
      ofpars = {v50:v50of,v84:v84of,v98:v98of,sig:sigof,f:fof,$
                v50r:v50ofr,v84r:v84ofr,v98r:v98ofr,sigr:sigofr,fr:fofr,$
                v50diff:v50ofdiff,v84diff:v84ofdiff,v98diff:v98ofdiff,$
                sigdiff:sigofdiff,fdiff:fofdiff}

   return,{vpk:vpk,v98:v98,v84:v84,v50:v50,v16:v16,v02:v02,sig:sig,$
           fpk:fpk,fv98:fv98,fv84:fv84,fv50:fv50,fv16:fv16,fv02:fv02,ftot:ftot}

end
