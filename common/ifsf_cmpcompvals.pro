; docformat = 'rst'
;
;+
;
; Compute velocities for emission lines based on component velocities and sigmas
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Hash EMLCOMPVEL[...
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
;      2017feb01, DSNR, created; copied from IFSF_CMPCVDFVALS
;    
; :Copyright:
;    Copyright (C) 2017 David S. N. Rupke
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
function ifsf_cmpcompvals,emlwav,emlsig,emlflx,zref,emlwaverr=emlwaverr,$
                          emlsigerr=emlsigerr

   bad = 1d99
   c_kms = 299792.458

   comps = emlwav->keys()
   ncomps = n_elements(comps)
   outlines = emlwav['c1']->keys()
   size_cube = size(emlwav['c1',outlines[0]])

   emlcompvel = hash()

   linelist = ifsf_linelist(outlines)
   foreach c,comps do begin
      emlcompvel['vsig'+c] = hash()
      emlcompvel['v%02'+c] = hash()
      emlcompvel['v%16'+c] = hash()
      emlcompvel['v%50'+c] = hash()
      emlcompvel['v%84'+c] = hash()
      emlcompvel['v%98'+c] = hash()
      if keyword_set(emlwaverr) then emlcompvel['v%50'+c+'err'] = hash()
      if keyword_set(emlsigerr) then emlcompvel['vsig'+c+'err'] = hash()
;      if keyword_set(emlwaverr) AND keyword_set(emlsigerr) then begin
;         emlcompvel['v%02'+c+'err'] = hash()
;         emlcompvel['v%98'+c+'err'] = hash()
;      endif

      foreach line,outlines do begin

         tmpwav = emlwav[c,line]
         if keyword_set(emlwaverr) then tmpwaverr = emlwaverr[c,line]
         tmpsig = emlsig[c,line]
         if keyword_set(emlsigerr) then tmpsigerr = emlsigerr[c,line]
         tmpflx = emlflx['f'+c,line]
         igd = where(tmpwav ne 0d AND tmpwav ne bad AND $
                     tmpsig ne 0d AND tmpsig ne bad AND $
                     tmpflx ne 0d AND tmpflx ne bad,ctgd)

         if ctgd gt 0 then begin
 
            zdiff = dblarr(size_cube[1],size_cube[2])+bad
            vel = dblarr(size_cube[1],size_cube[2])+bad
            zdiff[igd] = tmpwav[igd]/(linelist[line]*(1d + zref)) - 1d
            vel[igd] = c_kms * ((zdiff[igd]+1d)^2d - 1d) / ((zdiff[igd]+1d)^2d + 1d)
            if keyword_set(emlwaverr) then begin
               velerr = dblarr(size_cube[1],size_cube[2])
               velerr[igd] = $
                  c_kms * (4d/(linelist[line]*(1d + zref))*$
                  (zdiff[igd]+1d)/((zdiff[igd]+1d)^2d + 1d)^2d) * $
                  tmpwaverr[igd]
            endif

            emlcompvel['vsig'+c,line]=dblarr(size_cube[1],size_cube[2])+bad
            emlcompvel['v%02'+c,line]=dblarr(size_cube[1],size_cube[2])+bad
            emlcompvel['v%16'+c,line]=dblarr(size_cube[1],size_cube[2])+bad
            emlcompvel['v%50'+c,line]=dblarr(size_cube[1],size_cube[2])+bad
            emlcompvel['v%84'+c,line]=dblarr(size_cube[1],size_cube[2])+bad
            emlcompvel['v%98'+c,line]=dblarr(size_cube[1],size_cube[2])+bad
;            if keyword_set(emlwaverr) AND keyword_set(emlsigerr) then begin
;              emlcompvel['v%02'+c+'err',line] = dblarr(size_cube[1],size_cube[2])
;              emlcompvel['v%98'+c+'err',line] = dblarr(size_cube[1],size_cube[2])
;            endif

            xyindices = array_indices(tmpwav,igd)
            emlcompvel['vsig'+c,line,xyindices[0,*],xyindices[1,*]]=$
               emlsig[c,line,xyindices[0,*],xyindices[1,*]]
            emlcompvel['v%50'+c,line,xyindices[0,*],xyindices[1,*]]=vel[igd]
            if keyword_set(emlwaverr) then begin
               emlcompvel['v%50'+c+'err',line] = dblarr(size_cube[1],size_cube[2])
               emlcompvel['v%50'+c+'err',line,xyindices[0,*],xyindices[1,*]]=velerr[igd]
            endif
            if keyword_set(emlsigerr) then begin
               emlcompvel['vsig'+c+'err',line] = dblarr(size_cube[1],size_cube[2])
               emlcompvel['vsig'+c+'err',line,xyindices[0,*],xyindices[1,*]]=$
                  emlsigerr[c,line,xyindices[0,*],xyindices[1,*]]
            endif
;           figure out red and blue sides for using v%02/98
            redblue02 = dblarr(ctgd) + 2d
            redblue16 = dblarr(ctgd) + 1d
            redblue84 = dblarr(ctgd) - 1d
            redblue98 = dblarr(ctgd) - 2d
            ired = where(vel[igd] ge 0)
            redblue02[ired] *= -1d
            redblue16[ired] *= -1d
            redblue84[ired] *= -1d
            redblue98[ired] *= -1d
            emlcompvel['v%02'+c,line,xyindices[0,*],xyindices[1,*]]=$
               vel[igd] + redblue02*emlsig[c,line,xyindices[0,*],xyindices[1,*]]
            emlcompvel['v%16'+c,line,xyindices[0,*],xyindices[1,*]]=$
               vel[igd] + redblue16*emlsig[c,line,xyindices[0,*],xyindices[1,*]]
            emlcompvel['v%84'+c,line,xyindices[0,*],xyindices[1,*]]=$
               vel[igd] + redblue84*emlsig[c,line,xyindices[0,*],xyindices[1,*]]
            emlcompvel['v%98'+c,line,xyindices[0,*],xyindices[1,*]]=$
               vel[igd] + redblue98*emlsig[c,line,xyindices[0,*],xyindices[1,*]]

         endif
      endforeach

   endforeach

   return,emlcompvel

end
