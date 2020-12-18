; docformat = 'rst'
;
;+
;
; Compute "empirical" absorption and emission-line velocities for NaD.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Array of velocities.
;
; :Params:
;    wave: in, required, type=dblarr(N)
;      Wavelengths.
;    flux: in, required, type=dblarr(N)
;      Normalized fluxes.
;    err: in, required, type=dblarr(N)
;      Flux errors.
;    inds: in, required, type=intarr(4)
;      Indices into previous arrays giving locations of lower and
;      upper boundaries of absorption and emission lines.
;    zsys: in, required, type=double
;      Systemic redshift.
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
;      2014may15, DSNR, created
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
function ifsf_cmpnadvel,wave,flux,err,inds,zsys

   bad=1d99
   c_kms = 299792.458d
   linelist = ifsf_linelist(['NaD1','NaD2','HeI5876'],/quiet)
;  Average wavelength of NaD lines
   lamavgrest = ((linelist['NaD2']+linelist['NaD1'])/2d)

   if inds[0] ne -1l AND inds[1] ne -1l AND inds[0] ne inds[1] then begin
;     Average wavelength of absorption
      abs_lamavg = (wave[inds[1]]+wave[inds[0]])/2d
;     Width of absorption in wavelength and velocity spaces
      abs_lamwid = wave[inds[1]]-wave[inds[0]]
      abs_velwid = abs_lamwid/abs_lamavg*c_kms
;     Depth-weighted average absorption wavelength and velocity
      abs_lamwavg = total(wave[inds[0]:inds[1]]*(1d -flux[inds[0]:inds[1]]))/$
                    total(1d -flux[inds[0]:inds[1]])
      zdiff = abs_lamwavg/lamavgrest-1d - zsys
      abs_velwavg = c_kms * ((zdiff+1d)^2d - 1d) / ((zdiff+1d)^2d + 1d)
;     'Maximum' velocity: distance in velocity space between zsys and bluest 
;     wavelength, assuming the bluest wavelength corresponds to NaD2
      zdiff = wave[inds[0]]/linelist['NaD2']-1d - zsys
      abs_velmax = c_kms * ((zdiff+1d)^2d - 1d) / ((zdiff+1d)^2d + 1d)
   endif else begin
      abs_velwid=bad
      abs_velwavg=bad
      abs_velmax=bad
   endelse

   if inds[2] ne -1l AND inds[3] ne -1l AND inds[2] ne inds[3] then begin
;     Average wavelength of emission
      em_lamavg = (wave[inds[3]]+wave[inds[2]])/2d
;     Width of emission in wavelength and velocity spaces
      em_lamwid = wave[inds[3]]-wave[inds[2]]
      em_velwid = em_lamwid/em_lamavg*c_kms
;     Flux-weighted average emission wavelength and velocity
      em_lamwavg = total(wave[inds[2]:inds[3]]*flux[inds[2]:inds[3]])/$
                   total(flux[inds[2]:inds[3]])
      zdiff = em_lamwavg/lamavgrest-1d - zsys
      em_velwavg = c_kms * ((zdiff+1d)^2d - 1d) / ((zdiff+1d)^2d + 1d)
;     'Maximum' velocity: distance in velocity space between flux-weighted
;     average emission wavelength and reddest wavelength, assuming the reddest 
;     wavelength corresponds to NaD1
      zdiff = wave[inds[3]]/linelist['NaD1']-1d - zsys
      em_velmax = c_kms * ((zdiff+1d)^2d - 1d) / ((zdiff+1d)^2d + 1d)
   endif else begin
      em_velwid=bad
      em_velwavg=bad
      em_velmax=bad
   endelse

   return,[abs_velwid,abs_velwavg,abs_velmax,$
           em_velwid,em_velwavg,em_velmax]

end
