; docformat = 'rst'
;
;+
;
; Sort emission line components based on flux, wavelength, or linewidth.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;   New hash of emission line properties with sorting done.
;
; :Params:
;   dx: in, required, type=double
;     Number of columns in line maps.
;   dy: in, required, type=double
;     Number of rows in line maps.
;   linmaps: in, required, type=hash
;     Hash of input maps of emission line properties [each line consists of a
;     (dx,dy,maxncomp,5) array, with the last dimension being a line property].
;   linetie: in, required, type=hash
;     Hash indicating the reference line that each emission line is tied to.
;   sortlines: in, required, type=strarr(N)
;     Array of lines to sort. Lines that are tied to a reference line do not 
;     have to be included
;   sorttype: in, required, type=hash
;     Type of sorting for each emission line in sortlines: 
;     'flux', 'wave', or 'sigma'. Each hash element is one of these strings.
;     
; :Keywords:
;   sortdir: in, optional, type=hash
;     Default sort direction is ascending. This hash contains all of the lines
;     in sortlines; if a hash element is set to 'descending,' sort is reversed.
; 
; :Author:
;    David Rupke
;
; :History:
;    Change History::
;      2015feb06, DSNR, created
;      2015jun02, DSNR, added documentation, fixed bug
;
; :Copyright:
;    Copyright (C) 2015 David S. N. Rupke
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
function ifsf_sortcomp,dx,dy,linmaps,linetie,sortlines,sorttype,sortdir=sortdir

   bad = 1d99
   sortdown=0 ; default is to sort ascending
   slinmaps = linmaps

;  Loop through emission lines
   foreach line,sortlines do begin

      if sorttype[line] eq 'flux' then ilinmap=0
      if sorttype[line] eq 'wave' then ilinmap=2
      if sorttype[line] eq 'sigma' then ilinmap=3

      if keyword_set(sortdir) then $
         if sortdir[line] eq 'descending' then sortdown = 1 else sortdown = 0

      tiedlines = linetie.where(line,count=nties)
;     Remove reference line; otherwise it gets sorted twice. [I can't figure out
;     why ...]
      tiedlines.Remove,tiedlines.where(line)
      nties--

      for i=0,dx-1 do begin
         for j=0,dy-1 do begin
            tmpmap = linmaps[line,i,j,*,ilinmap]
            maxncomp = n_elements(tmpmap)
            igd = where(tmpmap ne bad AND tmpmap ne 0d,ctgd)
            if ctgd gt 0 then begin
               isort = sort(linmaps[line,i,j,igd,ilinmap])
               if sortdown then isort = reverse(isort)
               for k=0,4 do begin
                  slinmaps[line,i,j,0:ctgd-1,k] = linmaps[line,i,j,isort,k]
                  if ctgd lt maxncomp then $
                     slinmaps[line,i,j,ctgd:maxncomp-1,k] = bad
                  if nties gt 0 then begin
                     foreach tline,tiedlines do begin
                        slinmaps[tline,i,j,0:ctgd-1,k] = linmaps[tline,i,j,isort,k]
                        if ctgd lt maxncomp then $
                           slinmaps[tline,i,j,ctgd:maxncomp-1,k] = bad
                     endforeach
                  endif
               endfor
            endif
         endfor
      endfor
      
   endforeach

   return,slinmaps
end
