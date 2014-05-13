; docformat = 'rst'
;
;+
;
;
; :Categories:
;    IFSFIT
;
; :Returns:
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
;      2014apr30, DSNR, copied from GMOS_CHECKCOMP; rewrote, documented, 
;                       added copyright and license
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
;-
function ifsf_checkcomp,linepars,linetie,ncomp,newncomp,siglim,$
         sigcut=sigcut

   if ~ keyword_set(sigcut) then sigcut = 3d
   
   newncomp = hash()

;  Find unique anchors
   anchors = (linetie.values()).toarray()
   sanchors = anchors[sort(anchors)]
   uanchors = sanchors[uniq(sanchors)]

;  Find lines associated with each unique anchor. NEWLINETIE is a hash of lists,
;  with the keys being the lines that are tied TO and the list corresponding to 
;  each key consisting of the tied lines.
   newlinetie = hash(uanchors)
   foreach val,newlinetie,key do newlinetie[key] = list()
   foreach val,linetie,key do newlinetie[val].add,key

;  Loop through anchors   
   foreach tiedlist,newlinetie,key do begin
      if ncomp[key] gt 0 then begin
         goodcomp = intarr(ncomp[key])
;        Loop through lines tied to each anchor, looking for good components
         foreach line,tiedlist do begin
            igd = where((linepars.fluxpk)[line,0:ncomp[line]-1] ge $
                        sigcut*(linepars.fluxpkerr)[line,0:ncomp[line]-1] AND $
                        (linepars.fluxpkerr)[line,0:ncomp[line]-1] gt 0 AND $
                        (linepars.sigma)[line,0:ncomp[line]-1] gt siglim[0],$
                        ctgd)
            if ctgd gt 0 then goodcomp[igd]++
         endforeach
;        Find number of good components
         tmpncomp = 0
         for i=0,ncomp[key]-1 do if goodcomp[i] gt 0 then tmpncomp++
         if tmpncomp ne ncomp[key] then begin
            newncomp[key]=tmpncomp
;           Loop through lines tied to each anchor and set proper number of 
;           components
            foreach line,tiedlist do ncomp[line]=tmpncomp
         endif
      endif
   endforeach

end
