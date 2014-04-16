; docformat = 'rst'
;
;+
;
; Find appropriate number of color bar divisions, and adjust zrange array to 
; correspond to round numbers.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Number of divisions
;
; :Params:
;    zran: in, required, type=dblarr(2)
;      Initial z range.
;    zrandiv: in, required, type=double
;      Smallest separation in z space between divisions.
;    ncbdivmax: in, required, type=integer
;      Maximum number of colorbar divisions.
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
;      2014apr15, DSNR, created
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
function ifsf_cbdiv,zran,zrandiv,ncbdivmax

;  Make steps in zrandiv equal to half of it.
   deltazrandiv = zrandiv / 2d

;  get round numbers for ranges
   zran /= zrandiv
   zran=[floor(zran[0]),ceil(zran[1])]
;  Make sure (normalized) zran is not a prime; otherwise may not be able to 
;  subdivide it (depending on how deltazrandiv is defined above).
   lowprimes = primes(15)
   lowprimescut = lowprimes(where(lowprimes gt ncbdivmax AND lowprimes lt 30))
   if where(zran[1]-zran[0] eq lowprimescut) gt -1 then zran[1] += 1
   zran *= zrandiv
   dzran = zran[1]-zran[0]

;  get round number for colorbar divisions
   if dzran eq zrandiv then ncbdiv = 1 $
   else if dzran eq zrandiv*2d then ncbdiv = 2 $
   else begin
      zrandiv -= deltazrandiv
      ncbdiv = ncbdivmax+1
;     Outer loop steps through z range divisions to look for a factor
      while ncbdiv gt ncbdivmax AND zrandiv lt dzran/2d do begin
         ncbdiv = 2
         dzdivrem = 1
         zrandiv += deltazrandiv
;        Inner loop steps through number of divisions to look for a factor
         while dzdivrem ne 0 AND ncbdiv lt ncbdivmax do begin
            ncbdiv++
            dzdiv = dzran / double(ncbdiv) / zrandiv
            dzdivrem = dzdiv - floor(dzdiv)
         endwhile
      endwhile
   endelse

   return,ncbdiv
   
end