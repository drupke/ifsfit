; docformat = 'rst'
;
;+
;
; Find the median and 68% probability intervals above and below it (i.e., 
; distance from the mean in one direction that encompasses 34% of the prob,
; and the same in the other direction).
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Median, "low" sigma, "high" sigma.
;    
; :Params:
;    dat: in, required, type=dblarr
;      Data.
;      
; :Keywords:
;    mean: in, optional, type=byte
;      Use the mean instead of the median.
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
;      2014jul07, DSNR, created
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
function ifsf_lohisig,dat,mean=mean

   ponesig = 0.682689492137d

   n = n_elements(dat)
   if ~ keyword_set(mean) then med = median(dat) else med = mean(dat)
   sdat = dat[sort(dat)]
   isiglo = fix(double(n)/2d*(1d -ponesig))
   isighi = fix(double(n)/2d*(1d +ponesig))
   errlo = med - sdat[isiglo]
   errhi = sdat[isighi] - med
   
   return,[med,errlo,errhi]

end