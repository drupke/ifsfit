; docformat = 'rst'
;
;+
;
; Take outputs from MC error estimation in IFSF_FITLOOP, plot results, and find
; error bars in distributions. Returns hash where each tag is a stellar parameter
; that points to a two-component area of one-sided error bars.opto
;
; :Categories:
;    IFSF
;
; :Returns:
;    None.
;
; :Params:
;    mcdist: in, required, type=hash()
;      Hash of data to histogram. Each tag is a stellar parameter that points to
;      the array of simulated values of that parameter.
;    outfile: in, required, type=string
;      Name of output plot.
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
;      2018jun25, DSNR, created
;
; :Copyright:
;    Copyright (C) 2018 David S. N. Rupke
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
function ifsf_pltmcstel,mcdist,outfile

   errors = hash()

   cgps_open,outfile,charsize=1,/inches,xs=7.5,ys=7.5,/qui,/nomatch
   pos = cglayout([2,2],ixmar=[3,3],iymar=[4,4],oxmar=[0,0],oymar=[4,1],$
                  xgap=0,ygap=4,aspect=1d)

;  zstar
   key = 'zstar'
   label = 'z'
   if mcdist.haskey(key) then $
      errors[key] = ifsf_pltmcdist(mcdist[key],position=pos[*,0],xtitle=label)

;  sigma
   key = 'ct_ppxf_sigma'
   label = '$\sigma$ (km/s)'
   if mcdist.haskey(key) then $
      errors[key] = ifsf_pltmcdist(mcdist[key],position=pos[*,1],xtitle=label,$
                                   /noerase,mininput=0d,$
                                   min_value=0d,xrange=[0,max(mcdist[key])])

;  E(B-V)
   key = 'ct_ebv'
   label = 'E(B-V)'
   if mcdist.haskey(key) then $
      errors[key] = ifsf_pltmcdist(mcdist[key],position=pos[*,2],xtitle=label,$
                                   /noerase,mininput=0d,$
                                   min_value=0d,xrange=[0,max(mcdist[key])])

   cgtext,0.5,0.99,'Stellar Properties',/norm,align=0.5
   cgps_close,/pdf,/delete_ps
   
   return,errors

end