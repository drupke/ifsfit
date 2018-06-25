; docformat = 'rst'
;
;+
;
; Take outputs from MC error estimation in IFSF_FITLOOP, plot results, and find
; error bars in distributions. Returns hash where each tag is a stellar parameter
; that points to a two-component area of one-sided error bars
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
function ifsf_ploterrdist,dat,pos,xlab,noerase=noerase

   if ~ keyword_set(noerase) then noerase=0 else noerase=1

   bins=[]
   cghistoplot,dat,pos=pos,missing=bad,ytit='',$
               xtit=xlab,histdata=h,noerase=noerase,$
               locations=loc,binsize=bins,xticks=2
   binc = loc + (bins/2d)
   yfit = mpfitpeak(binc,h,a,nterms=3,errors=sqrt(h))
   cgoplot,binc,yfit
   stat = ifsf_lohisig(dat)
   cgoplot,[stat[0],stat[0]],[0,max(h)]
   cgoplot,[stat[0]-stat[1],stat[0]-stat[1]],[0,max(h)],linesty=2
   cgoplot,[stat[0]+stat[2],stat[0]+stat[2]],[0,max(h)],linesty=2
   
   return,[stat[1],stat[2]]

end
;- 
function ifsf_cterrdist,mcdist,outfile

   errors = hash()

   cgps_open,outfile,charsize=1,/inches,xs=7.5,ys=7.5,/qui,/nomatch
   pos = cglayout([2,2],ixmar=[3,3],iymar=[4,4],oxmar=[0,0],oymar=[4,1],$
                  xgap=0,ygap=4,aspect=1d)

;  zstar
   key = 'zstar'
   label = 'z'
   if mcdist.haskey(key) then $
      errors[key] = ifsf_ploterrdist(mcdist[key],pos[*,0],label)

;  sigma
   key = 'ct_ppxf_sigma'
   label = '$\sigma$ (km/s)'
   if mcdist.haskey(key) then $
      errors[key] = ifsf_ploterrdist(mcdist[key],pos[*,1],label,/noerase)

;  E(B-V)
   key = 'ct_ebv'
   label = 'E(B-V)'
   if mcdist.haskey(key) then $
      errors[key] = ifsf_ploterrdist(mcdist[key],pos[*,2],label,/noerase)

   cgtext,0.5,0.99,'Stellar Properties',/norm,align=0.5
   cgps_close,/pdf,/delete_ps
   
   return,errors

end