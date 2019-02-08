; docformat = 'rst'
;
;+
;
; Plot distribution of Monte Carlo simulation results and return error estimates.
;
; :Categories:
;    IFSF
;
; :Returns:
;    None.
;
; :Params:
;    dat: in, required, type=dblarr(N)
;    pos: in, required, type=dblarr(2)
;    xlab: in, required, type=string
;
; :Keywords:
;    noerase: in, optional, type=byte
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
;      2018jun26, DSNR, created
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
function ifsf_pltmcdist,dat,_REF_EXTRA=extra

   bad=1d99

   if ~ keyword_set(noerase) then noerase=0 else noerase=1

   bins=[]
;  Check that the data are not all a single value
   ione = where(dat eq dat[0],ctone)
   if ctone eq n_elements(dat) then begin
      errs=[0d,0d]
   endif else begin
      cghistoplot,dat,missing=bad,ytit='',histdata=h,locations=loc,binsize=bins,$
                  xticks=2,_EXTRA=extra
      binc = loc + (bins/2d)
      if n_elements(locations) ge 3 then begin
         yfit = mpfitpeak(binc,h,a,nterms=3,errors=sqrt(h))
         cgoplot,binc,yfit
      endif
      stat = ifsf_lohisig(dat)
      cgoplot,[stat[0],stat[0]],[0,max(h)]
      cgoplot,[stat[0]-stat[1],stat[0]-stat[1]],[0,max(h)],linesty=2
      cgoplot,[stat[0]+stat[2],stat[0]+stat[2]],[0,max(h)],linesty=2
      errs = [stat[1],stat[2]]
   endelse

   return,errs

end
