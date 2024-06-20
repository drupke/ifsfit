; docformat = 'rst'
;
;+
;
; Read single spectrum and output to structure QSOTEMPLATE in an IDL save file.
;
; :Categories:
;    IFSFIT
;
; :Returns:
; 
;    structure QSOTEMPLATE in an IDL save file
;
; :Params:
;    infits: in, required, type=string
;    outxdr: in, required, type=string
;
; :Keywords:
;    waveext: in, optional, type=integer
;    datext: in, optional, type=integer, default=1
;    varext: in, optional, type=integer, default=2
;    dqext: in, optional, type=integer, default=3
;      The extension numbers of a wavelength array, the data array, 
;      the variance array, and the dq array.
;    stelmod: in, optional, type=string
;    stelscale: in, optional, type=double, default=1d
;    stelshift: in, optional, type=integer, default=0
;      Parameters to control stellar model to subtract. STELMOD is the filename
;      of an XDR file output by IFSFA with continuum fit information. STELSCALE
;      is a multiplier to this. STELSHIFT shifts the stellar model up or down
;      in wavelength by some number of pixels.
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
;      2016dec14, DSNR, created
;      2018feb08, DSNR, added WAVEEXT keyword
;      2018may24, DSNR, added DATEXT keyword
;      2022jan05, DSNR, added VAREXT keyword; cleaned up documenation
;      2024apr05, DSNR, added option to remove stellar contribution
;
; :Copyright:
;    Copyright (C) 2016--2024 David S. N. Rupke
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
pro ifsf_makeqsotemplate,infits,outxdr,datext=datext,dqext=dqext,varext=varext,$
   waveext=waveext,stelmod=stelmod,stelscale=stelscale,stelshift=stelshift

   if ~ keyword_set(waveext) then waveext=0
   if ~ keyword_set(datext) then datext=1
   if datext eq -1 then datext=0
   if ~ keyword_set(varext) then varext=2
   if ~ keyword_set(dqext) then dqext=3
   if ~ keyword_set(stelshift) then stelshift=0

   spec = ifsr_readspec(infits,ext=datext,waveext=waveext)
   var = ifsr_readspec(infits,ext=varext)
   dq = ifsr_readspec(infits,ext=dqext)
   
   qsotemplate = {wave: double(spec[*,0]), flux: double(spec[*,1]), $
      var: double(var[*,1]), dq: dq[*,1]}

   if keyword_set(stelmod) AND keyword_set(stelscale) then begin

      restore,stelmod
      inz = where(contcube.stel_mod ne 0d)
      contcube.stel_mod /= (median(contcube.stel_mod[inz]) / $
         median(qsotemplate.flux[inz]))

      set_plot,'x'
      cgplot,qsotemplate.wave,qsotemplate.flux
      cgoplot,contcube.wave+1d,contcube.stel_mod*stelscale,color='magenta'

      stelmodshift = contcube.stel_mod
      if stelshift gt 0 then begin
         stelmodshift = [stelmodshift[0:abs(stelshift-1)],$
            stelmodshift[0:n_elements(stelmodshift)-(1+stelshift)]]
      endif else if stelshift lt 0 then begin
         stelmodshift = $
            [stelmodshift[abs(stelshift):n_elements(stelmodshift)-1],$
             stelmodshift[n_elements(stelmodshift)-abs(stelshift):n_elements(stelmodshift)-1]]
      endif

      cgoplot,qsotemplate.wave,qsotemplate.flux-stelmodshift*stelscale,$
         color='cyan'

      igd = where(qsotemplate.flux ne 0 AND qsotemplate.flux ne 1d99)
      ratio = total(stelmodshift[igd]*stelscale)/total(qsotemplate.flux[igd])
      print,'Average ratio of stellar model fraction = ',$
         string(ratio,format='(D0.2)')

      qsotemplate.flux -= stelmodshift*stelscale

   endif
   
   save,qsotemplate,file=outxdr

end
