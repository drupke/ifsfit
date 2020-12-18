; docformat = 'rst'
;
;+
;
; Compute MPFIT-type param structure from output structure of IFSF_FITNAD,
; which is stored as structure nadcube in *.nadfit.xdr. Essentially inverts operation
; performed in IFSF_FITNAD, for input to IFSF_NADFCN.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    param array
;
; :Params:
;    nadfit: in, required, type=structure
;    ix: in, required, type=int
;    iy: in, required, type=int
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
;      2020dec18, DSNR, created
;
; :Copyright:
;    Copyright (C) 2020 David S. N. Rupke
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
function ifsf_nadfit2param,nadfit,ix,iy

   bad = 1d99

   igdhei = where(nadfit.wavehei[ix,iy,*] ne bad,nhei)
   if nhei eq -1 then nhei=0
   igdnadabs = where(nadfit.waveabs[ix,iy,*] ne bad,nnadabs)
   if nnadabs eq -1 then nnadabs=0
   igdnadem = where(nadfit.waveem[ix,iy,*] ne bad,nnadem)
   if nnadem eq -1 then nnadem=0
   param = dblarr(3 + nhei*3 + nnadabs*4 + nnadem*4)
   param[0]=nhei
   param[1]=nnadabs
   param[2]=nnadem
   if nhei gt 0 then begin
      iarr = 3 + dindgen(nhei)*3
      param[iarr] = nadfit.wavehei[ix,iy,0:nhei-1]
      param[iarr+1] = nadfit.sigmahei[ix,iy,0:nhei-1]
      param[iarr+2] = nadfit.fluxhei[ix,iy,0:nhei-1]
   endif
   if nnadabs gt 0 then begin
      iarr = 3+nhei*3 + dindgen(nnadabs)*4
      param[iarr] = nadfit.cf[ix,iy,0:nnadabs-1]
      param[iarr+1] = nadfit.tau[ix,iy,0:nnadabs-1]
      param[iarr+2] = nadfit.waveabs[ix,iy,0:nnadabs-1]
      param[iarr+3] = nadfit.sigmaabs[ix,iy,0:nnadabs-1]
   endif
   if nnadem gt 0 then begin
      iarr = 3+nhei*3+nnadabs*4 + dindgen(nnadem)*4
      param[iarr] = nadfit.waveem[ix,iy,0:nnadem-1]
      param[iarr+1] = nadfit.sigmaem[ix,iy,0:nnadem-1]
      param[iarr+2] = nadfit.flux[ix,iy,0:nnadem-1]
      param[iarr+3] = nadfit.frat[ix,iy,0:nnadem-1]
   endif
   
   return,param

end