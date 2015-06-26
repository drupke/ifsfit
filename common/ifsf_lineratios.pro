; docformat = 'rst'
;
;+
;
; Compute emission line ratios.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Structure of extinction values and extinction-corrected line ratios.
;
; :Params:
;    linemaps: in, required, type=hash(lines,nx,ny,ncomp,4)
;      Hash of line fluxes for an (nx,ny) cube, with ncomp velocity components.
;      The first plane in the last dimension is flux, then flux error,
;      wavelength, and sigma.
;    linelist: in, required, type=hash(nlines)
;      Hash of wavelengths, as output by IFSF_LINELIST.
;
; :Keywords:
;    noerr: in, optional, type=byte
;      Turns off error processing.
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
;      2014apr17, DSNR, created
;      2015may14, DSNR, added option to not process errors
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
function ifsf_lineratios,linemaps,linelist,noerr=noerr

   bad = 1d99
   hahb_caseb = 2.86d
   loge = alog10(exp(1))

   lines = ['Halpha','Hbeta','[NII]6583','[OI]6300','[OIII]5007',$
            '[SII]6716','[SII]6731']

   if keyword_set(noerr) then doerr = 0b else doerr = 1b

   linemaps_ha_flux = linemaps['Halpha',*,*,*,0]
   linemaps_hb_flux = linemaps['Hbeta',*,*,*,0]
   linemaps_n2b_flux = linemaps['[NII]6583',*,*,*,0]
   linemaps_o1a_flux = linemaps['[OI]6300',*,*,*,0]
   linemaps_o3b_flux = linemaps['[OIII]5007',*,*,*,0]
   linemaps_s2a_flux = linemaps['[SII]6716',*,*,*,0]
   linemaps_s2b_flux = linemaps['[SII]6731',*,*,*,0]
      
   if doerr then begin
      linemaps_ha_fluxerr = linemaps['Halpha',*,*,*,1]
      linemaps_hb_fluxerr = linemaps['Hbeta',*,*,*,1]
      linemaps_n2b_fluxerr = linemaps['[NII]6583',*,*,*,1]
      linemaps_o1a_fluxerr = linemaps['[OI]6300',*,*,*,1]
      linemaps_o3b_fluxerr = linemaps['[OIII]5007',*,*,*,1]
      linemaps_s2a_fluxerr = linemaps['[SII]6716',*,*,*,1]
      linemaps_s2b_fluxerr = linemaps['[SII]6731',*,*,*,1]
   endif

   igdha = where(linemaps_ha_flux gt 0 AND $
                 linemaps_ha_flux ne bad AND $
                 finite(linemaps_ha_flux),ctgdha)
   igdhb = where(linemaps_hb_flux gt 0 AND $
                 linemaps_hb_flux ne bad AND $
                 finite(linemaps_hb_flux),ctgdhb)
   igdn2b = where(linemaps_n2b_flux gt 0 AND $
                  linemaps_n2b_flux ne bad AND $
                  finite(linemaps_n2b_flux),ctgdn2b)
   igdo1a = where(linemaps_o1a_flux gt 0 AND $
                  linemaps_o1a_flux ne bad AND $
                  finite(linemaps_o1a_flux),ctgdo1a)
   igdo3b = where(linemaps_o3b_flux gt 0 AND $
                  linemaps_o3b_flux ne bad AND $
                  finite(linemaps_o3b_flux),ctgdo3b)
   igds2a = where(linemaps_s2a_flux gt 0 AND $
                  linemaps_s2a_flux ne bad AND $
                  finite(linemaps_s2a_flux),ctgds2a)
   igds2b = where(linemaps_s2b_flux gt 0 AND $
                  linemaps_s2b_flux ne bad AND $
                  finite(linemaps_s2b_flux),ctgds2b)

   arrsize = size(linemaps_ha_flux)
   nx = arrsize[1]
   ny = arrsize[2]
   ncomp = arrsize[3]

;  Compute E(B-V) under Case B assumptions
   igdhahb = cgsetintersection(igdha,igdhb)

   ebv = dblarr(nx,ny,ncomp)+bad
   ebv_err = dblarr(nx,ny,ncomp)+bad
   if doerr then fluxerr = [[linemaps_ha_fluxerr[igdhahb]],$
                            [linemaps_hb_fluxerr[igdhahb]]] $
   else fluxerr=0b
                                       
   ebv_tmp = $
      ifsf_ebv_ccm([6562.80d,4861.32d],$
                   [[linemaps_ha_flux[igdhahb]],[linemaps_hb_flux[igdhahb]]],$
                   hahb_caseb,fluxerr=fluxerr)
   if doerr then begin
      ebv[igdhahb] = ebv_tmp[*,0]
      ebv_err[igdhahb] = ebv_tmp[*,1]
   endif else begin
      ebv[igdhahb] = ebv_tmp
   endelse
   
;  Compute de-reddened line ratios

;  [NII]/Halpha
   igdhahbn2b = cgsetintersection(igdhahb,igdn2b)
   if doerr then begin
      fluxerr = linemaps_ha_fluxerr[igdhahbn2b]
      ebverr = ebv_err[igdhahbn2b]
   endif else begin
      fluxerr = 0b
      ebverr = 0b
   endelse
   dustcor_ha = $
      ifsf_dustcor_ccm(linelist['Halpha'],linemaps_ha_flux[igdhahbn2b],$
                       ebv[igdhahbn2b],fluxerr=fluxerr,$
                       ebverr=ebverr)
   if doerr then begin
      fluxerr = linemaps_n2b_fluxerr[igdhahbn2b]
      ebverr = ebv_err[igdhahbn2b]
   endif else begin
      fluxerr = 0b
      ebverr = 0b
   endelse
   dustcor_n2b = $
      ifsf_dustcor_ccm(linelist['[NII]6583'],linemaps_n2b_flux[igdhahbn2b],$
                       ebv[igdhahbn2b],fluxerr=fluxerr,$
                       ebverr=ebverr)

   n2ha = dblarr(nx,ny,ncomp) + bad
   if doerr then begin
      n2ha_err = dblarr(nx,ny,ncomp) + bad
      n2ha[igdhahbn2b] = alog10(dustcor_n2b[*,0]/dustcor_ha[*,0])
      n2ha_err[igdhahbn2b] = loge*( dustcor_n2b[*,1]/dustcor_n2b[*,0] + $
                             dustcor_ha[*,1]/dustcor_ha[*,0] )
    endif else begin
      n2ha[igdhahbn2b] = alog10(dustcor_n2b/dustcor_ha) 
    endelse

;  [OIII]/Hbeta
   igdhahbo3b = cgsetintersection(igdhahb,igdo3b)
   if doerr then begin
      fluxerr = linemaps_hb_fluxerr[igdhahbo3b]
      ebverr = ebv_err[igdhahbo3b]
   endif else begin
      fluxerr = 0b
      ebverr = 0b
   endelse
   dustcor_hb = $
      ifsf_dustcor_ccm(linelist['Hbeta'],linemaps_hb_flux[igdhahbo3b],$
                       ebv[igdhahbo3b],fluxerr=fluxerr,$
                       ebverr=ebverr)
    if doerr then begin
      fluxerr = linemaps_o3b_fluxerr[igdhahbo3b]
      ebverr = ebv_err[igdhahbo3b]
   endif else begin
      fluxerr = 0b
      ebverr = 0b
   endelse
   dustcor_o3b = $
      ifsf_dustcor_ccm(linelist['[OIII]5007'],linemaps_o3b_flux[igdhahbo3b],$
                       ebv[igdhahbo3b],fluxerr=fluxerr,$
                       ebverr=ebverr)
   o3hb = dblarr(nx,ny,ncomp) + bad
   if doerr then begin
      o3hb_err = dblarr(nx,ny,ncomp) + bad
      o3hb[igdhahbo3b] = alog10(dustcor_o3b[*,0]/dustcor_hb[*,0])
      o3hb_err[igdhahbo3b] = loge*( dustcor_o3b[*,1]/dustcor_o3b[*,0] + $
                             dustcor_hb[*,1]/dustcor_hb[*,0] )
   endif else begin
      o3hb[igdhahbo3b] = alog10(dustcor_o3b/dustcor_hb)      
   endelse

   out = hash()
   out['ebv'] = ebv
   out['n2ha'] = n2ha
   out['o3hb'] = o3hb
         
   if doerr then begin
      out['ebv_err'] = ebv_err
      out['n2ha_err'] = n2ha_err
      out['o3hb_err'] = o3hb_err
   endif
                   
   return,out

end
