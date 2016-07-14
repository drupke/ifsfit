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
;    sigthresh: in, optional, type=byte
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
;    Copyright (C) 2014-2015 David S. N. Rupke
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
   hbhg_caseb = 2.18d
   loge = alog10(exp(1))

   lines = ['Halpha','Hbeta','Hgamma','[NII]6583','[OI]6300','[OIII]5007',$
            '[SII]6716','[SII]6731']
   linesgd = hash()
   foreach line,lines do $
      if linemaps.haskey(line) then linesgd[line] = 1b

   if keyword_set(noerr) then doerr = 0b else doerr = 1b

   arrsize = 0
   if linesgd.haskey('Halpha') then begin
      linemaps_ha_flux = linemaps['Halpha',*,*,*,0]
      igdha = where(linemaps_ha_flux gt 0 AND $
        linemaps_ha_flux ne bad AND $
        finite(linemaps_ha_flux),ctgdha)
      arrsize = size(linemaps_ha_flux)
   endif
   if linesgd.haskey('Hbeta') then begin
      linemaps_hb_flux = linemaps['Hbeta',*,*,*,0]
      igdhb = where(linemaps_hb_flux gt 0 AND $
        linemaps_hb_flux ne bad AND $
        finite(linemaps_hb_flux),ctgdhb)
      if arrsize[0] eq 0 then arrsize = size(linemaps_hb_flux)
   endif
   if linesgd.haskey('Hgamma') then begin
     linemaps_hg_flux = linemaps['Hgamma',*,*,*,0]
     igdhg = where(linemaps_hg_flux gt 0 AND $
       linemaps_hg_flux ne bad AND $
       finite(linemaps_hg_flux),ctgdhg)
   endif
   if linesgd.haskey('[NII]6583') then begin
      linemaps_n2b_flux = linemaps['[NII]6583',*,*,*,0]
      igdn2b = where(linemaps_n2b_flux gt 0 AND $
        linemaps_n2b_flux ne bad AND $
        finite(linemaps_n2b_flux),ctgdn2b)
   endif
   if linesgd.haskey('[OI]6300') then begin
      linemaps_o1a_flux = linemaps['[OI]6300',*,*,*,0]
      igdo1a = where(linemaps_o1a_flux gt 0 AND $
        linemaps_o1a_flux ne bad AND $
        finite(linemaps_o1a_flux),ctgdo1a)
   endif
   if linesgd.haskey('[OIII]5007') then begin
      linemaps_o3b_flux = linemaps['[OIII]5007',*,*,*,0]   
      igdo3b = where(linemaps_o3b_flux gt 0 AND $
        linemaps_o3b_flux ne bad AND $
        finite(linemaps_o3b_flux),ctgdo3b)
   endif
   if linesgd.haskey('[SII]6716') then begin
      linemaps_s2a_flux = linemaps['[SII]6716',*,*,*,0]
      igds2a = where(linemaps_s2a_flux gt 0 AND $
        linemaps_s2a_flux ne bad AND $
        finite(linemaps_s2a_flux),ctgds2a)
   endif
   if linesgd.haskey('[SII]6731') then begin
      linemaps_s2b_flux = linemaps['[SII]6731',*,*,*,0]
      igds2b = where(linemaps_s2b_flux gt 0 AND $
        linemaps_s2b_flux ne bad AND $
        finite(linemaps_s2b_flux),ctgds2b)
   endif

      
   if doerr then begin

      if linesgd.haskey('Halpha') then $
         linemaps_ha_fluxerr = linemaps['Halpha',*,*,*,1]
      if linesgd.haskey('Hbeta') then $
         linemaps_hb_fluxerr = linemaps['Hbeta',*,*,*,1]
       if linesgd.haskey('Hgamma') then $
         linemaps_hg_fluxerr = linemaps['Hgamma',*,*,*,1]
      if linesgd.haskey('[NII]6583') then $
         linemaps_n2b_fluxerr = linemaps['[NII]6583',*,*,*,1]
      if linesgd.haskey('[OI]6300') then $
         linemaps_o1a_fluxerr = linemaps['[OI]6300',*,*,*,1]
      if linesgd.haskey('[OIII]5007') then $
         linemaps_o3b_fluxerr = linemaps['[OIII]5007',*,*,*,1]
      if linesgd.haskey('[SII]6716') then $
         linemaps_s2a_fluxerr = linemaps['[SII]6716',*,*,*,1]
      if linesgd.haskey('[SII]6731') then $
         linemaps_s2b_fluxerr = linemaps['[SII]6731',*,*,*,1]

      if keyword_set(sigthresh) then begin

         if linesgd.haskey('Halpha') then begin
            igdha_st = $
               where(linemaps_ha_flux gt linemaps_ha_fluxerr*sigthresh,ctgdha_st)
            if ctgdha_st gt 0 AND ctgdha gt 0 then begin
               igdha = cgSETINTERSECTION(igdha,igdha_st)
               ctgdha = n_elements(igdha)
            endif else begin
               igdha = -1
               ctgdha = 0
            endelse
         endif

         if linesgd.haskey('Hbeta') then begin
            igdhb_st = $
               where(linemaps_hb_flux gt linemaps_hb_fluxerr*sigthresh,ctgdhb_st)
            if ctgdhb_st gt 0 AND ctgdhb gt 0 then begin
               igdhb = cgSETINTERSECTION(igdhb,igdhb_st)
               ctgdhb = n_elements(igdhb)
            endif else begin
               igdhb = -1
               ctgdhb = 0
            endelse
         endif

         if linesgd.haskey('Hgamma') then begin
           igdhg_st = $
             where(linemaps_hg_flux gt linemaps_hg_fluxerr*sigthresh,ctgdhg_st)
           if ctgdhg_st gt 0 AND ctgdhg gt 0 then begin
             igdhg = cgSETINTERSECTION(igdhg,igdhg_st)
             ctgdhg = n_elements(igdhg)
           endif else begin
             igdhg = -1
             ctgdhg = 0
           endelse
         endif

         if linesgd.haskey('[NII]6583') then begin
            igdn2b_st = $
               where(linemaps_n2b_flux gt linemaps_n2b_fluxerr*sigthresh,ctgdn2b_st)
            if ctgdn2b_st gt 0 AND ctgdn2b gt 0 then begin
               igdn2b = cgSETINTERSECTION(igdn2b,igdn2b_st)
               ctgdn2b = n_elements(igdn2b)
            endif else begin
               igdn2b = -1
               ctgdn2b = 0
            endelse
         endif

         if linesgd.haskey('[OI]6300') then begin
            igdo1b_st = $
               where(linemaps_o1b_flux gt linemaps_o1b_fluxerr*sigthresh,ctgdo1b_st)
            if ctgdo1b_st gt 0 AND ctgdo1b gt 0 then begin
               igdo1b = cgSETINTERSECTION(igdo1b,igdo1b_st)
               ctgdo1b = n_elements(igdo1b)
            endif else begin
               igdo1b = -1
               ctgdo1b = 0
            endelse
         endif

         if linesgd.haskey('[OIII]5007') then begin
            igdo3b_st = $
               where(linemaps_o3b_flux gt linemaps_o3b_fluxerr*sigthresh,ctgdo3b_st)
            if ctgdo3b_st gt 0 AND ctgdo3b gt 0 then begin
               igdo3b = cgSETINTERSECTION(igdo3b,igdo3b_st)
               ctgdo3b = n_elements(igdo3b)
            endif else begin
               igdo3b = -1
               ctgdo3b = 0
            endelse
         endif

         if linesgd.haskey('[SII]6716') then begin
            igds2a_st = $
               where(linemaps_s2a_flux gt linemaps_s2a_fluxerr*sigthresh,ctgds2a_st)
            if ctgds2a_st gt 0 AND ctgds2a gt 0 then begin
               igds2a = cgSETINTERSECTION(igds2a,igds2a_st)
               ctgds2a = n_elements(igds2a)
            endif else begin
               igds2a = -1
               ctgds2a = 0
            endelse
         endif

         if linesgd.haskey('[SII]6731') then begin
            igds2b_st = $
               where(linemaps_s2b_flux gt linemaps_s2b_fluxerr*sigthresh,ctgds2b_st)
            if ctgds2b_st gt 0 AND ctgds2b gt 0 then begin
               igds2b = cgSETINTERSECTION(igds2b,igds2b_st)
               ctgds2b = n_elements(igds2b)
            endif else begin
               igds2b = -1
               ctgds2b = 0
            endelse
         endif

      endif

   endif

   nx = arrsize[1]
   ny = arrsize[2]
   ncomp = arrsize[3]

;  Initialize output hash
   out = hash()

   doebv=0b
   if ((linesgd.haskey('Halpha') AND linesgd.haskey('Hbeta')) OR $
       (linesgd.haskey('Hbeta') AND linesgd.haskey('Hgamma'))) then begin

      doebv=1b
      
      dohahb=0b
      if linesgd.haskey('Halpha') AND linesgd.haskey('Hbeta') then dohahb=1b

;     Compute E(B-V) under Case B assumptions
      if dohahb then igdebv = cgsetintersection(igdha,igdhb) $
      else igdebv = cgsetintersection(igdhb,igdhg)

      ebv = dblarr(nx,ny,ncomp)+bad
      ebv_err = dblarr(nx,ny,ncomp)+bad
      if doerr then begin
         if dohahb then fluxerr = [[linemaps_ha_fluxerr[igdebv]],$
                                   [linemaps_hb_fluxerr[igdebv]]] $
         else           fluxerr = [[linemaps_hb_fluxerr[igdebv]],$
                                   [linemaps_hg_fluxerr[igdebv]]]
      endif else fluxerr=0b
                                       
      if dohahb then $
         ebv_tmp = ifsf_ebv_ccm([6562.80d,4861.32d],$
                                [[linemaps_ha_flux[igdebv]],$
                                 [linemaps_hb_flux[igdebv]]],$
                                hahb_caseb,fluxerr=fluxerr) $
      else $
         ebv_tmp = ifsf_ebv_ccm([4861.32d,4340.46d],$
                                [[linemaps_hb_flux[igdebv]],$
                                 [linemaps_hg_flux[igdebv]]],$
                                hbhg_caseb,fluxerr=fluxerr)
      if doerr then begin
         ebv[igdebv] = ebv_tmp[*,0]
         ebv_err[igdebv] = ebv_tmp[*,1]
      endif else begin
         ebv[igdebv] = ebv_tmp
      endelse

      out['ebv'] = ebv
      if doerr then out['ebv_err'] = ebv_err

   endif
   
;  Compute de-reddened line ratios

;  [NII]/Halpha
   if linesgd.haskey('Halpha') AND linesgd.haskey('[NII]6583') then begin
    
      if doebv then begin
         igdhan2 = cgsetintersection(igdebv,igdn2b)
         ebv_use = ebv[igdhan2]
         nocorr=0b
      endif else begin
         igdhan2 = cgsetintersection(igdha,igdn2b)
         ebv_use = 0d
         nocorr=1b
      endelse
      if doerr then begin
         fluxerr = linemaps_ha_fluxerr[igdhan2]
         if doebv then ebverr = ebv_err[igdhan2] else ebverr=0d
      endif else begin
         fluxerr = 0b
         ebverr = 0b
      endelse
      dustcor_ha = $
         ifsf_dustcor_ccm(linelist['Halpha'],linemaps_ha_flux[igdhan2],$
                          ebv_use,fluxerr=fluxerr,$
                          ebverr=ebverr,nocorr=nocorr)
      if doerr then begin
         fluxerr = linemaps_n2b_fluxerr[igdhan2]
         if doebv then ebverr = ebv_err[igdhan2] else ebverr=0d
      endif else begin
         fluxerr = 0b
         ebverr = 0b
      endelse
      dustcor_n2b = $
         ifsf_dustcor_ccm(linelist['[NII]6583'],linemaps_n2b_flux[igdhan2],$
                          ebv_use,fluxerr=fluxerr,$
                          ebverr=ebverr,nocorr=nocorr)
      n2ha = dblarr(nx,ny,ncomp) + bad
      if doerr then begin
         n2ha_err = dblarr(nx,ny,ncomp) + bad
         n2ha[igdhan2] = alog10(dustcor_n2b[*,0]/dustcor_ha[*,0])
         n2ha_err[igdhan2] = loge*( dustcor_n2b[*,1]/dustcor_n2b[*,0] + $
                                dustcor_ha[*,1]/dustcor_ha[*,0] )
      endif else begin
         n2ha[igdhan2] = alog10(dustcor_n2b/dustcor_ha) 
      endelse

      out['n2ha'] = n2ha
      if doerr then out['n2ha_err'] = n2ha_err
      
   endif

;  [OI]/Halpha
   if linesgd.haskey('Halpha') AND linesgd.haskey('[OI]6300') then begin

      if doebv then begin
         igdhao1 = cgsetintersection(igdebv,igdo1a)
         ebv_use = ebv[igdhao1]
         nocorr=0b
      endif else begin
         igdhao1 = cgsetintersection(igdha,igdo1a)
         ebv_use = 0d
         nocorr=1b
      endelse
      if doerr then begin
         fluxerr = linemaps_ha_fluxerr[igdhao1]
         if doebv then ebverr = ebv_err[igdhao1] else ebverr=0d
      endif else begin
         fluxerr = 0b
         ebverr = 0b
      endelse
      dustcor_ha = $
         ifsf_dustcor_ccm(linelist['Halpha'],linemaps_ha_flux[igdhao1],$
                          ebv_use,fluxerr=fluxerr,$
                          ebverr=ebverr,nocorr=nocorr)
      if doerr then begin
         fluxerr = linemaps_o1a_fluxerr[igdhao1]
         if doebv then ebverr = ebv_err[igdhao1] else ebverr=0d
      endif else begin
         fluxerr = 0b
         ebverr = 0b
      endelse
      dustcor_o1a = $
         ifsf_dustcor_ccm(linelist['[OI]6300'],linemaps_o1a_flux[igdhao1],$
                          ebv_use,fluxerr=fluxerr,$
                          ebverr=ebverr,nocorr=nocorr)

      o1ha = dblarr(nx,ny,ncomp) + bad
      if doerr then begin
         o1ha_err = dblarr(nx,ny,ncomp) + bad
         o1ha[igdhao1] = alog10(dustcor_o1a[*,0]/dustcor_ha[*,0])
         o1ha_err[igdhao1] = loge*( dustcor_o1a[*,1]/dustcor_o1a[*,0] + $
                             dustcor_ha[*,1]/dustcor_ha[*,0] )
      endif else begin
         o1ha[igdhao1] = alog10(dustcor_o1a/dustcor_ha) 
      endelse

      out['o1ha'] = o1ha
      if doerr then out['o1ha_err'] = o1ha_err
      
   endif

;  [OIII]/Hbeta
   if linesgd.haskey('Hbeta') AND linesgd.haskey('[OIII]5007') then begin
    
      if doebv then begin
         igdhbo3 = cgsetintersection(igdebv,igdo3b)
         ebv_use = ebv[igdhbo3]
         nocorr=0b
      endif else begin
         igdhbo3 = cgsetintersection(igdhb,igdo3b)
         ebv_use = 0b
         nocorr=1b
      endelse
      if doerr then begin
         fluxerr = linemaps_hb_fluxerr[igdhbo3]
         if doebv then ebverr = ebv_err[igdhbo3] else ebverr=0d
      endif else begin
         fluxerr = 0b
         ebverr = 0b
      endelse
      dustcor_hb = $
         ifsf_dustcor_ccm(linelist['Hbeta'],linemaps_hb_flux[igdhbo3],$
                          ebv_use,fluxerr=fluxerr,$
                          ebverr=ebverr,nocorr=nocorr)
      if doerr then begin
         fluxerr = linemaps_o3b_fluxerr[igdhbo3]
         if doebv then ebverr = ebv_err[igdhbo3] else ebverr=0d
      endif else begin
         fluxerr = 0b
         ebverr = 0b
      endelse
      dustcor_o3b = $
         ifsf_dustcor_ccm(linelist['[OIII]5007'],linemaps_o3b_flux[igdhbo3],$
                          ebv_use,fluxerr=fluxerr,$
                          ebverr=ebverr,nocorr=nocorr)
      o3hb = dblarr(nx,ny,ncomp) + bad
      if doerr then begin
         o3hb_err = dblarr(nx,ny,ncomp) + bad
         o3hb[igdhbo3] = alog10(dustcor_o3b[*,0]/dustcor_hb[*,0])
         o3hb_err[igdhbo3] = loge*( dustcor_o3b[*,1]/dustcor_o3b[*,0] + $
                             dustcor_hb[*,1]/dustcor_hb[*,0] )
      endif else begin
         o3hb[igdhbo3] = alog10(dustcor_o3b/dustcor_hb)      
      endelse

      out['o3hb'] = o3hb
      if doerr then out['o3hb_err'] = o3hb_err
      
   endif
                
   return,out

end
