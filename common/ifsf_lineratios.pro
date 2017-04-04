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
;    Hash of extinction values and/or line ratios, with errors:
;      OUT[keys,dx,dy]
;      keys = 'ebv', 'ebv_err'
;           = 'n2ha', 'n2ha_err'
;           = 'o1ha', 'o1ha_err'
;           = 'o3hb', 'o3hb_err'
;
; :Params:
;    flux: in, required, type=hash(lines,nx,ny)
;      Hash of line fluxes for an (nx,ny) cube.
;    fluxerr: in, required, type=hash(lines,nx,ny)
;    linelist: in, required, type=hash(nlines)
;      Hash of wavelengths, as output by IFSF_LINELIST.
;
; :Keywords:
;    ebvonly: in, optional, type=byte
;      Compute extinction values only.
;    lronly: in, optional, type=byte
;      Compute line ratios only.
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
;      2016oct03, DSNR, re-written to make input/output more transparent
;      2016oct10, DSNR, removed extinction correction from line ratios;
;                       set E(B-V) to 0 if Halpha/Hbeta too low; added
;                       [SII]/Ha
;
; :Copyright:
;    Copyright (C) 2014--2016 David S. N. Rupke
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
function ifsf_lineratios,flux,fluxerr,linelist,noerr=noerr,ebvonly=ebvonly,$
                         lronly=lronly,errlo=errlo,errhi=errhi

   bad = 1d99
   caseb = hash()
;  From Hummer & Storey 1987, T=10^4 K, n_e = 100 cm^-2
   caseb['HalphaHbeta'] = 2.86d
   caseb['HbetaHgamma'] = 2.13d
   caseb['HbetaHdelta'] = 3.86d
   loge = alog10(exp(1))

   inlines = flux.keys()
   arrsize = size(flux[inlines[0]])
   nx = arrsize[1]
   ny = arrsize[2]

   if keyword_set(noerr) then doerr = 0b else doerr = 1b

   igd = hash()
   ctgd = hash()
   foreach line,inlines do begin
      if flux.haskey(line) then $
         igd[line] = where(flux[line] gt 0 AND $
                           flux[line] ne bad AND $
                           finite(flux[line]),ctgdtmp)
      ctgd[line] = ctgdtmp
   endforeach

;   if keyword_set(sigthresh) then begin
;
;      if linesgd.haskey('Halpha') then begin
;            igdha_st = $
;               where(linemaps_ha_flux gt linemaps_ha_fluxerr*sigthresh,ctgdha_st)
;            if ctgdha_st gt 0 AND ctgdha gt 0 then begin
;               igdha = cgSETINTERSECTION(igdha,igdha_st)
;               ctgdha = n_elements(igdha)
;            endif else begin
;               igdha = -1
;               ctgdha = 0
;            endelse
;         endif
;
;         endif



;  Initialize output hashes
   out = hash()
   if doerr then begin
      errlo=hash()
      errhi=hash()
   endif

   doebv=0b
   if ((flux.haskey('Halpha') AND flux.haskey('Hbeta')) OR $
       (flux.haskey('Hbeta') AND flux.haskey('Hgamma')) OR $
       (flux.haskey('Hbeta') AND flux.haskey('Hdelta'))) AND $
       ~ keyword_set(lronly) then begin

      doebv=1b
      
      dohahb=0b
      if flux.haskey('Halpha') AND flux.haskey('Hbeta') then begin
         line1 = 'Halpha'
         line2 = 'Hbeta'
      endif else if flux.haskey('Hbeta') AND flux.haskey('Hgamma') then begin
         line1 = 'Hbeta'
         line2 = 'Hgamma'
      endif else if flux.haskey('Hbeta') AND flux.haskey('Hdelta') then begin
         line1 = 'Hbeta'
         line2 = 'Hdelta'
      endif

;     Compute E(B-V) under Case B assumptions
      igdebv = cgsetintersection(igd[line1],igd[line2])

      ebv = dblarr(nx,ny)+bad
      ebv_err = dblarr(nx,ny)+bad
      if doerr then fluxerr = [[fluxerr[line1,igdebv]],[fluxerr[line2,igdebv]]] $
      else fluxerr=0b
                                       
      ebv_tmp = ifsf_ebv_ccm([linelist[line1],linelist[line2]],$
                              [[flux[line1,igdebv]],$
                              [flux[line2,igdebv]]],$
                              caseb[line1+line2],fluxerr=fluxerr)
      if doerr then begin
         ibdebv = where(ebv_tmp[*,0] lt 0,ctbdebv)
         if ctbdebv gt 0 then begin
;            ebv_tmp[ibdebv,0] = bad
             ebv_tmp[ibdebv,0] = 0
;            ebv_tmp[ibdebv,1] = bad
         endif
         ebv[igdebv] = ebv_tmp[*,0]
         ebv_err[igdebv] = ebv_tmp[*,1]
      endif else begin
         ibdebv = where(ebv_tmp lt 0,ctbdebv)
         if ctbdebv gt 0 then $
;            ebv_tmp[ibdebv] = bad
            ebv_tmp[ibdebv] = 0
         ebv[igdebv] = ebv_tmp
      endelse

      out['ebv'] = ebv
      if doerr then errlo['ebv'] = ebv_err

   endif
   
;  Compute line ratios

   linrats = ['n2ha','o1ha','s2ha','o3hb']
   lrlines = [['[NII]6583','Halpha'],$
              ['[OI]6300','Halpha'],$
              ['[SII]6716+[SII]6731','Halpha'],$
              ['[OIII]5007','Hbeta']]

   if ~ keyword_set(ebvonly) then begin
      for i=0,n_elements(linrats)-1 do begin
         if flux.haskey(lrlines[0,i]) AND $
            flux.haskey(lrlines[1,i]) then begin
    
            igdlr = cgsetintersection(igd[lrlines[0,i]],igd[lrlines[1,i]])
            lr = dblarr(nx,ny) + bad
            lr_lin = flux[lrlines[0,i],igdlr]/$
                     flux[lrlines[1,i],igdlr]
            lr[igdlr] = alog10(lr_lin)
            out[linrats[i]] = lr
            if doerr then begin
;               lr_err = dblarr(nx,ny) + bad
;               lr_err[igdlr] = $
;                  loge*( fluxerr[lrlines[0,i],igdlr]/$
;                         flux[lrlines[0,i],igdlr] + $
;                         fluxerr[lrlines[1,i],igdlr]/$
;                         flux[lrlines[1,i],igdlr])
               lr_errlo = dblarr(nx,ny) + bad
               lr_errhi = dblarr(nx,ny) + bad
               lr_err_lin = $
                   lr_lin*$
                   sqrt((fluxerr[lrlines[0,i],igdlr]/$
                         flux[lrlines[0,i],igdlr])^2d + $
                        (fluxerr[lrlines[1,i],igdlr]/$
                         flux[lrlines[1,i],igdlr])^2d)
               lr_errlo[igdlr] = lr[igdlr] - alog10(lr_lin - lr_err_lin)
               lr_errhi[igdlr] = alog10(lr_lin + lr_err_lin) - lr[igdlr]
               errlo[linrats[i]] = lr_errlo
               errhi[linrats[i]] = lr_errhi
            endif
         endif
      endfor      
   endif
   
   return,out

end
