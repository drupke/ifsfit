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
;    Hash of extinction values and/or line ratios:
;      OUT[keys,dx,dy]
;      keys = 'ebv', 'o3hb', etc.
;    Now will also return generic (non-VO87) line ratio.
;
; :Params:
;    flux: in, required, type=hash(lines,nx,ny)
;      Hash of line fluxes for an (nx,ny) cube.
;    fluxerr: in, required, type=hash(lines,nx,ny)
;    linelist: in, required, type=hash(nlines)
;      Hash of wavelengths, as output by IFSF_LINELIST.
;
; :Keywords:
;    caseb: in, optional, type=hash
;      Hash with key-value pairs that are concatenated H line names and case B
;      flux ratios for those line names.
;    ebvonly: in, optional, type=byte
;      Compute extinction values only.
;    errlin: out, optional, type=hash
;    errlo: out, optional, type=hash
;    errhi: out, optional, type=hash
;      Output linear and log errors in quantities. (E(B-V) errors are all the same.)
;    fcnebv: in, optional, type=string, default='ifsf_ebv_ccm'
;      Function for computing selective extinction / colour excess
;    lrlist: in, optional, hash
;      Hash with key-value pairs that are line ratio labels and two-element
;      string arrays of line labels
;    lronly: in, optional, type=byte
;      Compute line ratios only.
;    noerr: in, optional, type=byte
;      Turns off error processing.
;    rv: in, optional, type=double, default=3.1
;      Normalized extinction value.
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
;      2021aug27, DSNR, added ability to import non-VO line ratio list
;      2021dec17, DSNR, added [OIII]/[OII] and [OII]3729/3726;
;                       optionally output linear errors
;
; :Copyright:
;    Copyright (C) 2014--2021 David S. N. Rupke
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
                         lronly=lronly,errlo=errlo,errhi=errhi,rv=rv,$
                         fcnebv=fcnebv,caseb=caseb,lrlist=lrlist,errlin=errlin


   if ~ keyword_set(rv) then rv=3.1d
   if ~ keyword_set(fcnebv) then fcnebv='ifsf_ebv_ccm'
   
   bad = 1d99
   if ~ keyword_set(caseb) then begin
      caseb = hash()
;     From Hummer & Storey 1987, T=10^4 K, n_e = 100 cm^-2
      caseb['HalphaHbeta'] = 2.86d
      caseb['HbetaHgamma'] = 2.13d
      caseb['HbetaHdelta'] = 3.86d
   endif
   
   loge = alog10(exp(1))

   inlines = flux.keys()
   arrsize = size(flux[inlines[0]])
   nx = arrsize[1]
   if arrsize[0] eq 2 then ny = arrsize[2] else ny = 1

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
      errlin=hash()
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
      if doerr then fluxerr_use = [[fluxerr[line1,igdebv]],[fluxerr[line2,igdebv]]] $
      else fluxerr_use=0b
                                       
      ebv_tmp = call_function(fcnebv,[linelist[line1],linelist[line2]],$
                              [[flux[line1,igdebv]],$
                              [flux[line2,igdebv]]],$
                              caseb[line1+line2],fluxerr=fluxerr_use,$
                              rv=rv)
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
      if doerr then begin
         errlo['ebv'] = ebv_err
         errhi['ebv'] = ebv_err
         errlin['ebv'] = ebv_err
      endif
         
   endif
   
;  Compute line ratios

   if not keyword_set(lrlist) then begin
      lrlabs = ['n2ha','o1ha','s2ha','o3hb','o3o2','s2','o2']
      lrlines = [['[NII]6583','Halpha'],$
         ['[OI]6300','Halpha'],$
         ['[SII]6716+[SII]6731','Halpha'],$
         ['[OIII]5007','Hbeta'],$
         ['[OIII]5007','[OII]3726+[OII]3729'],$
         ['[SII]6716','[SII]6731'],$
         ['[OII]3729','[OII]3726']]
      lrlist = hash()
      for i=0,n_elements(lrlabs)-1 do lrlist[lrlabs[i]]=lrlines[*,i]
   endif

   if ~ keyword_set(ebvonly) then begin
      foreach lrline, lrlist, lrlab do begin
         if flux.haskey(lrline[0]) AND flux.haskey(lrline[1]) then begin
    
            igdlr = cgsetintersection(igd[lrline[0]],igd[lrline[1]])
            lr = dblarr(nx,ny) + bad
            lr_lin = dblarr(nx,ny) + bad
            lr_lin[igdlr] = flux[lrline[0],igdlr]/$
               flux[lrline[1],igdlr]
            lr[igdlr] = alog10(lr_lin[igdlr])
            out[lrlab] = lr
            if doerr then begin
;               lr_err = dblarr(nx,ny) + bad
;               lr_err[igdlr] = $
;                  loge*( fluxerr[lrlines[0,i],igdlr]/$
;                         flux[lrlines[0,i],igdlr] + $
;                         fluxerr[lrlines[1,i],igdlr]/$
;                         flux[lrlines[1,i],igdlr])
               lr_errlo = dblarr(nx,ny) + bad
               lr_errhi = dblarr(nx,ny) + bad
               lr_errlin = dblarr(nx,ny) + bad
               lr_errlin[igdlr] = $
                   lr_lin[igdlr]*$
                   sqrt((fluxerr[lrline[0],igdlr]/$
                         flux[lrline[0],igdlr])^2d + $
                        (fluxerr[lrline[1],igdlr]/$
                         flux[lrline[1],igdlr])^2d)
               lr_errlo[igdlr] = lr[igdlr] - alog10(lr_lin[igdlr] - lr_errlin[igdlr])
               lr_errhi[igdlr] = alog10(lr_lin[igdlr] + lr_errlin[igdlr]) - lr[igdlr]
               errlo[lrlab] = lr_errlo
               errhi[lrlab] = lr_errhi
               errlin[lrlab] = lr_errlin
            endif
         endif
      endforeach   
   endif
   
   return,out

end
