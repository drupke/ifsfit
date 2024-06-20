; docformat = 'rst'
;
;+
;
; Plot emission line fit and output to JPG.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    None.
;
; :Params:
;    instr: in, required, type=structure
;      Contains results of fit.
;    pltpar: in, required, type=structure
;      Contains parameters to control plot. Tags:
;      label: type=strarr(Nlines), line labels for plot
;      wave: type=dblarr(Nlines), rest wavelengths of lines
;      linoth: type=dblarr(Notherlines,Ncomp), wavelengths of other
;              lines to plot
;      nx: type=int, # of plot columns
;      ny: type=int, # of plot rows
;      
;    outfile: in, required, type=string
;      Full path and name of output plot.
;
; :Keywords:
;    boxsmooth: in, optional, type=double
;      Boxcar smooth the data. Set keyword to box width.
;    emlsig: in, optional, type=hash
;      Print emission-line sigmas on plot.
;    micron: in, optional, type=byte
;      Label output plots in um rather than A. Input wavelengths still assumed
;      to be in A.
;    ps: in, optional, type=bool
;      Output in eps rather than jpg.
;    resid: in, optional, type=string, default='subemlmod'
;      type of residual to show
;      options: 'subemlmod', 'subtotmod'
;    title: in, optional, type=string
;    velinset: in, optional, type=list
;      Will plot an inset with velocity profiles of a line in one panel. First
;      element is line label; second is 4-element array of [x0,y0,x1,y1] for
;      inset boundaries in coordinates normalized to the data+residual panel.
;      Third element is two-element array for velocity range in km/s.
;    xmajtickint: in, optional, type=double
;    xmintickint: in, optional, type=double
;    yranscl: in, optional, type=double
;      Default is to scale y-axis using data and model. Use this multiplier to 
;      scale the y-axis by this fraction of the model max-min or by this
;      fraction of 6*mean(error).
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
;      2009, DSNR, created
;      13sep12, DSNR, re-written
;      2013oct, DSNR, documented
;      2013nov21, DSNR, renamed, added license and copyright 
;      2015may13, DSNR, switched from using LAYOUT keyword to using CGLAYOUT
;                       procedure to fix layout issues
;      2016aug31, DSNR, added overplotting of continuum ranges masked during
;                       continuum fit with thick cyan line
;      2016sep13, DSNR, added MICRON keyword
;      2022oct18, DSNR, added various capabilities, incl. VELINSET and BOXSMOOTH
;      2024mar18, DSNR, added some plot options
;    
; :Copyright:
;    Copyright (C) 2013--2024 David S. N. Rupke
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
pro ifsf_pltlin,instr,pltpar,outfile,emlz,ps=ps,emlsig=emlsig,$
   velinset=velinset,yranscl=yranscl,boxsmooth=boxsmooth,title=title,$
   xmajtickint=xmajtickint,xmintickint=xmintickint,resid=resid

   bad = 1d99
   c_kms = 299792.458

   if keyword_set(ps) then dops=1 else dops=0
   if ~ keyword_set(resid) then resid='subemlmod'

   if dops then begin
      cgps_open,filename=outfile+'.eps',/encapsulated,xsize=10.5,ysize=7.5,$
         /inches,/nomatch,charsize=1,default_thick=2
      backgcol='White'
      linecol='Black'
      colors = ['Magenta','Teal','Orange']
      pos = cglayout([pltpar.nx,pltpar.ny],$ ; ixmar=[5d,0d],iymar=[-5d,0d],$
         oxmar=[8,1],oymar=[6,3],xgap=5,ygap=3)
      linelabsize = 1.
      linelabcol = 'Blue'
   endif else begin
      set_plot,'Z'
      device,decomposed=0,set_resolution=[1575,900],set_pixel_depth=24
      !P.charsize=1.5
      !P.charthick=1
      !P.thick=4
      erase
      backgcol='Black'
      linecol='White'
      colors = ['Magenta','Green','Orange','Teal']
      linelabcol = 'White'
      pos = cglayout([pltpar.nx,pltpar.ny],$ ; ixmar=[5d,0d],iymar=[-5d,0d],$
         oxmar=[22,0],oymar=[10,2],xgap=20,ygap=6)
      linelabsize=1.5
   endelse

   defaultXtickint=!X.tickinterval
   defaultXminor=!X.minor
   if keyword_set(xmajtickint) then !X.tickinterval = xmajtickint $
   else !X.tickinterval=50
   if keyword_set(xmintickint) then !X.minor = xmintickint $
   else !X.minor=10
   if tag_exist(pltpar,'micron') then begin
      !X.tickinterval /= 0.5d4
      !X.minor /= 0.5d4
   endif
   if tag_exist(pltpar,'meter') then begin
      !X.tickinterval /= 0.5d10
      !X.minor /= 0.5d10
   endif

  maxncomp = instr.param[1]

  wave = instr.wave
  spectot = instr.spec
  specstars = instr.cont_dat
  speclines = instr.emlin_dat
  specerr = instr.spec_err
  modstars = instr.cont_fit
  modlines = instr.emlin_fit
  modtot = modstars + modlines

;  norm = max(modstars)
;  spectot /= norm
;  specstars /= norm
;  speclines /= norm
;  modtot /= norm
;  modstars /= norm
;  modlines /= norm

  zbase = instr.zstar

   if tag_exist(pltpar,'micron') then wave /= 1d4
;   if tag_exist(pltpar,'meter') then wave /= 1d10

   ; check that continuum was fit
   if n_elements(instr.ct_indx) gt 1 then begin
;  Find masked regions during continuum fit
      nmasked = 0 ; # of masked regions
;  Find consecutive unmasked regions
      consec,instr.ct_indx,lct,hct,nct
;  Set interior masked regions
      if nct gt 1 then begin
         nmasked = nct-1
         masklam = dblarr(2,nmasked)
         for i=0,nmasked-1 do $
            masklam[*,i] = [wave[instr.ct_indx[hct[i]]],$
               wave[instr.ct_indx[lct[i+1]]]]
      endif
;  Set masked region if it occurs at beginning of lambda array
      if instr.ct_indx[0] ne 0 then begin
         nmasked++
         masklam = [[wave[0],wave[instr.ct_indx[lct[0]-1]]],[masklam]]
      endif
;  Set masked region if it occurs at end of lambda array
      if instr.ct_indx[n_elements(instr.ct_indx)-1] ne $
         n_elements(instr.wave)-1 then begin
         nmasked++
         masklam = [[masklam],$
            [wave[instr.ct_indx[hct[nct-1]]],$
            wave[n_elements(instr.wave)-1]]]
      endif
   endif else nmasked=0

  nlin = n_elements(pltpar.label)
  linlab = pltpar.label
  linwav = pltpar.wave
  off = pltpar.off
  if tag_exist(pltpar,'linoth') then linoth = pltpar.linoth $
  else linoth = strarr(1,nlin)
  for i=0,nlin-1 do begin

     ; flag for velocity inset plot
     dovelinset = 0b

     linwavtmp = linwav[i]
     xran = (linwavtmp[0] + off[*,i]) * (1d + zbase)
     ind = where(wave gt xran[0] AND wave lt xran[1],ct)

     cgplot,[0],/nodata,xsty=4,ysty=4,pos=pos[*,i],noerase=i ne 0,backg=backgcol
     xwin = [pos[0,i],pos[2,i]]
     ywin = [pos[1,i],pos[3,i]]
     dxwin = xwin[1]-xwin[0]
     dywin = ywin[1]-ywin[0]

     if ct gt 0 then begin
        pos_fit = [xwin[0],ywin[0]+0.3*dywin,$
                   xwin[1],ywin[1]]
        ydat = spectot
        ymod = modtot
        if keyword_set(yranscl) then begin
           sigdat = mean(specerr[ind])
           yran = [min(ymod[ind]),max(ymod[ind])]
           dymod = max(ymod[ind])-min(ymod[ind])
           if 2d*sigdat > dymod then begin
              yran += yranscl*6d*[-sigdat,sigdat]
           endif else begin
              yran += yranscl*[-dymod,dymod]
              if yran[0] lt 0 then yran[0] = 0d
           endelse
        endif else begin
           yran = [min([ydat[ind],ymod[ind]]),max([ydat[ind],ymod[ind]])]
        endelse
        icol = double(i)/double(pltpar.nx)
        if icol eq fix(icol) then ytit = 'Fit' else ytit = ''
        if keyword_set(boxsmooth) then ydat_use = smooth(ydat,boxsmooth) $
        else ydat_use = ydat
        cgplot,wave,ydat_use,xran=xran,yran=yran,pos=pos_fit,$
               xtickn=replicate(' ',60),ytit=ytit,/noerase,$
               axiscol=linecol,col=linecol,/norm,/xsty,/ysty,thick=1
        cgoplot,wave,ymod,color='Red',thick=4
        for j=1,maxncomp do begin
           flux = ifsf_cmplin(instr,linlab[i],j,/velsig)
           cgoplot,wave,yran[0]+flux,color=colors[j-1],linesty=2,thick=2
           ; check for velocity inset plot
           if keyword_set(velinset) then $
              if velinset[0] eq linlab[i] then dovelinset=1b
           if linoth[0,i] ne '' then begin
              for k=0,n_elements(linoth[*,i])-1 do begin
                 if linoth[k,i] ne '' then begin
                    flux = ifsf_cmplin(instr,linoth[k,i],j,/velsig)
                    cgoplot,wave,yran[0]+flux,color=colors[j-1],linesty=2,$
                            thick=2
                    if keyword_set(velinset) then $
                       if velinset[0] eq linoth[k,i] then dovelinset=1b
                 endif
              endfor
           endif
        endfor
        cgtext,xran[0]+(xran[1]-xran[0])*0.05d,$
               yran[0]+(yran[1]-yran[0])*0.90d,$
               linlab[i],charsize=linelabsize,charthick=2,/dat,$
               col=linelabcol
        if linoth[0,i] ne '' then $
           for k=0,n_elements(linoth[*,i])-1 do $
              if linoth[k,i] ne '' then $
                 cgtext,xran[0]+(xran[1]-xran[0])*0.05d,$
                    yran[0]+(yran[1]-yran[0])*(0.90d - double(k+1)*0.1d),$
                    linoth[k,i],charsize=linelabsize,charthick=2,/dat,$
                    col=linelabcol
; z labels
        if emlz.haskey(linlab[i]) then $
           for j=0,maxncomp-1 do $
              if emlz[linlab[i],j] ne bad then $
                 cgtext,xran[0]+(xran[1]-xran[0])*0.7d,$
                    yran[0]+(yran[1]-yran[0])*(0.90d - double(j)*0.1d),$
                    'z$\down'+string(j+1,format='(I0)')+'$='+$
                    string(emlz[linlab[i],j],format='(D0.4)'),$
                    charsize=linelabsize,charthick=2,/dat,color=colors[j]
        if keyword_set(emlsig) then $
           if emlsig.haskey(linlab[i]) then $
              for j=0,maxncomp-1 do $
                 if emlsig[linlab[i],j] ne bad then $
                    cgtext,xran[0]+(xran[1]-xran[0])*0.7d,$
                       yran[0]+(yran[1]-yran[0])*(0.90d - double(j)*0.1d),$
                       '$\sigma\down'+string(j+1,format='(I0)')+'$='+$
                       string(emlsig[linlab[i],j],format='(D0.4)'),$
                       charsize=linelabsize,charthick=2,/dat,color=colors[j]
        if nmasked gt 0 then $
           for j=0,nmasked-1 do $
              cgoplot,[masklam[0,j],masklam[1,j]],[yran[0],yran[0]],thick=8,$
                      color='Cyan'
        ; residual
        pos_res = [xwin[0],ywin[0],$
                   xwin[1],ywin[0]+0.3*dywin]
        if resid eq 'subtotmod' then begin
           ydat = spectot - modtot
           yran = [min(ydat[ind]),max(ydat[ind])]
           ymod = dblarr(n_elements(wave))
        endif else begin
           ydat = specstars
           yran = [min([ydat[ind],ymod[ind]]),max([ydat[ind],ymod[ind]])]
           ymod = modstars
        endelse
        if icol eq fix(icol) then ytit = 'Residual' else ytit = ''
        cgplot,wave,ydat,xran=xran,yran=yran,/noerase,ytit=ytit,$
               axiscol=linecol,col=linecol,/norm,pos=pos_res,/xsty,/ysty,thick=1
        cgoplot,wave,ymod,color='Red',thick=4
        cgoplot,wave,specerr,color='Cyan',thick=2
        ; velocity inset
        if dovelinset then begin
           ; position of inset
           pos_vi = [xwin[0] + dxwin*velinset[1,0],$
              ywin[0] + dywin*velinset[1,1],$
              xwin[0] + dxwin*velinset[1,2],$
              ywin[0] + dywin*velinset[1,3]]
           xraninset=velinset[2,*]
           yraninset=[0d,1d]
           ; color in background
           fillyscale = 0.1d
           cgpolygon,[xwin[0] + dxwin*velinset[1,0],$
              xwin[0] + dxwin*velinset[1,2],$
              xwin[0] + dxwin*velinset[1,2],$
              xwin[0] + dxwin*velinset[1,0],$
              xwin[0] + dxwin*velinset[1,0]],$
              [ywin[0] + dywin*velinset[1,1] - dywin*fillyscale,$
               ywin[0] + dywin*velinset[1,1] - dywin*fillyscale,$
               ywin[0] + dywin*velinset[1,3],$
               ywin[0] + dywin*velinset[1,3],$
               ywin[0] + dywin*velinset[1,1] - dywin*fillyscale],$
              /fill,/normal
           cgplot,[0],[0],xran=xraninset,yran=yraninset,/noerase,$
              pos=pos_vi,xsty=1,ysty=1,xtickint=1000d,xminor=2,$
              chars=!P.charsize*0.75,ytickf='(A1)',xtit='vel (km/s)',xticklen=0.05
           ; plot components
           for j=1,maxncomp do begin
              flux = ifsf_cmplin(instr,velinset[0],j,/velsig)
              zdiff = wave/(instr.linelist[velinset[0]]*(1d + instr.zstar)) - 1d
              vel = c_kms * ((zdiff+1d)^2d - 1d) / ((zdiff+1d)^2d + 1d)
              if j eq 1 then fluxnorm=max(flux)
              cgoplot,vel,flux/fluxnorm,color=colors[j-1],linesty=2,thick=2
           endfor              
           ; label emission line
           ;cgtext,xraninset[0]+(xraninset[1]-xraninset[0])*0.05d,$
           ;   yraninset[0]+(yraninset[1]-yraninset[0])*0.85d,$
           ;   velinset[0],charsize=linelabsize*0.75,charthick=2,/dat,$
           ;   col=linelabcol
        endif

     endif

  endfor

  if tag_exist(pltpar,'micron') then $
    xtit = 'Observed Wavelength ($\mu$m)' $
  else if tag_exist(pltpar,'meter') then $
    xtit = 'Observed Wavelength (m)' $
  else $
     xtit = 'Observed Wavelength (!3' + STRING(197B) + '!X)'
  cgtext,0.5,0.01,xtit,/norm,align=0.5,charsize=1.5,charthick=2
  if keyword_set(title) then $
     cgtext,0.5,0.97,title,/norm,align=0.5,charsize=1.5,charthick=2

  tmpfile = outfile
  if dops then cgps_close $
  else img = cgsnapshot(filename=tmpfile,/jpeg,/nodialog,quality=100)

  !X.tickinterval=defaultXtickint
  !X.minor=defaultXminor
  
end
