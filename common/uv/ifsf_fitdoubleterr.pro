; docformat = 'rst'
;
;+
;
; Monte Carlo the errors in a doublet fit. The input is a doublet fit and the 
; associated data. The error spectrum is varied as a randomized Gaussian and added to the
; best-fit model at each wavelength, and the result is refit. The resulting
; best-fit distributions are parameterized by the median and 68% confidence
; intervals to each side of the median.
;
; The routine can be sped up significantly by using the NSPLIT keyword to
; split the Monte Carlo realizations into NSPLIT separate computations, run in
; parallel. In the workbench on Mac OS X, this cannot exceed the number of
; cores (or the related IDLBrigde objects cannot be destroyed). NSPLIT can
; exceed the number of cores when the routine is run from the command line; this
; will further speed up the code in the case of hyper-threaded cores.
;
; The multi-core functionality requires split_for.pro and its associated
; subroutines, written by Robert da Silva. These are now bundled with DRTOOLS:
; https://github.com/drupke/drtools
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Nparams x 2 array of parameter errors (below and above the median),
;    PDF file of histograms of parameter distributions, and XDR file of
;    Monte Carlo-ed distributions. In output array, order of parameters is
;    HeI emission, doublet absorption, and doublet emission.
;
; :Params:
;    ncomp: in, required, type=dblarr(3)
;      Number of doublet absorption and doublet emission compoennts.
;    wave: in, required, type=dblarr(N)
;      Wavelengths to fit.
;    modflux: in, required, type=dblarr(N)
;      Best-fit model fluxes.
;    err: in, required, type=dblarr(N)
;      Data errors.
;    cont: in, required, type=dblarr(N)
;      Continuum fit; used for computing emission-line fluxes.
;    parinit: in, required, type=structure
;      PARINIT structure for MPFIT; same as used in fit to data.
;    outplot: in, required, type=string
;      Path and filename for PDF file of output parameter distribution
;      histograms. Extension should be '.ps,' and file is then converted to PDF.
;    outfile: in, required, type=string
;      Path and filename for IDL save file (XDR) to hold Monte Carlo results.
;
; :Keywords:
;    dofirstemfit: in, optional, type=byte
;      Do a first fit of doublet emission only before fitting absorption and/or
;      HeI emission.
;    first_modflux: in, optional, type=dblarr(N)
;      Best-fit model fluxes for first fit.
;    first_parinit: in, optional, type=structure
;      PARINIT for first fit.
;    fitfcn: in, optional, type=string, default='ifsf_doubletfcn'
;      Function to fit; input to MPFIT.
;    niter: in, optional, type=integer, default=1000
;      Number of Monte Carlo realizations.
;    noplot: in, optional, type=byte
;      Do not produce plots.
;    nsplit: in, optional, type=integer, default=1
;      Number of independent child processes into which the Monte Carlo
;      computation is split.
;    doubletemfix: in, optional, type=bytarr(Ncomp_doubletem,3)
;      Array of fix/free flags for doublet emission parameters; same as that used
;      in the fit to the data.
;    plotonly: in, optional, type=boolean
;      Do not re-run MC simulations; grab previous outputs from file.
;    quiet: in, optional, type=byte
;      Suppress error and progress messages. Preently controls only SPLIT_FOR
;      messages; MPFIT messages are suppressed b/c of number of times it is
;      called ...
;    weqerr: out, optional, type=dblarr(2,2)
;      Emission and absorption line equivalent width errors (low and high).
;      First element of first dimension is absorption, second is emission.
;    doubletemfluxerr: out, optional, type=dblarr(2)
;      Low and high emission line flux errors.
;    specres: in, required, type=double, def=0.64d
;      Estimated spectral resolution in wavelength units (sigma).
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
;      2021aug26, DSNR, adapted from IFSF_FITNADERR
;
; :Copyright:
;    Copyright (C) 2021 David S. N. Rupke
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
function ifsf_fitdoubleterr,doubletname,ncomp,wave,vels,modflux,err,cont,parinit,$
   outplot,outfile,dofirstemfit=dofirstemfit,first_modflux=first_modflux,$
   first_parinit=first_parinit,fitfcn=fitfcn,niter=niter,nsplit=nsplit,$
   doubletabsfix=doubletabsfix,doubletemfix=doubletemfix,quiet=quiet,$
   noplot=noplot,weqerr=weqerr,doubletemfluxerr=doubletemfluxerr,$
   plotonly=plotonly,specres=specres,vwtavgerr=vwtavgerr,vwtrmserr=vwtrmserr

   bad = 1d99
   plotquantum = 2.5 ; in inches

   if ~ keyword_set(fitfcn) then fitfcn = 'ifsf_doubletfcn'
   if ~ keyword_set(niter) then niter = 1000
   if ~ keyword_set(nsplit) then nsplit = 1
   if ~ keyword_set(quiet) then quiet=0
   if ~ keyword_set(noplot) then noplot=0
   if ncomp[0] gt 0 AND ~ keyword_set(doubletabsfix) then $
      doubletabsfix=bytarr(ncomp[0],4)
   if ncomp[1] gt 0 AND ~ keyword_set(doubletemfix) then $
      doubletemfix=bytarr(ncomp[1],4)
   nwave = n_elements(wave)
   if ~ keyword_set(specres) then specres = 0b
   argsfitdoublet = {doubletname: doubletname, specres: specres}

   ;  Outputs
   errors = dblarr(ncomp[0]*4+ncomp[1]*4,2)
   weqerr=dblarr(2,2)
   vwtavgerr=dblarr(2)
   vwtrmserr=dblarr(2)
   doubletemfluxerr=dblarr(2)

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; Monte Carlo
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if ~ keyword_set(plotonly) then begin

      if keyword_set(dofirstemfit) then nfits=2 else nfits=1
      for j=0,nfits-1 do begin

         if keyword_set(dofirstemfit) then begin
            if j eq 0 then begin
               parinituse = first_parinit
               modfluxinit = first_modflux
               ncompuse = [0,ncomp[2]]
            endif else begin
               parinituse = parinit
               modfluxinit = modflux
               ncompuse = [ncomp[0],0]
            endelse
         endif else begin
            parinituse = parinit
            modfluxinit = modflux
            ncompuse = ncomp
         endelse

         rans = randomn(seed,nwave,niter,/double)

         ;     This block processes on one core...
         if nsplit eq 1 then begin

            if ncompuse[0] gt 0 then abspar = dblarr(niter,4,ncompuse[0])+bad $
            else abspar=0
            if ncompuse[1] gt 0 then empar = dblarr(niter,4,ncompuse[1])+bad $
            else empar=0
            weq=dblarr(niter,2)
            vwtavg=dblarr(niter)
            vwtrms=dblarr(niter)
            doubletemflux=dblarr(niter)
            for k=0,niter-1 do begin

               modflux_use = modfluxinit + rans[*,k]*err

               param = Mpfitfun(fitfcn,wave,modflux_use,err,$
                  parinfo=parinituse,perror=perror,maxiter=100,$
                  bestnorm=chisq,covar=covar,yfit=specfit,dof=dof,$
                  nfev=nfev,niter=niter_mpfit,status=status,quiet=quiet,$
                  npegged=npegged,ftol=1D-6,errmsg=errmsg,$
                  functargs=argsfitdoublet)
               if status eq 5 then message,'Max. iterations reached.',/cont
               if status eq 0 OR status eq -16 then $
                  message,'Error in MPFIT.'

               weq_i=1
               doubletemflux_i=1
               vwtabs_i=1
               dumy = call_function(fitfcn,wave,param,weq=weq_i,$
                  doubletemflux=doubletemflux_i,cont=cont,_extra=argsfitdoublet,$
                  vwtabs=vwtabs_i,vels=vels)

               if ncompuse[0] gt 0 then begin
                  ilo = 2
                  abspar[k,*,*] = reform(param[ilo:ilo+4*ncompuse[0]-1],4,ncompuse[0])
               endif
               if ncompuse[1] gt 0 then begin
                  ilo = 2+4*ncompuse[0]
                  empar[k,*,*] = reform(param[ilo:ilo+4*ncompuse[1]-1],4,ncompuse[1])
               endif
               weq[k,*] = [weq_i.abs[0],weq_i.em[0]]
               vwtavg[k] = vwtabs_i[0]
               vwtrms[k] = vwtabs_i[1]
               doubletemflux[k] = doubletemflux_i[0]

               if k mod 100 eq 0 then $
                  message,'Iteration '+string(k,format='(I0)'),/cont
            endfor

            ;     This block splits to nsplit children...speedup is a factor of 5.4 for 7
            ;     children on my 8-core Xeon Mac Pro (measured with 1000 iterations,
            ;     spaxel (14,14) of F05189-2524 cube).
         endif else begin

            split_for,0,niter-1,commands=[$
               'resolve_routine, '+string(39B)+'mpfit'+string(39B)+', /EITHER, /NO_RECOMPILE',$
               'modflux_use = modfluxinit + rans[*,i]*err',$
               'param = Mpfitfun(fitfcn,wave,modflux_use,err,'+$
               '                 parinfo=parinituse,perror=perror,maxiter=100,'+$
               '                 bestnorm=chisq,covar=covar,yfit=specfit,dof=dof,'+$
               '                 nfev=nfev,niter=niter_mpfit,status=status,quiet=quiet,'+$
               '                 functargs=argsfitdoublet,'+$
               '                 npegged=npegged,ftol=1D-6,errmsg=errmsg)',$
               'weq_i=1',$
               'vwtabs_i=1',$
               'doubletemflux_i=1',$
               'dumy = call_function(fitfcn,wave,param,weq=weq_i,doubletemflux=doubletemflux_i,'+$
               '                    cont=cont,_extra=argsfitdoublet,vwtabs=vwtabs_i,vels=vels)',$
               'if ncompuse[0] gt 0 then begin',$
               '   ilo = 2',$
               '   if n_elements(abspar) eq 0 then'+$
               '      abspar = reform(param[ilo:ilo+4*ncompuse[0]-1],1,4,ncompuse[0])'+$
               '   else abspar = [abspar,reform(param[ilo:ilo+4*ncompuse[0]-1],1,4,ncompuse[0])]',$
               'endif else abspar=0',$
               'if ncompuse[1] gt 0 then begin',$
               '   ilo = 2+4*ncompuse[0]',$
               '   if n_elements(empar) eq 0 then'+$
               '      empar = reform(param[ilo:ilo+4*ncompuse[1]-1],1,4,ncompuse[1])'+$
               '   else empar = [empar,reform(param[ilo:ilo+4*ncompuse[1]-1],1,4,ncompuse[1])]',$
               'endif else empar=0',$
               'if n_elements(weq) eq 0 then begin',$
               '   weq=reform([weq_i.abs[0],weq_i.em[0]],1,2)',$
               '   vwtavg=vwtabs_i[0]',$
               '   vwtrms=vwtabs_i[1]',$
               'endif else begin',$
               '   weq=[weq,reform([weq_i.abs[0],weq_i.em[0]],1,2)]',$
               '   vwtavg=[vwtavg,vwtabs_i[0]]',$
               '   vwtrms=[vwtrms,vwtabs_i[1]]',$
               'endelse',$
               'if n_elements(doubletemflux) eq 0 then'+$
               '   doubletemflux=doubletemflux_i[0]'+$
               'else doubletemflux=[doubletemflux,doubletemflux_i[0]]'],$
               varnames=['modfluxinit','rans','err','fitfcn','wave','vels','ncompuse',$
               'cont','quiet'],$
               struct2pass1=parinituse,struct2pass2=argslinefit,$
               struct2pass3=argsfitdoublet,$
               outvar=['abspar','empar','weq','doubletemflux','vwtavg','vwtrms'],$
               nsplit=nsplit,silent=quiet,quietchild=quiet
            abspar=abspar0
            empar=empar0
            weq=weq0
            vwtavg=vwtavg0
            vwtrms=vwtrms0
            doubletemflux=doubletemflux0
            for k=1,nsplit-1 do begin
               abspar=[abspar,scope_varfetch('abspar'+string(k,format='(I0)'))]
               empar=[empar,scope_varfetch('empar'+string(k,format='(I0)'))]
               weq=[weq,scope_varfetch('weq'+string(k,format='(I0)'))]
               vwtavg=[vwtavg,scope_varfetch('vwtavg'+string(k,format='(I0)'))]
               vwtrms=[vwtrms,scope_varfetch('vwtrms'+string(k,format='(I0)'))]
               doubletemflux=[doubletemflux,scope_varfetch('doubletemflux'+string(k,format='(I0)'))]
            endfor

         endelse

         if keyword_set(dofirstemfit) then begin
            if j eq 0 then begin
               first_empar = empar
               first_weq = weq
               first_vwtavg = vwtavg
               first_vwtrms = vwtrms
               first_doubletemflux = doubletemflux
            endif
            if j eq 1 then begin
               empar = first_empar
               weq[*,1] = first_weq[*,1]
               vwtavg = first_vwtavg
               vwtrms = first_vwtrms
               doubletemflux = first_doubletemflux
            endif
         endif

      endfor

      ;  Save MC results to a file.
      doubletmc = {doubletabs: abspar, doubletem: empar, doubletweq: weq, $
         doubletemflux: doubletemflux, doubletvwtavg: vwtavg, $
         doubletvwtrms: vwtrms}
      save,doubletmc,file=outfile

   endif else begin

      restore,file=outfile
      abspar=doubletmc.doubletabs
      empar=doubletmc.doubletem
      weq=doubletmc.doubletweq
      vwtavg=doubletmc.doubletvwtavg
      vwtrms=doubletmc.doubletvwtrms
      doubletemflux=doubletmc.doubletemflux

   endelse

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; Compute errors from distributions and plot results
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if noplot then set_plot,'x' $
   else cgps_open,outplot,charsize=1,/inches,xs=7.5,ys=7.5,/qui,/nomatch
   pos = cglayout([2,2],ixmar=[3,3],iymar=[4,4],oxmar=[0,0],oymar=[4,1],$
      xgap=0,ygap=4,aspect=1d)
   ;
   ;  doublet abs
   ;
   for i=0,ncomp[0]-1 do begin
      ;     covering factor x optical depth
      if (~ doubletabsfix[i,0] OR ~ doubletabsfix[i,1]) then begin
         dat = abspar[*,0,i]*abspar[*,1,i]
         errtmp = $
            ifsf_pltmcdist(dat,position=pos[*,0],xtitle='C$\down$f x $\tau$',$
               mininput=0d,min_value=0d,xrange=[0,max(dat)])
         errors[i*4,*]=errtmp/2d
      endif else $
         cghistoplot,abspar[*,0,i]*abspar[*,1,i],pos=pos[*,0],missing=bad*bad,$
            ytit='',xtit='C$\down$f x $\tau$',xticks=2,nbins=1,$
            bins=0.1d,mininput=0d,min_value=0d
      ;     optical depth
      isame = uniq(float(abspar[*,1,i]))
      if ~ doubletabsfix[i,1] AND n_elements(isame) ne 1  then begin
         dat_tmp = abspar[*,1,i]
         igd_tmp = where(dat_tmp ne 0 AND dat_tmp ne bad)
         dat = alog10(dat_tmp[igd_tmp])
         errors[i*4+1,*]=$
            ifsf_pltmcdist(dat,position=pos[*,1],xtitle='log($\tau$)',$
               /noerase)
      endif else $
         cghistoplot,dat,pos=pos[*,1],/noerase,missing=bad,$
            ytit='',xtit='log($\tau$)',xticks=2,nbins=1,bins=0.1d,$
            mininput=0d,min_value=0d
      ;     wavelength
      if ~ doubletabsfix[i,2] then $
         errors[i*4+2,*] = $
            ifsf_pltmcdist(abspar[*,2,i],position=pos[*,2],$
               xtitle='$\lambda$ ($\Angstrom$)',/noerase) $
      else cghistoplot,dat,pos=pos[*,2],missing=bad,/noerase,ytit='',$
         xtit='$\lambda$ ($\Angstrom$)',xticks=2,nbins=1,bins=1d
      ;     sigma
      if ~ doubletabsfix[i,3] then $
         errors[i*4+3,*] = $
            ifsf_pltmcdist(abspar[*,3,i],position=pos[*,3],$
               xtitle='$\sigma$ (km/s)',/noerase,$
               mininput=0d,min_value=0d,xran=[0,max(abspar[*,3,i])]) $
      else cghistoplot,dat,pos=pos[*,3],missing=bad,/noerase,ytit='',$
         xtit='$\sigma$ (km/s)',xticks=2,nbins=1,bins=10d,$
         mininput=0d,min_value=0d
      ;
      cgtext,0.5,0.99,'doublet abs., c'+string(i+1,format='(I0)'),/norm,align=0.5
   endfor
   ;
   ;  doublet em
   ;
   for i=0,ncomp[1]-1 do begin
      ;     wavelength
      if ~ doubletemfix[i,0] then $
         errors[ncomp[0]*4+i*4,*] = $
            ifsf_pltmcdist(empar[*,0,i],position=pos[*,0],$
               xtitle='$\lambda$ ($\Angstrom$)') $
      else cghistoplot,empar[*,0,i],pos=pos[*,0],missing=bad,ytit='',$
         xtit='$\lambda$ ($\Angstrom$)',xticks=2,nbins=1,binsize=1d
      ;     sigma
      isame = uniq(float(empar[*,1,i])) ; float deals with machine precision issues
      if ~ doubletemfix[i,1] AND n_elements(isame) ne 1 then $
         errors[ncomp[0]*4+i*4+1,*] = $
            ifsf_pltmcdist(empar[*,1,i],position=pos[*,1],$
               xtitle='$\sigma$ (km/s)',/noerase,$
               mininput=0d,min_value=0d,xran=[0,max(empar[*,1,i])]) $
      else cghistoplot,empar[*,1,i],pos=pos[*,1],missing=bad,/noerase,ytit='',$
         xtit='$\sigma$ (km/s)',xticks=2,nbins=1,binsize=10d,$
         mininput=0d,min_value=0d
      ;     peak flux
      inonzero = where(empar[*,2,i] ne 0d,ctnonzero)
      if ~ doubletemfix[i,2] AND ctnonzero gt 0 then $
         errors[ncomp[0]*4+i*4+2,*] = $
            ifsf_pltmcdist(empar[*,2,i],position=pos[*,2],$
               xtitle='F$\down\\lambda$ (peak)',/noerase,$
               mininput=0d,min_value=0d,xran=[0,max(empar[*,2,i])]) $
      else cghistoplot,empar[*,2,i],pos=pos[*,2],missing=bad,/noerase,ytit='',$
         xtit='F$\down\\lambda$ (peak)',xticks=3,nbins=1,binsize=0.1d,$
         mininput=0d,min_value=0d
      ;     flux ratio
      isame = uniq(float(empar[*,3,i])) ; float deals with machine precision issues
      if ~ doubletemfix[i,3] AND n_elements(isame) ne 1 then $
         errors[ncomp[0]*4+i*4+3,*] = $
            ifsf_pltmcdist(empar[*,3,i],position=pos[*,3],$
               xtitle='F$\down2$/F$\down1$',/noerase) $
      else cghistoplot,empar[*,3,i],pos=pos[*,3],missing=bad,/noerase,ytit='',$
         xtit='F$\down2$/F$\down1$',xticks=2,nbins=1,binsize=0.1d
      cgtext,0.5,0.99,'doublet em., c'+string(i+1,format='(I0)'),/norm,align=0.5
   endfor
   ;
   ; doublet equivalent widths and fluxes
   ;
   if ncomp[0] gt 0 OR ncomp[1] gt 0 then begin
      ;     absorption
      if ncomp[0] gt 0 then begin
         em_noerase=1
         ;        equivalent width
         weqerr[0,*] = $
            ifsf_pltmcdist(weq[*,0],position=pos[*,0],$
               xtitle='W$\downeq$ (abs, $\Angstrom$)',$
               mininput=0d,min_value=0d,xran=[0,max(weq[*,0])])
      endif else em_noerase=0
      ;     emission
      if ncomp[1] gt 0 then begin
         ;        equivalent width
         dat = weq[*,1]
         inonzero = where(dat ne 0d,ctnonzero)
         if ctnonzero gt 0 then $
            weqerr[1,*] = $
               ifsf_pltmcdist(weq[*,1],position=pos[*,2],$
                  xtitle='W$\downeq$ (em, $\Angstrom$)',$
                  maxinput=0d,min_value=0d,xran=[min(weq[*,1]),0],$
                  noerase=em_noerase) $
         else cghistoplot,dat,pos=pos[*,2],missing=bad,noerase=em_noerase,ytit='',$
            xtit='W$\downeq$ (em, $\Angstrom$)',xticks=2,nbins=1,binsize=0.1d,$
            maxinput=0d,max_value=0d
         ;        flux
         dat = doubletemflux
         inonzero = where(dat ne 0d,ctnonzero)
         if ctnonzero gt 0 then $
            doubletemfluxerr = $
               ifsf_pltmcdist(doubletemflux,position=pos[*,3],xtitle='Flux (em)',$
                  mininput=0d,min_value=0d,xran=[0,max(doubletemflux)],$
                  /noerase) $
         else cghistoplot,dat,pos=pos[*,3],missing=bad,/noerase,ytit='',$
            xtit='Flux (em)',xticks=2,nbins=1,binsize=0.1d,$
            mininput=0d,min_value=0d
      endif
      cgtext,0.5,0.99,'doublet equivalent widths'+$
         string(i+1,format='(I0)'),/norm,align=0.5

   endif
   ;
   ; absorption doublet velocities
   ;
   if ncomp[0] gt 0 then begin
      ; vwtavg
      vwtavgerr = $
         ifsf_pltmcdist(vwtavg,position=pos[*,0],$
            xtitle='Average depth-weighted velocity (km/s)',$
            xran=[min(vwtavg),max(vwtavg)])
      vwtrmserr = $
         ifsf_pltmcdist(vwtrms,position=pos[*,1],$
            xtitle='RMS depth-weighted velocity (km/s)',$
            xran=[min(vwtrms),max(vwtrms)],/noerase)

   endif

   if ~ noplot then cgps_close,/pdf,/delete_ps

   return,errors

end
