; docformat = 'rst'
;
;+
;
; Monte Carlo the errors in a UV multiplet fit. The input is a multiplet fit and the 
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
;    Monte Carlo-ed distributions.
;
; :Params:
;    ncomp: in, required, type=double
;      Number of multiplet absorption components
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
;    fitfcn: in, optional, type=string, default='ifsf_multipletfcn'
;      Function to fit; input to MPFIT.
;    niter: in, optional, type=integer, default=1000
;      Number of Monte Carlo realizations.
;    noplot: in, optional, type=byte
;      Do not produce plots.
;    nsplit: in, optional, type=integer, default=1
;      Number of independent child processes into which the Monte Carlo
;      computation is split.
;    plotonly: in, optional, type=boolean
;      Do not re-run MC simulations; grab previous outputs from file.
;    quiet: in, optional, type=byte
;      Suppress error and progress messages. Preently controls only SPLIT_FOR
;      messages; MPFIT messages are suppressed b/c of number of times it is
;      called ...
;    weqerr: out, optional, type=dblarr(2)
;      Absorption line equivalent width errors (low and high).
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
function ifsf_fitmultipleterr,argsfit,ncomp,wave,vels,modflux,err,cont,parinit,$
   outplot,outfile,fitfcn=fitfcn,niter=niter,nsplit=nsplit,$
   absfix=absfix,quiet=quiet,noplot=noplot,weqerr=weqerr,plotonly=plotonly,$
   vwtavgerr=vwtavgerr,vwtrmserr=vwtrmserr
   
   bad = 1d99
   plotquantum = 2.5 ; in inches

   if ~ keyword_set(fitfcn) then fitfcn = 'ifsf_multipletfcn'
   if ~ keyword_set(niter) then niter = 1000
   if ~ keyword_set(nsplit) then nsplit = 1
   if ~ keyword_set(quiet) then quiet=0
   if ~ keyword_set(noplot) then noplot=0
   if ncomp gt 0 AND ~ keyword_set(absfix) then $
      absfix=bytarr(ncomp,4)
   nwave = n_elements(wave)

   ;  Outputs
   errors = dblarr(ncomp*4,2)
   weqerr=dblarr(2)
   vwtavgerr=dblarr(2)
   vwtrmserr=dblarr(2)

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; Monte Carlo
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if ~ keyword_set(plotonly) then begin

      parinituse = parinit
      modfluxinit = modflux
      ncompuse = ncomp

      ; Use a Gaussian random seed with sigma = uncertainty. RANDOMN produces
      ; a draw from a Gaussian with sigma=1, so when we multiply by our error
      ; this scales the x-value of the distribution by using the 
      ; x=1 value from the RANDOMN draw. (It doesn't scale the y-value of the
      ; Gaussian distribution in fluxes....)
      rans = randomn(seed,nwave,niter,/double)

      ;     This block processes on one core...
      if nsplit eq 1 then begin

         if ncompuse gt 0 then abspar = dblarr(niter,4,ncompuse)+bad $
         else abspar=0
         weq=dblarr(niter)
         vwtavg=dblarr(niter)
         vwtrms=dblarr(niter)
         for k=0,niter-1 do begin
            modflux_use = modfluxinit + rans[*,k]*err
            param = Mpfitfun(fitfcn,wave,modflux_use,err,$
               parinfo=parinituse,perror=perror,maxiter=100,$
               bestnorm=chisq,covar=covar,yfit=specfit,dof=dof,$
               nfev=nfev,niter=niter_mpfit,status=status,quiet=quiet,$
               npegged=npegged,ftol=1D-6,errmsg=errmsg,$
               functargs=argsfit)
            if status eq 5 then message,'Max. iterations reached.',/cont
            if status eq 0 OR status eq -16 OR status eq -18 then $
               message,'Error in MPFIT.'
            weq_i=1
            vwtabs_i=1
            dumy = call_function(fitfcn,wave,param,weq=weq_i,$
               cont=cont,_extra=argsfit,vwtabs=vwtabs_i,vels=vels,$
               specres=argsfit.specres,upsample=argsfit.upsample)
            if ncompuse gt 0 then begin
               ilo = 1
               abspar[k,*,*] = reform(param[ilo:ilo+4*ncompuse-1],4,ncompuse)
            endif
            if k mod 100 eq 0 then $
               message,'Iteration '+string(k,format='(I0)'),/cont
            weq[k]=weq_i.abs[0]
            vwtavg[k]=vwtabs_i[0]
            vwtrms[k]=vwtabs_i[1]
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
            '                 functargs=argsfit,'+$
            '                 npegged=npegged,ftol=1D-6,errmsg=errmsg)',$
            'weq_i=1',$
            'vwtabs_i=1',$
            'dumy = call_function(fitfcn,wave,param,weq=weq_i,'+$
            '                    cont=cont,_extra=argsfit,vwtabs=vwtabs_i,vels=vels,'+$
            '                    specres=argsfit.specres,upsample=argsfit.upsample)',$
            'if ncompuse gt 0 then begin',$
            '   ilo = 1',$
            '   if n_elements(abspar) eq 0 then'+$
            '      abspar = reform(param[ilo:ilo+4*ncompuse-1],1,4,ncompuse)'+$
            '   else abspar = [abspar,reform(param[ilo:ilo+4*ncompuse[0]-1],1,4,ncompuse[0])]',$
            'endif else abspar=0',$
            'if n_elements(weq) eq 0 then begin',$
            '   weq=weq_i.abs[0]',$
            '   vwtavg=vwtabs_i[0]',$
            '   vwtrms=vwtabs_i[1]',$
            'endif else begin',$
            '   weq=[weq,weq_i.abs[0]]',$
            '   vwtavg=[vwtavg,vwtabs_i[0]]',$
            '   vwtrms=[vwtrms,vwtabs_i[1]]',$
            'endelse'],$
            varnames=['modfluxinit','rans','err','fitfcn','wave','vels','ncompuse',$
               'cont','quiet'],$
            struct2pass1=parinituse,struct2pass2=argslinefit,$
            struct2pass3=argsfit,$
            outvar=['abspar','weq','vwtavg','vwtrms'],$
            nsplit=nsplit,silent=quiet,quietchild=quiet
            abspar=abspar0
            weq=weq0
            vwtavg=vwtavg0
            vwtrms=vwtrms0
         for k=1,nsplit-1 do begin
            abspar=[abspar,scope_varfetch('abspar'+string(k,format='(I0)'))]
            weq=[weq,scope_varfetch('weq'+string(k,format='(I0)'))]
            vwtavg=[vwtavg,scope_varfetch('vwtavg'+string(k,format='(I0)'))]
            vwtrms=[vwtrms,scope_varfetch('vwtrms'+string(k,format='(I0)'))]
         endfor

      endelse

      ;  Save MC results to a file.
      multipletmc = {abs: abspar, weq: weq, vwtavg: vwtavg, vwtrms: vwtrms}
      save,multipletmc,file=outfile

   endif else begin

      restore,file=outfile
      abspar=multipletmc.abs
      weq=multipletmc.weq
      vwtavg=multipletmc.vwtavg
      vwtrms=multipletmc.vwtrms

   endelse

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; Compute errors from distributions and plot results
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if noplot then set_plot,'x' $
   else cgps_open,outplot,charsize=1,/inches,xs=7.5,ys=7.5,/qui,/nomatch
   pos = cglayout([2,2],ixmar=[3,3],iymar=[4,4],oxmar=[0,0],oymar=[4,1],$
      xgap=0,ygap=4,aspect=1d)
   ;
   ;  multiplet abs
   ;
   for i=0,ncomp-1 do begin
      ;     covering factor x optical depth
      if (~ absfix[i,0] OR ~ absfix[i,1]) then begin
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
      if ~ absfix[i,1] AND n_elements(isame) ne 1  then begin
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
      if ~ absfix[i,2] then $
         errors[i*4+2,*] = $
            ifsf_pltmcdist(abspar[*,2,i],position=pos[*,2],$
               xtitle='$\lambda$ ($\Angstrom$)',/noerase) $
      else cghistoplot,dat,pos=pos[*,2],missing=bad,/noerase,ytit='',$
         xtit='$\lambda$ ($\Angstrom$)',xticks=2,nbins=1,bins=1d
      ;     sigma
      if ~ absfix[i,3] then $
         errors[i*4+3,*] = $
            ifsf_pltmcdist(abspar[*,3,i],position=pos[*,3],$
               xtitle='$\sigma$ (km/s)',/noerase,$
               mininput=0d,min_value=0d,xran=[0,max(abspar[*,3,i])]) $
      else cghistoplot,dat,pos=pos[*,3],missing=bad,/noerase,ytit='',$
         xtit='$\sigma$ (km/s)',xticks=2,nbins=1,bins=10d,$
         mininput=0d,min_value=0d
      ;
      cgtext,0.5,0.99,'multiplet abs., c'+string(i+1,format='(I0)'),/norm,align=0.5
   endfor
   ;
   ; multiplet equivalent widths and velocities
   ;
   if ncomp gt 0 then begin
   ;        equivalent width
      weqerr = $
         ifsf_pltmcdist(weq,position=pos[*,0],$
            xtitle='W$\downeq$ (abs, $\Angstrom$)',$
            mininput=0d,min_value=0d,xran=[0,max(weq[*,0])])
   ; vwtavg
      vwtavgerr = $
         ifsf_pltmcdist(vwtavg,position=pos[*,1],$
            xtitle='Average depth-weighted velocity (km/s)',$
            xran=[min(vwtavg),max(vwtavg)],/noerase)
      vwtrmserr = $
         ifsf_pltmcdist(vwtrms,position=pos[*,2],$
            xtitle='RMS depth-weighted velocity (km/s)',$
            xran=[min(vwtrms),max(vwtrms)],/noerase)
      cgtext,0.5,0.99,'multiplet equivalent widths + vels',/norm,align=0.5
   endif

   if ~ noplot then cgps_close,/pdf,/delete_ps

   return,errors

end
