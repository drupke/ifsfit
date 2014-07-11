; docformat = 'rst'
;
;+
;
; Monte Carlo the errors in a NaD fit. The input is a NaD fit and the associated
; data. The error spectrum is varied as a randomized Gaussian and added to the 
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
; subroutines, written by Robert da Silva and available at
;   http://slugidl.pbworks.com/w/page/29199259/Child%20Processes
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Nparams x 2 array of parameter errors (below and above the median), 
;    PDF file of histograms of parameter distributions, and XDR file of 
;    Monte Carlo-ed distributions. In output array, order of parameters is
;    HeI emission, NaD absorption, and NaD emission.
;    
; :Params:
;    ncomp: in, required, type=dblarr(3)
;      Number of HeI emission, NaD absorption, and NaD emission compoennts.
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
;      Do a first fit of NaD emission only before fitting absorption and/or
;      HeI emission.
;    first_modflux: in, optional, type=dblarr(N)
;      Best-fit model fluxes for first fit.
;    first_parinit: in, optional, type=structure
;      PARINIT for first fit.
;    fitfcn: in, optional, type=string, default='ifsf_nadfcn'
;      Function to fit; input to MPFIT.
;    niter: in, optional, type=integer, default=1000
;      Number of Monte Carlo realizations.
;    nsplit: in, optional, type=integer, default=1
;      Number of independent child processes into which the Monte Carlo 
;      computation is split.
;    heifix: in, optional, type=bytarr(Ncomp_HeI,3)
;      Array of fix/free flags for HeI parameters; same as that used in the fit
;      to the data.
;    nademfix: in, optional, type=bytarr(Ncomp_NaDem,3)
;      Array of fix/free flags for NaD emission parameters; same as that used 
;      in the fit to the data.
;    quiet: in, optional, type=byte
;      Suppress error and progress messages. Preently controls only SPLIT_FOR
;      messages; MPFIT messages are suppressed b/c of number of times it is
;      called ...
;    weqerr: out, optional, type=dblarr(2,2)
;      Emission and absorption line equivalent width errors (low and high).
;      First element of first dimension is absorption, second is emission.
;    nademfluxerr: out, optional, type=dblarr(2)
;      Low and high emission line flux errors.
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
;      2014jun19, DSNR, created
;      2014jul03, DSNR, added multi-core functionality
;      2014jul07, DSNR, documented
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
function ifsf_fitnaderr,ncomp,wave,modflux,err,cont,parinit,outplot,outfile,$
                        dofirstemfit=dofirstemfit,$
                        first_modflux=first_modflux,$
                        first_parinit=first_parinit,$
                        fitfcn=fitfcn,niter=niter,nsplit=nsplit,$
                        heifix=heifix,nademfix=nademfix,quiet=quiet,$
                        weqerr=weqerr,nademfluxerr=nademfluxerr

   bad = 1d99
   plotquantum = 2.5 ; in inches
   nadlinelist = ifsf_linelist(['NaD1','NaD2','HeI5876'])

   if ~ keyword_set(fitfcn) then fitfcn = 'ifsf_nadfcn'
   if ~ keyword_set(niter) then niter = 1000
   if ~ keyword_set(nsplit) then nsplit = 1
   if ~ keyword_set(quiet) then quiet=0
   if ncomp[0] gt 0 AND ~ keyword_set(heifix) then heifix=bytarr(ncomp[0],3)  
   if ncomp[2] gt 0 AND ~ keyword_set(nademfix) then nademfix=bytarr(ncomp[2],4)
   nwave = n_elements(wave)

;  Outputs
   errors = dblarr(ncomp[0]*3+ncomp[1]*4+ncomp[2]*4,2)
   weqerr=dblarr(2,2)
   nademfluxerr=dblarr(2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Monte Carlo
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if keyword_set(dofirstemfit) then nfits=2 else nfits=1
   for j=0,nfits-1 do begin

      if keyword_set(dofirstemfit) then begin
         if j eq 0 then begin
            parinituse = first_parinit
            modfluxinit = first_modflux
            ncompuse = [0,0,ncomp[2]]
          endif else begin
            parinituse = parinit
            modfluxinit = modflux
            ncompuse = [ncomp[0],ncomp[1],0]
          endelse
      endif else begin
         parinituse = parinit
         modfluxinit = modflux
         ncompuse = ncomp
      endelse

      rans = randomn(seed,nwave,niter,/double)

;     This block processes on one core...
      if nsplit eq 1 then begin
         
         if ncompuse[0] gt 0 then heipar = dblarr(niter,3,ncompuse[0])+bad $
            else heipar=0
         if ncompuse[1] gt 0 then abspar = dblarr(niter,4,ncompuse[1])+bad $
            else abspar=0
         if ncompuse[2] gt 0 then empar = dblarr(niter,4,ncompuse[2])+bad $
            else empar=0
         weq=dblarr(niter,2)
         nademflux=dblarr(niter)
         for k=0,niter-1 do begin

            modflux_use = modfluxinit + rans[*,k]*err

            param = Mpfitfun(fitfcn,wave,modflux_use,err,$
                             parinfo=parinituse,perror=perror,maxiter=100,$
                             bestnorm=chisq,covar=covar,yfit=specfit,dof=dof,$
                             nfev=nfev,niter=niter_mpfit,status=status,/quiet,$
                             npegged=npegged,ftol=1D-6,errmsg=errmsg)
            if status eq 5 then print,'IFSF_FITNADERR: Max. iterations reached.'
            if status eq 0 OR status eq -16 then begin
               print,'IFSF_FITNADERR: Error in MPFIT. Aborting.'
               goto,finish
            endif

            weq_i=1
            nademflux_i=1
            dumy = ifsf_nadfcn(wave,param,weq=weq_i,nademflux=nademflux_i,$
                               cont=cont)

            if ncompuse[0] gt 0 then begin
               ilo = 3
               heipar[k,*,*] = reform(param[ilo:ilo+3*ncompuse[0]-1],3,ncompuse[0])
            endif
            if ncompuse[1] gt 0 then begin
               ilo = 3+3*ncompuse[0]
               abspar[k,*,*] = reform(param[ilo:ilo+4*ncompuse[1]-1],4,ncompuse[1])
            endif
            if ncompuse[2] gt 0 then begin
               ilo = 3+3*ncompuse[0]+4*ncompuse[1]
               empar[k,*,*] = reform(param[ilo:ilo+4*ncompuse[2]-1],4,ncompuse[2])
            endif
            weq[k,*] = [weq_i.abs[0],weq_i.em[0]]
            nademflux[k] = nademflux_i[0]

            if k mod 100 eq 0 then print,'IFSF_FITNADERR: Iteration ',$
                                      string(k,format='(I0)')
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
            '                 nfev=nfev,niter=niter_mpfit,status=status,/quiet,'+$
            '                 npegged=npegged,ftol=1D-6,errmsg=errmsg)',$
            'weq_i=1',$
            'nademflux_i=1',$
            'dumy = ifsf_nadfcn(wave,param,weq=weq_i,nademflux=nademflux_i,'+$
            '                   cont=cont)',$
            'if ncompuse[0] gt 0 then begin',$
            '   ilo = 3',$
            '   if n_elements(heipar) eq 0 then'+$
            '      heipar = reform(param[ilo:ilo+3*ncompuse[0]-1],1,3,ncompuse[0])'+$
            '   else heipar = [heipar,reform(param[ilo:ilo+3*ncompuse[0]-1],1,3,ncompuse[0])]',$
            'endif else heipar=0',$
            'if ncompuse[1] gt 0 then begin',$
            '   ilo = 3+3*ncompuse[0]',$
            '   if n_elements(abspar) eq 0 then'+$
            '      abspar = reform(param[ilo:ilo+4*ncompuse[1]-1],1,4,ncompuse[1])'+$
            '   else abspar = [abspar,reform(param[ilo:ilo+4*ncompuse[1]-1],1,4,ncompuse[1])]',$
            'endif else abspar=0',$
            'if ncompuse[2] gt 0 then begin',$
            '  ilo = 3+3*ncompuse[0]+4*ncompuse[1]',$
            '  if n_elements(empar) eq 0 then'+$
            '     empar = reform(param[ilo:ilo+4*ncompuse[2]-1],1,4,ncompuse[2])'+$
            '  else empar = [empar,reform(param[ilo:ilo+4*ncompuse[2]-1],1,4,ncompuse[2])]',$
            'endif else empar=0',$
            'if n_elements(weq) eq 0 then'+$
            '   weq=reform([weq_i.abs[0],weq_i.em[0]],1,2)'+$
            'else weq=[weq,reform([weq_i.abs[0],weq_i.em[0]],1,2)]',$
            'if n_elements(nademflux) eq 0 then'+$
            '   nademflux=nademflux_i[0]'+$
            'else nademflux=[nademflux,nademflux_i[0]]'],$
            varnames=['modfluxinit','rans','err','fitfcn','wave','ncompuse',$
                      'cont'],$
            struct2pass1=parinituse,struct2pass2=argslinefit,$
            outvar=['heipar','abspar','empar','weq','nademflux'],$
            nsplit=nsplit,silent=quiet,$
            quietchild=quiet
         heipar=heipar0
         abspar=abspar0
         empar=empar0
         weq=weq0
         nademflux=nademflux0
         for k=1,nsplit-1 do begin
            heipar=[heipar,scope_varfetch('heipar'+string(k,format='(I0)'))]
            abspar=[abspar,scope_varfetch('abspar'+string(k,format='(I0)'))]
            empar=[empar,scope_varfetch('empar'+string(k,format='(I0)'))]
            weq=[weq,scope_varfetch('weq'+string(k,format='(I0)'))]
            nademflux=[nademflux,scope_varfetch('nademflux'+string(k,format='(I0)'))]
         endfor
      
      endelse

      if keyword_set(dofirstemfit) then begin
         if j eq 0 then begin
            first_empar = empar
            first_weq = weq
            first_nademflux = nademflux
         endif
         if j eq 1 then begin
            empar = first_empar
            weq[*,1] = first_weq[*,1]
            nademflux = first_nademflux
         endif
      endif

   endfor

;  Save MC results to a file.
   nadmc = {hei: heipar, nadabs: abspar, nadem: empar}
   save,nadmc,file=outfile

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Compute errors from distributions and plot results
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   cgps_open,outplot,charsize=1,/inches,xs=7.5,ys=7.5,/qui,$
             /nomatch
   pos = cglayout([2,2],ixmar=[3,3],iymar=[4,4],oxmar=[0,0],oymar=[4,1],$
                   xgap=0,ygap=4,aspect=1d)
;
;  HeI
;
   for i=0,ncomp[0]-1 do begin
;     wavelength
      if ~ heifix[i,0] then begin
         bins=[]
         dat = heipar[*,0,i]
         cghistoplot,dat,pos=pos[*,0],missing=bad,ytit='',$
                     xtit='$\lambda$ ($\Angstrom$)',histdata=h,$
                     locations=loc,binsize=bins,xticks=2
         binc = loc + (bins/2d)
         yfit = mpfitpeak(binc,h,a,nterms=3,errors=sqrt(h))
         cgoplot,binc,yfit
         stat = ifsf_lohisig(dat)
         cgoplot,[stat[0],stat[0]],[0,max(h)]
         cgoplot,[stat[0]-stat[1],stat[0]-stat[1]],[0,max(h)],linesty=2
         cgoplot,[stat[0]+stat[2],stat[0]+stat[2]],[0,max(h)],linesty=2
         errors[i*3,*]=[stat[1],stat[2]]      
      endif else cghistoplot,heipar[*,0,i],pos=pos[*,0],missing=bad,ytit='',$
                             xtit='$\lambda$ ($\Angstrom$)',nbins=1,binsize=1d,$
                             xticks=2
;     sigma
      if ~ heifix[i,1] then begin
         bins=[]
         dat = heipar[*,1,i]
         cghistoplot,dat,pos=pos[*,1],missing=bad,/noerase,ytit='',$
                     xtit='$\sigma$ (km/s)',histdata=h,$
                     locations=loc,binsize=bins,xticks=2
         binc = loc + (bins/2d)
         yfit = mpfitpeak(binc,h,a,nterms=3,errors=sqrt(h))
         cgoplot,binc,yfit
         stat = ifsf_lohisig(dat)
         cgoplot,[stat[0],stat[0]],[0,max(h)]
         cgoplot,[stat[0]-stat[1],stat[0]-stat[1]],[0,max(h)],linesty=2
         cgoplot,[stat[0]+stat[2],stat[0]+stat[2]],[0,max(h)],linesty=2
         errors[i*3+1,*]=[stat[1],stat[2]]
      endif else cghistoplot,heipar[*,1,i],pos=pos[*,1],missing=bad,/noerase,$
                             ytit='',xtit='$\sigma$ (km/s)',nbins=1,binsize=1d,$
                             xticks=2
;     peak flux
      if ~ heifix[i,2] then begin
         bins=[]
         dat = heipar[*,2,i]
         cghistoplot,dat,pos=pos[*,2],missing=bad,/noerase,ytit='',$
                     xtit='F$\down\\lambda$ (peak)',histdata=h,$
                     locations=loc,binsize=bins,xticks=3
         binc = loc + (bins/2d)
         yfit = mpfitpeak(binc,h,a,nterms=3,errors=sqrt(h))
         cgoplot,binc,yfit
         stat = ifsf_lohisig(dat)
         cgoplot,[stat[0],stat[0]],[0,max(h)]
         cgoplot,[stat[0]-stat[1],stat[0]-stat[1]],[0,max(h)],linesty=2
         cgoplot,[stat[0]+stat[2],stat[0]+stat[2]],[0,max(h)],linesty=2
         errors[i*3+2,*]=[stat[1],stat[2]]
      endif else $
         cghistoplot,heipar[*,2,i],pos=pos[*,2],missing=bad,/noerase,ytit='',$
                     xtit='F$\down\\lambda$ (peak)',nbins=1,binsize=1d,xticks=3
      cgtext,0.5,0.98,'HeI em., c'+string(i+1,format='(I0)'),/norm
   endfor
;
;  NaD abs
;
   for i=0,ncomp[1]-1 do begin
;     covering factor x optical depth
      dat = abspar[*,0,i]*abspar[*,1,i]
      cghistoplot,abspar[*,0,i]*abspar[*,1,i],pos=pos[*,0],missing=bad*bad,$
                  ytit='',xtit='C$\down$f x $\tau$',xticks=2,histdata=h
      stat = ifsf_lohisig(dat)
      cgoplot,[stat[0],stat[0]],[0,max(h)]
      cgoplot,[stat[0]-stat[1],stat[0]-stat[1]],[0,max(h)],linesty=2
      cgoplot,[stat[0]+stat[2],stat[0]+stat[2]],[0,max(h)],linesty=2
      errors[ncomp[0]*3+i*4,*]=[stat[1]/2d,stat[2]/2d]
;     optical depth
      dat = abspar[*,1,i]
      cghistoplot,dat,pos=pos[*,1],/noerase,missing=bad,$
                  ytit='',xtit='$\tau$',xticks=2,histdata=h
      stat = ifsf_lohisig(dat)
      cgoplot,[stat[0],stat[0]],[0,max(h)]
      cgoplot,[stat[0]-stat[1],stat[0]-stat[1]],[0,max(h)],linesty=2
      cgoplot,[stat[0]+stat[2],stat[0]+stat[2]],[0,max(h)],linesty=2
      errors[ncomp[0]*3+i*4+1,*]=[stat[1],stat[2]]
;     wavelength
      bins=[]
      dat=abspar[*,2,i]
      cghistoplot,dat,pos=pos[*,2],missing=bad,/noerase,ytit='',$
                  xtit='$\lambda$ ($\Angstrom$)',xticks=2,histdata=h,$
                  locations=loc,binsize=bins
      binc = loc + (bins/2d)
      yfit = mpfitpeak(binc,h,a,nterms=3,errors=sqrt(h))
      cgoplot,binc,yfit
      stat = ifsf_lohisig(dat)
      cgoplot,[stat[0],stat[0]],[0,max(h)]
      cgoplot,[stat[0]-stat[1],stat[0]-stat[1]],[0,max(h)],linesty=2
      cgoplot,[stat[0]+stat[2],stat[0]+stat[2]],[0,max(h)],linesty=2
      errors[ncomp[0]*3+i*4+2,*]=[stat[1],stat[2]]
;     sigma
      bins=[]
      dat = abspar[*,3,i]
      cghistoplot,dat,pos=pos[*,3],missing=bad,/noerase,ytit='',$
                  xtit='$\sigma$ (km/s)',xticks=2,histdata=h,$
                  locations=loc,binsize=bins
      binc = loc + (bins/2d)
      yfit = mpfitpeak(binc,h,a,nterms=3,errors=sqrt(h))
      cgoplot,binc,yfit
      stat = ifsf_lohisig(dat)
      cgoplot,[stat[0],stat[0]],[0,max(h)]
      cgoplot,[stat[0]-stat[1],stat[0]-stat[1]],[0,max(h)],linesty=2
      cgoplot,[stat[0]+stat[2],stat[0]+stat[2]],[0,max(h)],linesty=2
      errors[ncomp[0]*3+i*4+3,*]=[stat[1],stat[2]]
;
      cgtext,0.5,0.98,'NaI D abs., c'+string(i+1,format='(I0)'),/norm
   endfor
;
;  NaD em
;
   for i=0,ncomp[2]-1 do begin
;     wavelength
      bins=[]
      if ~ nademfix[i,0] then begin
         dat = empar[*,0,i]
         cghistoplot,dat,pos=pos[*,0],missing=bad,ytit='',$
                     xtit='$\lambda$ ($\Angstrom$)',xticks=2,histdata=h,$
                     locations=loc,binsize=bins
         binc = loc + (bins/2d)
         yfit = mpfitpeak(binc,h,a,nterms=3,errors=sqrt(h))
         cgoplot,binc,yfit
         stat = ifsf_lohisig(dat)
         cgoplot,[stat[0],stat[0]],[0,max(h)]
         cgoplot,[stat[0]-stat[1],stat[0]-stat[1]],[0,max(h)],linesty=2
         cgoplot,[stat[0]+stat[2],stat[0]+stat[2]],[0,max(h)],linesty=2
         errors[ncomp[0]*3+ncomp[1]*4+i*4,*]=[stat[1],stat[2]]
      endif else $
         cghistoplot,empar[*,0,i],pos=pos[*,0],missing=bad,ytit='',$
                     xtit='$\lambda$ ($\Angstrom$)',xticks=2,nbins=1,binsize=1d
;     sigma
      bins=[]
      if ~ nademfix[i,1] then begin
         dat = empar[*,1,i]
         cghistoplot,dat,pos=pos[*,1],missing=bad,/noerase,ytit='',$
                     xtit='$\sigma$ (km/s)',xticks=2,histdata=h,$
                     locations=loc,binsize=bins
         binc = loc + (bins/2d)
         yfit = mpfitpeak(binc,h,a,nterms=3,errors=sqrt(h))
         cgoplot,binc,yfit
         stat = ifsf_lohisig(dat)
         cgoplot,[stat[0],stat[0]],[0,max(h)]
         cgoplot,[stat[0]-stat[1],stat[0]-stat[1]],[0,max(h)],linesty=2
         cgoplot,[stat[0]+stat[2],stat[0]+stat[2]],[0,max(h)],linesty=2
         errors[ncomp[0]*3+ncomp[1]*4+i*4+1,*]=[stat[1],stat[2]]
      endif else $
         cghistoplot,empar[*,1,i],pos=pos[*,1],missing=bad,/noerase,ytit='',$
                     xtit='$\sigma$ (km/s)',xticks=2,nbins=1,binsize=10d
;     peak flux
      bins=[]
      if ~ nademfix[i,2] then begin
         dat = empar[*,2,i]
         cghistoplot,dat,pos=pos[*,2],missing=bad,/noerase,ytit='',$
                     xtit='F$\down\\lambda$ (peak)',histdata=h,$
                     locations=loc,binsize=bins,xticks=3
         binc = loc + (bins/2d)
         yfit = mpfitpeak(binc,h,a,nterms=3,errors=sqrt(h))
         cgoplot,binc,yfit
         stat = ifsf_lohisig(dat)
         cgoplot,[stat[0],stat[0]],[0,max(h)]
         cgoplot,[stat[0]-stat[1],stat[0]-stat[1]],[0,max(h)],linesty=2
         cgoplot,[stat[0]+stat[2],stat[0]+stat[2]],[0,max(h)],linesty=2
         errors[ncomp[0]*3+ncomp[1]*4+i*4+2,*]=[stat[1],stat[2]]
      endif else $
         cghistoplot,empar[*,2,i],pos=pos[*,2],missing=bad,/noerase,ytit='',$
                     xtit='F$\down\\lambda$ (peak)',xticks=3,nbins=1,binsize=0.1d
;     flux ratio
      bins=[]
      if ~ nademfix[i,3] then begin
         dat = empar[*,3,i]
         cghistoplot,dat,pos=pos[*,3],missing=bad,/noerase,ytit='',$
                     xtit='F$\down2$/F$\down1$',xticks=2,histdata=h,$
                     locations=loc,binsize=bins
         binc = loc + (bins/2d)
         yfit = mpfitpeak(binc,h,a,nterms=3,errors=sqrt(h))
         cgoplot,binc,yfit
         stat = ifsf_lohisig(dat)
         cgoplot,[stat[0],stat[0]],[0,max(h)]
         cgoplot,[stat[0]-stat[1],stat[0]-stat[1]],[0,max(h)],linesty=2
         cgoplot,[stat[0]+stat[2],stat[0]+stat[2]],[0,max(h)],linesty=2
         errors[ncomp[0]*3+ncomp[1]*4+i*4+3,*]=[stat[1],stat[2]]
      endif else $
         cghistoplot,empar[*,3,i],pos=pos[*,3],missing=bad,/noerase,ytit='',$
                     xtit='F$\down2$/F$\down1$',xticks=2,nbins=1,binsize=0.1d
   endfor
;
; NaD equivalent widths and fluxes
;
   if ncomp[1] gt 0 OR ncomp[2] gt 0 then begin
;     absorption
      if ncomp[1] gt 0 then begin
         em_noerase=1
;        equivalent width
         bins=[]
         dat = weq[*,0]
         cghistoplot,dat,pos=pos[*,0],missing=bad,ytit='',$
                  xtit='W$\downeq$ (abs, $\Angstrom$)',xticks=2,histdata=h,$
                  locations=loc,binsize=bins
         binc = loc + (bins/2d)
         yfit = mpfitpeak(binc,h,a,nterms=3,errors=sqrt(h))
         cgoplot,binc,yfit
         stat = ifsf_lohisig(dat)
         cgoplot,[stat[0],stat[0]],[0,max(h)]
         cgoplot,[stat[0]-stat[1],stat[0]-stat[1]],[0,max(h)],linesty=2
         cgoplot,[stat[0]+stat[2],stat[0]+stat[2]],[0,max(h)],linesty=2
         weqerr[0,*]=[stat[1],stat[2]]
      endif else em_noerase=0
;     emission
      if ncomp[2] gt 0 then begin
;        equivalent width
         bins=[]
         dat = weq[*,1]
         cghistoplot,dat,pos=pos[*,2],missing=bad,noerase=em_noerase,ytit='',$
                  xtit='W$\downeq$ (em, $\Angstrom$)',xticks=2,histdata=h,$
                  locations=loc,binsize=bins
         binc = loc + (bins/2d)
         yfit = mpfitpeak(binc,h,a,nterms=3,errors=sqrt(h))
         cgoplot,binc,yfit
         stat = ifsf_lohisig(dat)
         cgoplot,[stat[0],stat[0]],[0,max(h)]
         cgoplot,[stat[0]-stat[1],stat[0]-stat[1]],[0,max(h)],linesty=2
         cgoplot,[stat[0]+stat[2],stat[0]+stat[2]],[0,max(h)],linesty=2
         weqerr[1,*]=[stat[1],stat[2]]
;        flux
         bins=[]
         dat = nademflux
         cghistoplot,dat,pos=pos[*,3],missing=bad,/noerase,ytit='',$
                  xtit='Flux (em)',xticks=2,histdata=h,$
                  locations=loc,binsize=bins
         binc = loc + (bins/2d)
         yfit = mpfitpeak(binc,h,a,nterms=3,errors=sqrt(h))
         cgoplot,binc,yfit
         stat = ifsf_lohisig(dat)
         cgoplot,[stat[0],stat[0]],[0,max(h)]
         cgoplot,[stat[0]-stat[1],stat[0]-stat[1]],[0,max(h)],linesty=2
         cgoplot,[stat[0]+stat[2],stat[0]+stat[2]],[0,max(h)],linesty=2
         nademfluxerr=[stat[1],stat[2]]
      endif
      
   endif
   
   cgps_close,/pdf,/delete_ps


finish:

   return,errors

end
