; docformat = 'rst'
;
;+
;
; Fit Na D absorption line. Plot fit and write fit parameters to a file.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Produces a plot ([gal_col_row]_nad_fit.jpg) for each fit, and a
;    parameter file ([gal].nad.dat) containing the fit parameters for
;    all spaxels.
;
; :Params:
;    initproc: in, required, type=string
;      Name of procedure to initialize the fit.
;
; :Keywords:
;    cols: in, optional, type=intarr, default=all
;      Columns to fit, in 1-offset format. Either a scalar or a
;      two-element vector listing the first and last columns to fit.
;    rows: in, optional, type=intarr, default=all
;      Rows to fit, in 1-offset format. Either a scalar or a
;      two-element vector listing the first and last rows to fit.
;    noerr: in, optional, type=byte
;      Do not Monte Carlo the errors.
;    nomc: in, optional, type=byte
;      Do not re-run Monte Carlo error simulations, but do grab old outputs
;      and append them to the output structure.
;    noplot: in, optional, type=byte
;      Do not produce plots.
;    noxdr: in, optional, type=byte
;      Do no write XDR file with fit parameters.
;    nsplit: in, optional, type=integer, default=1
;      Number of independent child processes into which the Monte Carlo 
;      computation of errors is split.
;    verbose: in, optional, type=byte
;      Print error and progress messages. Propagates to most/all
;      subroutines.
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
;      2010jul21, DSNR, created
;      2013nov22, DSNR, renamed, added license and copyright
;      2014may09, DSNR, completely re-written to use MPFIT 
;      2014may30, DSNR, updated to allow use without previous emission-line fit
;                       with IFSF
;      2014jun10, DSNR, compute equivalent widths and print NaD parameters to 
;                       XDR file
;      2014jul07, DSNR, added error computation
;      2014jul22, DSNR, added option to not re-run MC simulations but grab old
;                       outputs
;    
; :Copyright:
;    Copyright (C) 2013-2014 David S. N. Rupke
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
pro ifsf_fitnad,initproc,cols=cols,rows=rows,nsplit=nsplit,verbose=verbose,$
                noerr=noerr,noplot=noplot,noxdr=noxdr,nomc=nomc

   bad = 1d99
   tratio = 2.005d ; NaD optical depth ratio (blue to red)
   nad_emrat_init = 1.5d

   starttime = systime(1)
   time = 0
   if ~ keyword_set(nsplit) then nsplit=1
   if keyword_set(verbose) then quiet=0 else quiet=1
   if ~ keyword_set(noplot) then noplot=0 else noplot=1

   ; Get fit initialization
   initnad={dumy: 1}
   initdat = call_function(initproc,initnad=initnad)
   maxncomp = initnad.maxncomp


   ; Get linelist
   linelist = ifsf_linelist(['NaD1','NaD2','HeI5876'])

   if tag_exist(initnad,'taumax') then taumax=initnad.taumax $
   else taumax = 5d

   ifsf_printnadpar,nadparlun,outfile=initdat.outdir+initdat.label+'.nad.dat'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if ~ tag_exist(initnad,'noemlinfit') then $
      restore,file=initdat.outdir+initdat.label+'.lin.xdr'
   restore,file=initdat.outdir+initdat.label+'.nadspec.xdr'
   nadsize = size(nadcube.wave)

   ncols = nadsize[1]
   nrows = nadsize[2]
   nz = nadsize[3]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Loop through spaxels
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  Switch to track when first fit done
   firstfit=1

   if not keyword_set(cols) then cols=[1,ncols] $
   else if n_elements(cols) eq 1 then cols = [cols,cols]
   cols = fix(cols)
   for i=cols[0]-1,cols[1]-1 do begin

      print,'Column ',i+1,' of ',ncols,format='(A,I0,A,I0)'

      if not keyword_set(rows) then rows=[1,nrows] $
      else if n_elements(rows) eq 1 then rows = [rows,rows]
      rows = fix(rows)
      for j=rows[0]-1,rows[1]-1 do begin

         print,'  Row ',j+1,' of ',nrows,format='(A,I0,A,I0)'

         if ~ tag_exist(initnad,'noemlinfit') then begin
;           Restore continuum + emission-line fit
            lab = string(i+1,'_',j+1,format='(I04,A,I04)')
            infile = initdat.outdir+initdat.label+'_'+lab+'.xdr'
            outfile = initdat.outdir+initdat.label+'_'+lab
            if ~ file_test(infile) then begin
               print,'IFSF_FITNAD: No XDR file for ',i+1,', ',j+1,'.',$
                     format='(A0,I4,A0,I4,A0)'
               goto,nofit
            endif
            restore,file=infile
         endif

;        Get HeI parameters
         refline = initnad.heitie[i,j]
         heifix=0
;        Cases of fitting HeI
         if refline ne '' then begin
;           Case of constraining with another line
            if refline ne 'HeI5876' AND $
               ~ tag_exist(initnad,'noemlinfit') then begin
;              Get reference line parameters
               reflinelist = ifsf_linelist([refline])
;              Case of constraining with another line from another spaxel
               if tag_exist(initnad,'heitiecol') AND $
                  tag_exist(initnad,'heitierow') then begin
                  struct_tmp = struct
;                 Restore continuum + emission-line fit for reference spaxel
                  refi = initnad.heitiecol[i,j]-1
                  refj = initnad.heitierow[i,j]-1
                  lab = string(refi+1,'_',refj+1,format='(I04,A,I04)')
                  infile = initdat.outdir+initdat.label+'_'+lab+'.xdr'
                  if ~ file_test(infile) then begin
                     print,'IFSF_FITNAD: No XDR file for ',i+1,', ',j+1,'.',$
                           format='(A0,I4,A0,I4,A0)'
                     goto,nofit
                  endif
                  restore,file=infile
                  linepars = ifsf_sepfitpars(reflinelist,struct.param,$
                                             struct.perror,struct.parinfo)
                  struct = struct_tmp
;              Case of constraining with another line from the same spaxel
               endif else begin
                  linepars = ifsf_sepfitpars(reflinelist,struct.param,$
                                             struct.perror,struct.parinfo)                  
               endelse
               iem = where(linepars.flux[refline,*] ne 0d,nhei)
               if nhei gt 0 then $
                  inithei = [[(linepars.wave)[refline,0:nhei-1]/$
                              reflinelist[refline]*linelist['HeI5876']],$
                             [(linepars.sigma)[refline,0:nhei-1]],$
                             [dblarr(nhei)+0.1d]] $
               else inithei=0
               heifix = bytarr(nhei,3)
               heifix[*,0:1] = 1b              
;           Case of allowing the line to vary freely
            endif else if tag_exist(initnad,'hei_zinit') AND $
                          tag_exist(initnad,'hei_siginit') AND $
                          tag_exist(initnad,'nhei') then begin
               nhei=initnad.nhei[i,j]
               inithei = [[((initnad.hei_zinit)[i,j,0:nhei-1]+1d)*$
                           linelist['HeI5876']],$
                          [(initnad.hei_siginit)[i,j,0:nhei-1]],$
                          [dblarr(nhei)+0.1d]]
               heifix = bytarr(nhei,3)
            endif else begin
               print,'IFSF_FITNAD: HeI5876 initialization parameters'+$
                     ' not properly specified.'
               goto,nofit
            endelse
;        Case of no HeI fit
         endif else begin
            nhei=0
            inithei=0
         endelse

;        Get NaD absorption parameters
         nnadabs = initnad.nnadabs[i,j]
         if nnadabs gt 0 then begin
            if tag_exist(initnad,'nadabs_cfinit') then $
               cfinit = reform((initnad.nadabs_cfinit)[i,j,0:nnadabs-1],nnadabs) $
            else cfinit = dblarr(nnadabs)+0.5d
            if tag_exist(initnad,'nadabs_tauinit') then $
               tauinit = reform((initnad.nadabs_tauinit)[i,j,0:nnadabs-1],nnadabs) $
            else tauinit = dblarr(nnadabs)+0.5d
            winit = reform(((initnad.nadabs_zinit)[i,j,0:nnadabs-1]$
                            +1d)*linelist['NaD1'],nnadabs)
            siginit = reform((initnad.nadabs_siginit)$
                              [i,j,0:nnadabs-1],nnadabs)
            initnadabs = [[cfinit],[tauinit],[winit],[siginit]]
         endif else initnadabs=0

;        Get NaD emission parameters
         nnadem = initnad.nnadem[i,j]
;        placeholders for case of separately fitting emission and absorption
         dofirstemfit=0b
         first_nademfix=0d
         first_parinit=0d
         first_modflux=0d
         if nnadem gt 0 then begin
               
            winit = reform(((initnad.nadem_zinit)[i,j,0:nnadem-1]+1d)$
                           *linelist['NaD1'],nnadem)
            siginit = reform((initnad.nadem_siginit)[i,j,0:nnadem-1],nnadem)
            if tag_exist(initnad,'nadem_finit') then $
               finit = (initnad.nadem_finit)[i,j,0:nnadem-1] $
            else finit = dblarr(nnadem)+0.1d
            if tag_exist(initnad,'nadem_rinit') then $
               rinit = (initnad.nadem_rinit)[i,j,0:nnadem-1] $
            else rinit = dblarr(nnadem)+nad_emrat_init
            initnadem = [[winit],[siginit],[finit],[rinit]]
            if tag_exist(initnad,'nadem_fix') then $
               nademfix=reform((initnad.nadem_fix)[i,j,0:nnadem-1,*],nnadem,4) $
            else nademfix=0b
            
;           Fit the emission line region only if requested, and if NaD emission 
;           and absorption (or NaD emission and HeI emission) are to be fit.
            if tag_exist(initnad,'nadem_fitinit') AND $
               (nnadabs gt 0 OR nhei gt 0) then begin

               dofirstemfit=1b
               
;              Fill out parameter structure with initial guesses and constraints
               if tag_exist(initnad,'argsinitpar') then parinit = $
                  call_function(initnad.fcninitpar,0,0,initnadem,$
                                initnad.nadabs_siglim,initnad.nadem_siglim,$
                                heifix=0,nademfix=nademfix,$
                                _extra=initnad.argsinitpar) $
               else parinit = $
                  call_function(initnad.fcninitpar,0,0,initnadem,$
                                initnad.nadabs_siglim,initnad.nadem_siglim,$
                                heifix=0,nademfix=nademfix)

;              This block looks for automatically-determined absorption and 
;              emission line indices (from IFSF_CMPNADWEQ, invoked by IFSFA)
;              and uses these to only fit the emission line region by setting
;              anything blueward to 1. Note that if only an emission line was 
;              found, then the index used is shifted blueward slightly to make
;              sure the entire line is included.   
               tmpdat = (nadcube.dat)[i,j,*]
               if (nadcube.iweq)[i,j,1] ne -1 then $
                  tmpind = (nadcube.iweq)[i,j,1] $
               else if (nadcube.iweq)[i,j,2] ne -1 then $
                  tmpind = (nadcube.iweq)[i,j,2]-1 $
               else begin
                  print,'IFSF_FITNAD: No absorption or emission indices. Aborting.'                  
                  goto,finish
               endelse
               tmpdat[0:tmpind] = 1d
                                
               param = Mpfitfun(initnad.fcnfitnad,$
                                (nadcube.wave)[i,j,*],$
                                tmpdat,$
                                (nadcube.err)[i,j,*],$
                                parinfo=parinit,perror=perror,maxiter=100,$
                                bestnorm=chisq_emonly,covar=covar,$
                                yfit=specfit,dof=dof_emonly,$
                                nfev=nfev,niter=niter_emonly,status=status,$
                                quiet=quiet,$
                                npegged=npegged,ftol=1D-6,errmsg=errmsg)
               if status eq 5 then print,'IFSF_FITNAD: Max. iterations reached.'
               if status eq 0 OR status eq -16 then begin
                  print,'IFSF_FITNAD: Error in MPFIT. Aborting.'
                  goto,finish
               endif

               first_nademfix = nademfix
               first_parinit = parinit
               first_modflux = tmpdat
               initnadem = reform(param[3:2+nnadem*4],nnadem,4)
               nademfix=rebin(transpose([1b,1b,1b,1b]),nnadem,4)
                           
            endif
         endif else begin
            initnadem=0
            nademfix=0b
         endelse

;        Fill out parameter structure with initial guesses and constraints

         if (nnadem eq 0 AND nnadabs eq 0 AND nhei eq 0) then begin
            print,'IFSF_FITNAD: No components specified. Skipping this spaxel.'
            goto,nofit
         endif

         if nnadabs gt 0 AND tag_exist(initnad,'nadabs_fix') then $
            nadabsfix=reform((initnad.nadabs_fix)[i,j,0:nnadabs-1,*],nnadabs,4) $
         else nadabsfix=0b

         if tag_exist(initnad,'argsinitpar') then parinit = $
            call_function(initnad.fcninitpar,inithei,initnadabs,initnadem,$
                          initnad.nadabs_siglim,initnad.nadem_siglim,$
                          heifix=heifix,nadabsfix=nadabsfix,nademfix=nademfix,$
                          _extra=initnad.argsinitpar) $
         else parinit = $
            call_function(initnad.fcninitpar,inithei,initnadabs,initnadem,$
                          initnad.nadabs_siglim,initnad.nadem_siglim,$
                          heifix=heifix,nadabsfix=nadabsfix,nademfix=nademfix)

         param = Mpfitfun(initnad.fcnfitnad,$
                          (nadcube.wave)[i,j,*],$
                          (nadcube.dat)[i,j,*],$
                          (nadcube.err)[i,j,*],$
                          parinfo=parinit,perror=perror,maxiter=100,$
                          bestnorm=chisq,covar=covar,yfit=specfit,dof=dof,$
                          nfev=nfev,niter=niter,status=status,quiet=quiet,$
                          npegged=npegged,ftol=1D-6,errmsg=errmsg)
         if status eq 5 then print,'IFSF_FITNAD: Max. iterations reached.'
         if status eq 0 OR status eq -16 then begin
            print,'IFSF_FITNAD: Error in MPFIT. Aborting.'
            goto,finish
         endif

         if ~ tag_exist(initnad,'noemlinfit') then zuse = struct.zstar $
         else zuse = initdat.zsys_gas
         if ~ noplot then begin
            if tag_exist(initnad,'argspltfitnad') then $
               ifsf_pltnadfit,(nadcube.wave)[i,j,*],$
                              (nadcube.dat)[i,j,*],$
                              param,outfile+'_nad_fit',struct.zstar,$
                              _extra=initnad.argspltfitnad $
            else $
               ifsf_pltnadfit,(nadcube.wave)[i,j,*],$
                              (nadcube.dat)[i,j,*],$
                              param,outfile+'_nad_fit',struct.zstar
         endif

;        Compute model equivalent widths
         weq=1
         nademflux=1
         modspec = ifsf_nadfcn((nadcube.wave)[i,j,*],param,weq=weq,$
                               nademflux=nademflux,cont=(nadcube.cont)[i,j,*])

;        Compute errors in fit
         if ~ keyword_set(noerr) then begin
            if dofirstemfit then nademfix_use = first_nademfix $
            else nademfix_use = nademfix
            if keyword_set(nomc) then plotonly=1b else plotonly=0b
            errors = ifsf_fitnaderr([nhei,nnadabs,nnadem],(nadcube.wave)[i,j,*],$
                                    modspec,(nadcube.err)[i,j,*],$
                                    (nadcube.cont)[i,j,*],parinit,$
                                    outfile+'_nad_errs.ps',outfile+'_nad_mc.xdr',$
                                    dofirstemfit=dofirstemfit,$
                                    first_parinit=first_parinit,$
                                    first_modflux=first_modflux,$
                                    heifix=heifix,nadabsfix=nadabsfix,$
                                    nademfix=nademfix_use,$
                                    niter=initnad.mcniter,$
                                    nsplit=nsplit,quiet=quiet,weqerr=weqerr,$
                                    nademfluxerr=nademfluxerr,noplot=noplot,$
                                    plotonly=plotonly)
         endif
         
         ifsf_printnadpar,nadparlun,i+1,j+1,param

         if firstfit then begin
            nadfit = $
               {chisq: dblarr(ncols,nrows)+bad,$
                chisq_emonly: dblarr(ncols,nrows)+bad,$
                dof: dblarr(ncols,nrows)+bad,$
                niter: dblarr(ncols,nrows)+bad,$
;                
                weqabs: dblarr(ncols,nrows,1+maxncomp)+bad,$
                weqabserr: dblarr(ncols,nrows,2)+bad,$
                weqem: dblarr(ncols,nrows,1+maxncomp)+bad,$
                weqemerr: dblarr(ncols,nrows,2)+bad,$
                totfluxem: dblarr(ncols,nrows,1+maxncomp)+bad,$
                totfluxemerr: dblarr(ncols,nrows,2)+bad,$
                cf: dblarr(ncols,nrows,maxncomp)+bad,$
                cferr: dblarr(ncols,nrows,maxncomp,2)+bad,$
                tau: dblarr(ncols,nrows,maxncomp)+bad,$
                tauerr: dblarr(ncols,nrows,maxncomp,2)+bad,$
                waveabs: dblarr(ncols,nrows,maxncomp)+bad,$
                waveabserr: dblarr(ncols,nrows,maxncomp,2)+bad,$
                sigmaabs: dblarr(ncols,nrows,maxncomp)+bad,$
                sigmaabserr: dblarr(ncols,nrows,maxncomp,2)+bad,$
                waveem: dblarr(ncols,nrows,maxncomp)+bad,$
                waveemerr: dblarr(ncols,nrows,maxncomp,2)+bad,$
                sigmaem: dblarr(ncols,nrows,maxncomp)+bad,$
                sigmaemerr: dblarr(ncols,nrows,maxncomp,2)+bad,$
                flux: dblarr(ncols,nrows,maxncomp)+bad,$
                fluxerr: dblarr(ncols,nrows,maxncomp,2)+bad,$
                frat: dblarr(ncols,nrows,maxncomp)+bad,$
                fraterr: dblarr(ncols,nrows,maxncomp,2)+bad}
            firstfit = 0
         endif
         nadfit.chisq[i,j]=chisq
         if dofirstemfit then nadfit.chisq_emonly[i,j]=chisq_emonly
         nadfit.dof[i,j]=dof
         nadfit.niter[i,j]=niter
         nadfit.weqabs[i,j,0:nnadabs]=weq.abs
         nadfit.weqem[i,j,0:nnadem]=weq.em
         nadfit.totfluxem[i,j,0:nnadem]=nademflux
         if ~ keyword_set(noerr) then begin
            nadfit.weqabserr[i,j,*]=weqerr[0,*]
            nadfit.weqemerr[i,j,*]=weqerr[1,*]
            nadfit.totfluxemerr[i,j,*]=nademfluxerr
         endif
         if nnadabs gt 0 then begin
            iarr = 3+nhei*3 + dindgen(nnadabs)*4
            nadfit.cf[i,j,0:nnadabs-1]=param[iarr]
            nadfit.tau[i,j,0:nnadabs-1]=param[iarr+1]
            nadfit.waveabs[i,j,0:nnadabs-1]=param[iarr+2]
            nadfit.sigmaabs[i,j,0:nnadabs-1]=param[iarr+3]
            if ~ keyword_set(noerr) then begin
               nadfit.cferr[i,j,0:nnadabs-1,*]=errors[iarr-3,*]
               nadfit.tauerr[i,j,0:nnadabs-1,*]=errors[iarr-3+1,*]
               nadfit.waveabserr[i,j,0:nnadabs-1,*]=errors[iarr-3+2,*]
               nadfit.sigmaabserr[i,j,0:nnadabs-1,*]=errors[iarr-3+3,*]
            endif
         endif
         if nnadem gt 0 then begin
            iarr = 3+nhei*3+nnadabs*4 + dindgen(nnadem)*4
            nadfit.waveem[i,j,0:nnadem-1]=param[iarr]
            nadfit.sigmaem[i,j,0:nnadem-1]=param[iarr+1]
            nadfit.flux[i,j,0:nnadem-1]=param[iarr+2]
            nadfit.frat[i,j,0:nnadem-1]=param[iarr+3]
            if ~ keyword_set(noerr) then begin
               nadfit.waveemerr[i,j,0:nnadem-1,*]=errors[iarr-3,*]
               nadfit.sigmaemerr[i,j,0:nnadem-1,*]=errors[iarr-3+1,*]
               nadfit.fluxerr[i,j,0:nnadem-1,*]=errors[iarr-3+2,*]
               nadfit.fraterr[i,j,0:nnadem-1,*]=errors[iarr-3+3,*]               
            endif
         endif

nofit:

      endfor
      
   endfor


finish:

   if ~ keyword_set(noxdr) then begin
      if tag_exist(initnad,'outxdr') then outxdr=initnad.outxdr $
      else outxdr=initdat.outdir+initdat.label+'.nadfit.xdr'
      save,nadfit,file=outxdr
   endif

   free_lun,nadparlun

   print,'Total runtime: ',systime(1)-starttime,' s.',$
         format='(/,A0,I0,A0,/)'

end
