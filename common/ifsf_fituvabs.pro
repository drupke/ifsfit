; docformat = 'rst'
;
;+
;
; Fit O VI and N V absorption lines. Plot fit and write fit parameters to a file.
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
;    initfile: in, required, type=string
;      Name of file with two columns; first column is a list of galaxies 
;      and second column is a list of the corresponding redshift.
;      
;    directoryname: in, required, type=string
;      Directory path where spectra data is located and where results are to outputted.
;      
;    galaxyname: in, required, type=string
;      Name of target galaxy to fit.
;      
;    doublet: in, required, type=string
;      Doublet profiles to fit. Currently can fit OVI and NV doublets.
;      Captilization is required.
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
;    weights: in, optional, type=string
;      Changes fitting procedure to automatically weight each value. 
;      If not used, instead uses the error spectra to produce weights
; 
; :Authors:
;    David S. N. Rupke::
;      Rhodes College
;      Department of Physics
;      2000 N. Parkway
;      Memphis, TN 38104
;      drupke@gmail.com
;      
;    Anthony To::
;      Rhodes College
;      Department of Physics
;      2000 N. Parkway
;      Memphis, TN 38104
;      andto94@gmail.com
;
; :History:
;    ChangeHistory::
;      2015jul01, Begin work on adapting from David Rupke's ifs_fitnad procedure.
;    
; :Copyright:
;    Copyright (C) 2013-2015 David S. N. Rupke
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
pro ifsf_fituvabs,initfile,directoryname,galaxyname,doublet,cols=cols,rows=rows,nsplit=nsplit,verbose=verbose,$
                noerr=noerr,noplot=noplot,noxdr=noxdr,nomc=nomc,weights=weights

   bad = 1d99
   tratio = 2.005d ; UVAbs optical depth ratio (blue to red)
   nad_emrat_init = 1.5d
   c = 299792.458d

   starttime = systime(1)
   time = 0
   if ~ keyword_set(nsplit) then nsplit=1
   if keyword_set(verbose) then quiet=0 else quiet=1
   if ~ keyword_set(noplot) then noplot=0 else noplot=1

   ; Get fit initialization
   initnad={dumy: 1}
   directoryname+='/'
   readcol, directoryname+initfile, galaxynamelist,redshiftlist, SKIPLINE=1,format = '(A,D)'
   Galaxiestoloop = FILE_LINES(directoryname+initfile)-1
   selectionparameter=WHERE(galaxynamelist eq galaxyname)
   galaxyname=galaxynamelist[selectionparameter[0]]
   redshift=redshiftlist[selectionparameter]
   readcol, directoryname+galaxyname+'/'+galaxyname+doublet+'param', profileshifts, profilesig, coveringfactor, opticaldepth, FORMAT='(A,D,D,D,D)'
   initproc = 'ifsf_'+galaxyname+doublet
   comps=N_ELEMENTS(profilesig)
   initdat = call_function(initproc,directoryname, galaxyname, redshift, profileshifts, profilesig, coveringfactor, $
    opticaldepth, initnad=initnad)

   maxncomp = initnad.maxncomp
   
   ; Get linelist
   linelist = ifsf_linelist(['[OVI1]1032','[OVI2]1038','[LyB]1026','[LyA]1216','[NV1]1239','[NV2]1243'])

   if tag_exist(initnad,'taumax') then taumax=initnad.taumax $
   else taumax = 5d

   ext='.txt'
   
   if tag_exist(initnad,'nadem_siglim') then $fr
      nadem_siglim = initnad.nadem_siglim $
   else nadem_siglim = 0d

;   ifsf_printnadpar,nadparlun,outfile=initdat.outdir+initdat.label+'.txt'
;    print, initdat.outdir+initdat.label+'.txt'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;   if ~ tag_exist(initnad,'noemlinfit') then $
;      restore,file=initdat.outdir+initdat.label+'.txt'
;   restore,file=initdat.outdir+initdat.label+'.txt'
;   nadsize = size(nadcube.wave)



   ncols = 1
   nrows = 1
   oned=1
   nz = 1
   nadcube = {err: initnad.error, $
    dat : initnad.relativeflux, $
    wave : initnad.wavelength, $ 
    cont: initnad.continuum, $
    flux: initnad.flux $   
    }
;   nadcube.err=[1,1,initnad.error]
;   nadcube.dat=[1,1,initnad.relativeflux]
;   nadcube.wave=[1,1,initnad.wavelength]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Loop through spaxels
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  Switch to track when first fit done
   firstfit=1
   nhei = 1
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
            if oned then lab = string(i+1,format='(I04)') $
            else lab = string(i+1,'_',j+1,format='(I04,A,I04)')
            infile = initdat.outdir+initdat.label+ext
            print, infile
            outfile = initdat.outdir+initdat.label+'_'+lab
            if ~ file_test(infile) then begin
               print,'IFSF_FITUVABS: No file for ',i+1,', ',j+1,'.',$
                     format='(A0,I4,A0,I4,A0)'
               goto,nofit
            endif
;            restore,file=infile
         endif

;        Get UVAbs absorption parameters
         nnadabs = initnad.nnadabs[i,j]
         IF (doublet eq 'OVI') THEN BEGIN
          if nnadabs gt 0 then begin
            if tag_exist(initnad,'nadabs_cfinit') then BEGIN
              cfinit=initnad.nadabs_cfinit
;               cfinit = reform((initnad.nadabs_cfinit)[i,j,0:nnadabs-1],nnadabs) $
            ENDIF else cfinit = dblarr(nnadabs)+0.5d
            if tag_exist(initnad,'nadabs_tauinit') then BEGIN
              tauinit=initnad.nadabs_tauinit
;               tauinit = reform((initnad.nadabs_tauinit)[i,j,0:nnadabs-1],nnadabs) $
            ENDIF else tauinit = dblarr(nnadabs)+0.5d
            winit = reform(((initnad.nadabs_zinit)[i,j,0:nnadabs-1]$
                            +1d)*linelist['[OVI1]1032'],nnadabs)
            siginit = reform((initnad.nadabs_siginit)$
                                 [i,j,0:nnadabs-1],nnadabs)
            initnadabs = [[cfinit],[tauinit],[winit],[siginit]]
          endif else initnadabs=0
         ENDIF
         IF (doublet eq 'NV') THEN BEGIN
           if nnadabs gt 0 then begin
             if tag_exist(initnad,'nadabs_cfinit') then BEGIN
               cfinit=initnad.nadabs_cfinit
               ;               cfinit = reform((initnad.nadabs_cfinit)[i,j,0:nnadabs-1],nnadabs) $
             ENDIF else cfinit = dblarr(nnadabs)+0.5d
             if tag_exist(initnad,'nadabs_tauinit') then BEGIN
               tauinit=initnad.nadabs_tauinit
               ;               tauinit = reform((initnad.nadabs_tauinit)[i,j,0:nnadabs-1],nnadabs) $
             ENDIF else tauinit = dblarr(nnadabs)+0.5d
             winit = reform(((initnad.nadabs_zinit)[i,j,0:nnadabs-1]$
               +1d)*linelist['[NV1]1239'],nnadabs)
             siginit = reform((initnad.nadabs_siginit)$
               [i,j,0:nnadabs-1],nnadabs)
             initnadabs = [[cfinit],[tauinit],[winit],[siginit]]
           endif else initnadabs=0
         ENDIF
         
         IF (doublet eq 'NV') THEN BEGIN
          linename = '[NV1]1239'
         ENDIF
         IF (doublet eq 'OVI') THEN BEGIN
          linename = '[OVI1]1032'
         ENDIF
;        Get UVAbs emission parameters
         nnadem = initnad.nnadem[i,j]
;        placeholders for case of separately fitting emission and absorption
         dofirstemfit=0b
         first_nademfix=0d
         first_parinit=0d
         first_modflux=0d
         if nnadem gt 0 then begin
               
            if nadem_siglim[0] eq 0 then print,'IFSF_FITUVABS: ERROR: Emission '+$
               'line sigma limits not set (UVABSEM_SIGLIM) in INITUVABS '+$
               'structure, but emission lines need to be fit.'
               
            winit = reform(((initnad.nadem_zinit)[i,j,0:nnadem-1]+1d)$
                           *linelist[linename],nnadem)
            siginit = reform((initnad.nadem_siginit)[i,j,0:nnadem-1],nnadem)
            if tag_exist(initnad,'nadem_finit') then $
               finit = reform((initnad.nadem_finit)[i,j,0:nnadem-1],nnadem) $
            else finit = dblarr(nnadem)+0.1d
            if tag_exist(initnad,'nadem_rinit') then $
               rinit = reform((initnad.nadem_rinit)[i,j,0:nnadem-1],nnadem) $
            else rinit = dblarr(nnadem)+nad_emrat_init
            initnadem = [[winit],[siginit],[finit],[rinit]]
            if tag_exist(initnad,'nadem_fix') then $
               nademfix=reform((initnad.nadem_fix)[i,j,0:nnadem-1,*],nnadem,4) $
            else nademfix=0b
            
;           Fit the emission line region only if requested, and if UVAbs emission 
;           and absorption (or UVAbs emission and HeI emission) are to be fit.
;            if tag_exist(initnad,'nadem_fitinit') AND $
;               (nnadabs gt 0 OR nhei gt 0) then begin
;
;               dofirstemfit=1b
;               
;;              Fill out parameter structure with initial guesses and constraints
;               if tag_exist(initnad,'argsinitpar') then parinit = $
;                  call_function(initnad.fcninitpar,0,0,initnadem,$
;                                initnad.nadabs_siglim,nadem_siglim,$
;                                heifix=0,nademfix=nademfix,$
;                                _extra=initnad.argsinitpar) $
;               else parinit = $
;                  call_function(initnad.fcninitpar,0,0,initnadem,$
;                                initnad.nadabs_siglim,nadem_siglim,$
;                                heifix=0,nademfix=nademfix)
;
;;              This block looks for automatically-determined absorption and 
;;              emission line indices (from IFSF_CMPUVABSWEQ, invoked by IFSFA)
;;              and uses these to only fit the emission line region by setting
;;              anything blueward to 1. Note that if only an emission line was 
;;              found, then the index used is shifted blueward slightly to make
;;              sure the entire line is included.   
;               tmpdat = (nadcube.dat)[i,j,*]
;;               if (nadcube.iweq)[i,j,1] ne -1 then $
;;                  tmpind = (nadcube.iweq)[i,j,1] $
;;               else if (nadcube.iweq)[i,j,2] ne -1 then $
;;                  tmpind = (nadcube.iweq)[i,j,2]-1 $
;;               else begin
;;                  print,'IFSF_FITUVABS: No absorption or emission indices. Aborting.'                  
;;                  goto,finish
;;               endelse
;;               tmpdat[0:tmpind] = 1d
;                                
;               param = Mpfitfun(initnad.fcnfitnad,$
;                                (nadcube.wave)[i,j,*],$
;                                tmpdat,$
;                                (nadcube.err)[i,j,*],$
;                                parinfo=parinfo,perror=perror,maxiter=100,$
;                                bestnorm=chisq_emonly,covar=covar,$
;                                yfit=specfit,dof=dof_emonly,$
;                                nfev=nfev,niter=niter_emonly,status=status,$
;                                quiet=quiet,$
;                                npegged=npegged,ftol=1D-6,errmsg=errmsg)
;               if status eq 5 then print,'IFSF_FITuvabs: Max. iterations reached.'
;               if status eq 0 OR status eq -16 then begin
;                  print,'IFSF_FITuvabs: Error in MPFIT. Aborting.'
;                  goto,finish
;               endif
;
;               first_UVABSemfix = nademfix
;               first_parinit = parinit
;               first_modflux = tmpdat
;               initnadem = transpose(reform(param[3:2+nnadem*4],4,nnadem))
;               nademfix=rebin(transpose([1b,1b,1b,1b]),nnadem,4)
;                           
;            endif
         endif else begin
            initnadem=0
            nademfix=0b
         endelse

;        Fill out parameter structure with initial guesses and constraints

         if nnadabs gt 0 AND tag_exist(initnad,'nadabs_fix') then $
            nadabsfix=reform((initnad.nadabs_fix)[i,j,0:nnadabs-1,*],nnadabs,4) $
         else nadabsfix=0b

         if tag_exist(initnad,'argsinitpar') then parinit = $
            call_function(initnad.fcninitpar,initnadabs,initnadem,$
                          initnad.nadabs_siglim,nadem_siglim,$
                          nadabsfix=nadabsfix,nademfix=nademfix) $
         else parinit = $
            call_function(initnad.fcninitpar,initnadabs,initnadem,$
                          initnad.nadabs_siglim,nadem_siglim,$
                          nadabsfix=nadabsfix,nademfix=nademfix)
         if (~ keyword_set(weights)) then begin
           param = Mpfitfun(initnad.fcnfitnad,$
                            (nadcube.wave)[initnad.fitindex[0]:initnad.fitindex[1]],$
                            (nadcube.dat)[initnad.fitindex[0]:initnad.fitindex[1]],$
                            (nadcube.err)[initnad.fitindex[0]:initnad.fitindex[1]],$
                            parinfo=parinit,perror=perror,maxiter=500,$
                            bestnorm=chisq,covar=covar,yfit=specfit,dof=dof,$
                            nfev=nfev,niter=niter,status=status,quiet=quiet,$
                            npegged=npegged,ftol=1D-6,errmsg=errmsg)
         endif else begin
           param = Mpfitfun(initnad.fcnfitnad,$
                            (nadcube.wave)[initnad.fitindex[0]:initnad.fitindex[1]],$
                            (nadcube.dat)[initnad.fitindex[0]:initnad.fitindex[1]],$
                            (nadcube.err)[initnad.fitindex[0]:initnad.fitindex[1]],$
                            parinfo=parinit,perror=perror,maxiter=500,$
                            bestnorm=chisq,covar=covar,yfit=specfit,dof=dof,$
                            nfev=nfev,niter=niter,status=status,quiet=quiet,$
                            npegged=npegged,ftol=1D-6,errmsg=errmsg,/WEIGHTS)          
         endelse
         if status eq 5 then print,'IFSF_FITUVABS: Max. iterations reached.'
         if status eq 0 OR status eq -16 then begin
            print,'IFSF_FITUVABS: Error in MPFIT. Aborting.'
            goto,finish
         endif

;        Plot fit
;        
;        If the data was not first processed with IFSF, then set the redshift
;        for plotting to be the systemic redshift. Otherwise, use the stellar 
;        redshift determined from the fit.
;         if ~ tag_exist(initnad,'noemlinfit') then zuse = struct.zstar $
;         else 
         zuse = initdat.zsys_gas
         if ~ noplot then begin
          print, 'Plotting...'
            if tag_exist(initnad,'argspltfitnad') then $
               ifsf_pltuvabsfit,galaxyname,(nadcube.wave)[initnad.plotindex[0]:initnad.plotindex[1]],$
                              (nadcube.dat)[initnad.plotindex[0]:initnad.plotindex[1]],$
                              (nadcube.cont)[initnad.plotindex[0]:initnad.plotindex[1]],$
                              (nadcube.flux)[initnad.plotindex[0]:initnad.plotindex[1]],$
                              param,doublet,directoryname,outfile+'_uvabs'+doublet+'_fit',zuse,$
                              _extra=initnad.argspltfitnad $
            else $
               ifsf_pltuvabsfit,galaxyname,(nadcube.wave)[initnad.plotindex[0]:initnad.plotindex[1]],$
                              (nadcube.dat)[initnad.plotindex[0]:initnad.plotindex[1]],$
                              (nadcube.cont)[initnad.plotindex[0]:initnad.plotindex[1]],$
                              (nadcube.flux)[initnad.plotindex[0]:initnad.plotindex[1]],$
                              param,doublet,directoryname,outfile+'_uvabs'+doublet+'_fit',zuse
         endif

;        Compute model equivalent widths
         weq=1
         nademflux=1
         IF (doublet eq 'OVI') THEN BEGIN
          modspec = ifsf_uvabsfcnOVI((nadcube.wave)[*],param,weq=weq,$
                               nademflux=nademflux)
         ENDIF
         IF (doublet eq 'NV') THEN BEGIN
          modspec = ifsf_uvabsfcnNV((nadcube.wave)[*],param,weq=weq,$
                               nademflux=nademflux)
         ENDIF

;;        Compute errors in fit
;         if ~ keyword_set(noerr) then begin
;            if dofirstemfit then nademfix_use = first_nademfix $
;            else nademfix_use = nademfix
;            if keyword_set(nomc) then plotonly=1b else plotonly=0b
;            errors = ifsf_fitnaderr([nhei,nnadabs,nnadem],(nadcube.wave)[*],$
;                                    modspec,(nadcube.err)[*],$
;                                    parinit,$
;                                    outfile+'_nad_errs.ps',outfile+'_nad_mc.xdr',$
;                                    dofirstemfit=dofirstemfit,$
;                                    first_parinit=first_parinit,$
;                                    first_modflux=first_modflux,$
;                                    heifix=heifix,nadabsfix=nadabsfix,$
;                                    nademfix=nademfix_use,$
;                                    niter=initnad.mcniter,$
;                                    nsplit=nsplit,quiet=quiet,weqerr=weqerr,$
;                                    nademfluxerr=nademfluxerr,noplot=noplot,$
;                                    plotonly=plotonly)
;         endif
;         
;         ifsf_printnadpar,nadparlun,i+1,j+1,param

;        Initialize cubes to hold physical quantities
         if firstfit then begin
            nadfit = $
;               Fit properties
               {chisq: dblarr(ncols,nrows)+bad,$
                chisq_emonly: dblarr(ncols,nrows)+bad,$
                dof: dblarr(ncols,nrows)+bad,$
                niter: dblarr(ncols,nrows)+bad,$
;               Equivalent widths and fluxes
                weqabs: dblarr(ncols,nrows,1+maxncomp)+bad,$
                weqabserr: dblarr(ncols,nrows,2)+bad,$
                weqem: dblarr(ncols,nrows,1+maxncomp)+bad,$
                weqemerr: dblarr(ncols,nrows,2)+bad,$
                totfluxem: dblarr(ncols,nrows,1+maxncomp)+bad,$
                totfluxemerr: dblarr(ncols,nrows,2)+bad,$
;               HeI line parameters
                wavehei: dblarr(ncols,nrows,maxncomp)+bad,$
                sigmahei: dblarr(ncols,nrows,maxncomp)+bad,$
                fluxhei: dblarr(ncols,nrows,maxncomp)+bad,$
;               UVAbs absorption line parameters
                cf: dblarr(ncols,nrows,maxncomp)+bad,$
                cferr: dblarr(ncols,nrows,maxncomp,2)+bad,$
                tau: dblarr(ncols,nrows,maxncomp)+bad,$
                tauerr: dblarr(ncols,nrows,maxncomp,2)+bad,$
                waveabs: dblarr(ncols,nrows,maxncomp)+bad,$
                waveabserr: dblarr(ncols,nrows,maxncomp,2)+bad,$
                sigmaabs: dblarr(ncols,nrows,maxncomp)+bad,$
                sigmaabserr: dblarr(ncols,nrows,maxncomp,2)+bad,$
;               UVAbs emission line parameters
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
;        Populate cubes of physical properties
         nadfit.chisq[i,j]=chisq
         if dofirstemfit then nadfit.chisq_emonly[i,j]=chisq_emonly
         nadfit.dof[i,j]=dof
         nadfit.niter[i,j]=niter
         nadfit.weqabs[i,j,0:nnadabs]=weq.abs
         nadfit.weqem[i,j,0:nnadem]=weq.em
         nadfit.totfluxem[i,j,0:nnadem]=nademflux
;         if ~ keyword_set(noerr) then begin
;            nadfit.weqabserr[i,j,*]=weqerr[0,*]
;            nadfit.weqemerr[i,j,*]=weqerr[1,*]
;            nadfit.totfluxemerr[i,j,*]=nademfluxerr
;         endif
;         if nhei gt 0 then begin
;            iarr = 3 + dindgen(nhei)*3
;            nadfit.wavehei[i,j,0:nhei-1]=param[iarr]
;            nadfit.sigmahei[i,j,0:nhei-1]=param[iarr+1]
;            nadfit.fluxhei[i,j,0:nhei-1]=param[iarr+2]
;         endif
         if nnadabs gt 0 then begin
            iarr = 3+ dindgen(nnadabs)*4
;            iarr = 3+nhei*3 + dindgen(nnadabs)*4
            nadfit.cf[i,j,0:nnadabs-1]=param[iarr]
            nadfit.tau[i,j,0:nnadabs-1]=param[iarr+1]
            nadfit.waveabs[i,j,0:nnadabs-1]=param[iarr+2]
            nadfit.sigmaabs[i,j,0:nnadabs-1]=param[iarr+3]
;            if ~ keyword_set(noerr) then begin
;               nadfit.cferr[i,j,0:nnadabs-1,*]=errors[iarr-3,*]
;               nadfit.tauerr[i,j,0:nnadabs-1,*]=errors[iarr-3+1,*]
;               nadfit.waveabserr[i,j,0:nnadabs-1,*]=errors[iarr-3+2,*]
;               nadfit.sigmaabserr[i,j,0:nnadabs-1,*]=errors[iarr-3+3,*]
;            endif
         endif
         if nnadem gt 0 then begin
            iarr = 3+nnadabs*4 + dindgen(nnadem)*4
;            iarr = 3+nhei*3+nnadabs*4 + dindgen(nnadem)*4
            nadfit.waveem[i,j,0:nnadem-1]=param[iarr]
            nadfit.sigmaem[i,j,0:nnadem-1]=param[iarr+1]
            nadfit.flux[i,j,0:nnadem-1]=param[iarr+2]
            nadfit.frat[i,j,0:nnadem-1]=param[iarr+3]
;            if ~ keyword_set(noerr) then begin
;               nadfit.waveemerr[i,j,0:nnadem-1,*]=errors[iarr-3,*]
;               nadfit.sigmaemerr[i,j,0:nnadem-1,*]=errors[iarr-3+1,*]
;               nadfit.fluxerr[i,j,0:nnadem-1,*]=errors[iarr-3+2,*]
;               nadfit.fraterr[i,j,0:nnadem-1,*]=errors[iarr-3+3,*]               
;            endif
         endif

nofit:

      endfor
      
   endfor

;  Velocity calculations
  zcomp=MAKE_ARRAY(2+4*comps)
  delz=MAKE_ARRAY(comps)
  velocity=MAKE_ARRAY(comps)
  FOR M = 0, comps-1 DO BEGIN
    zcomp[M] = param[4+4*M]/1242.804d - 1d
    delz[M]= zcomp[M]-redshift
    velocity[M] = c*((delz[M]+1)^2-1)/((delz[M]+1)^2+1)
  ENDFOR
    

finish:
  
   if ~ keyword_set(noxdr) then begin
      output=initdat.outdir+initdat.label+'_uvabs'+doublet+'.txt'
      openw, lun, output, /GET_LUN
      printf, lun,'Covering Factor','Optical Depth','Wavelength ($\Angstrom$)','Sigma','Velocity (km/s)', FORMAT='(A-20,A-20,A-20,A-20,A-20)'
      FOR M = 0, comps-1 DO BEGIN
        printf,lun, param[2+4*M:2+4*M+3],velocity[M], FORMAT='(F-20.4,F-20.4,F-20.4,F-20.4,F-20.4)'
      ENDFOR
      printf,lun,'============'
      printf, lun, 'Number of Absorption Components:',param[0], FORMAT='(A-30,I)'
      printf, lun, 'Number of Emission Components:',param[1], FORMAT='(A-30,I)'
      close, lun
      FREE_LUN, lun
   endif

;   free_lun,nadparlun
   print,'Total runtime: ',systime(1)-starttime,' s.',$
         format='(/,A0,I0,A0,/)'

end
