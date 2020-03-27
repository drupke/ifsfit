; docformat = 'rst'
;
;+
;
; Fit O VI and N V doublets. Plot fit and write fit parameters to a file.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Produces a plot ([gal_col_row]_doublet_fit.jpg) for each fit, and a
;    parameter file ([gal].doublet.dat) containing the fit parameters for
;    all spaxels.
;
; :Params:
;    table: in, required, type=string
;      Table containing input redshifts.
;    dir: in, required, type=string
;      Directory path where spectra data is located and where results 
;      are to be output.
;    galshort: in, required, type=string
;      Shorthand for galaxy to fit.
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
;    init: in, optional, type=byte
;      Plot initial guess rather than fitting.
;    noerr: in, optional, type=byte
;      Do not Monte Carlo the errors.
;      [No error computation done, presently.]
;    nomc: in, optional, type=byte
;      Do not re-run Monte Carlo error simulations, but do grab old outputs
;      and append them to the output structure.
;      [No error computation done, presently.]
;    noplot: in, optional, type=byte
;      Do not produce plots.
;    noxdr: in, optional, type=byte
;      Do no write XDR file with fit parameters.
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
;      2015jul01, AT, began work on adapting from ifs_fitdoublet
;      2016aug08, DSNR, added option to examine initial guess
;      2018jun05, DSNR, added measurements of depth-weighted velocity and
;                       RMS
;    
; :Copyright:
;    Copyright (C) 2015--2018 Anthony To, David S. N. Rupke
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
pro ifsf_fitdoublet,dir,galshort,doublet,fcngalinfo,$
                    cols=cols,rows=rows,verbose=verbose,$
                    noxdr=noxdr,noplot=noplot,weights=weights,$
                    noerr=noerr,nomc=nomc,init=init,$
                    argsgalinfo=argsgalinfo

   bad = 1d99
   doublet_emrat_init = 1.5d
   c = 299792.458d

   IF (doublet eq 'MgII') THEN linename = 'MgII2802'
   IF (doublet eq 'NV') THEN linename = 'NV1242'
   IF (doublet eq 'OVI') THEN linename = 'OVI1037'
   IF (doublet eq 'PV') THEN linename = 'PV1128'
   IF (doublet eq 'FeIIUV1') THEN linename = 'FeII2585'
   IF (doublet eq 'FeIIUV2') THEN linename = 'FeII2373'

   starttime = systime(1)
   time = 0
   if keyword_set(verbose) then quiet=0 else quiet=1
   if ~ keyword_set(noplot) then noplot=0 else noplot=1

   ; Get initial fit parameters and galaxy properties
   redshift = 0d
   if keyword_set(argsgalinfo) then begin
      initstr = call_function(fcngalinfo,redshift,argsgalinfo=argsgalinfo)
   endif else begin
      initstr = call_function(fcngalinfo,redshift)
   endelse

   maxncomp = initstr.maxncomp
   
   ; Get linelist
   linelist = $
      ifsf_linelist(['OVI1031','OVI1037','Lyalpha','Lybeta',$
                     'NV1238','NV1242','PV1117','PV1128','MgII2795','MgII2802',$
                     'FeII2585','FeII2599','FeII2373','FeII2382'])

   if tag_exist(initstr,'taumax') then taumax=initstr.taumax $
   else taumax = 5d

   ext='.txt'
   
   if tag_exist(initstr,'doubletem_siglim') then $
      doubletem_siglim = initstr.doubletem_siglim $
   else doubletem_siglim = 0d

   ncols = 1
   nrows = 1
   oned=1
   nz = 1
   doubletcube = {err: initstr.error, $
                  dat : initstr.relativeflux, $
                  wave : initstr.wavelength, $ 
                  cont: initstr.continuum, $
                  flux: initstr.flux $   
                 }

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

         if ~ tag_exist(initstr,'noemlinfit') then begin
;           Restore continuum + emission-line fit
            if oned then lab = string(i+1,format='(I04)') $
            else lab = string(i+1,'_',j+1,format='(I04,A,I04)')
;            infile = initstr.outdir+initstr.galaxy+ext
;            print, infile
            outfile = initstr.outdir+initstr.galaxy
;            if ~ file_test(infile) then begin
;               print,'IFSF_FITDOUBLET: No file for ',i+1,', ',j+1,'.',$
;                     format='(A0,I4,A0,I4,A0)'
;               goto,nofit
;            endif
;            restore,file=infile
         endif

;        Get UVAbs absorption parameters
         ndoubletabs = initstr.ndoubletabs[i,j]
         if ndoubletabs gt 0 then begin
            if tag_exist(initstr,'doubletabs_cfinit') then BEGIN
               cfinit=initstr.doubletabs_cfinit
             ENDIF else cfinit = dblarr(ndoubletabs)+0.5d
             if tag_exist(initstr,'doubletabs_tauinit') then BEGIN
               tauinit=initstr.doubletabs_tauinit
             ENDIF else tauinit = dblarr(ndoubletabs)+0.5d
             winit = reform(((initstr.doubletabs_zinit)[i,j,0:ndoubletabs-1]$
               +1d)*linelist[linename],ndoubletabs)
             siginit = reform((initstr.doubletabs_siginit)$
               [i,j,0:ndoubletabs-1],ndoubletabs)
             initabs = [[cfinit],[tauinit],[winit],[siginit]]
         endif else initabs=0
         
;        Get UVAbs emission parameters
         ndoubletem = initstr.ndoubletem[i,j]
;        placeholders for case of separately fitting emission and absorption
         dofirstemfit=0b
         first_doubletemfix=0d
         first_parinit=0d
         first_modflux=0d
         if ndoubletem gt 0 then begin
               
            if doubletem_siglim[0] eq 0 then print,'IFSF_FITUVABS: ERROR: Emission '+$
               'line sigma limits not set (UVABSEM_SIGLIM) in INITUVABS '+$
               'structure, but emission lines need to be fit.'
               
            winit = reform(((initstr.doubletem_zinit)[i,j,0:ndoubletem-1]+1d)$
                           *linelist[linename],ndoubletem)
            siginit = reform((initstr.doubletem_siginit)[i,j,0:ndoubletem-1],ndoubletem)
            if tag_exist(initstr,'doubletem_finit') then $
               finit = reform((initstr.doubletem_finit)[i,j,0:ndoubletem-1],ndoubletem) $
            else finit = dblarr(ndoubletem)+0.1d
            if tag_exist(initstr,'doubletem_rinit') then $
               rinit = reform((initstr.doubletem_rinit)[i,j,0:ndoubletem-1],ndoubletem) $
            else rinit = dblarr(ndoubletem)+doublet_emrat_init
            initem = [[winit],[siginit],[finit],[rinit]]
            if tag_exist(initstr,'doubletem_fix') then $
               doubletemfix=reform((initstr.doubletem_fix)[i,j,0:ndoubletem-1,*],ndoubletem,4) $
            else doubletemfix=0b
         endif else begin
            initem=0
            doubletemfix=0b
         endelse

;        Fill out parameter structure with initial guesses and constraints

         if ndoubletabs gt 0 AND tag_exist(initstr,'doubletabs_fix') then $
            doubletabsfix=reform((initstr.doubletabs_fix)[i,j,0:ndoubletabs-1,*],ndoubletabs,4) $
         else doubletabsfix=0b

         parinit = $
            call_function(initstr.fcninitpar,doublet,initabs,initem,$
                          initstr.doubletabs_siglim,doubletem_siglim,$
                          doubletabsfix=doubletabsfix,doubletemfix=doubletemfix,$
                          taumax=taumax)

;        Plot initial guess if requested
         if keyword_set(init) then begin
            zuse = initstr.zsys_gas
            ilo = initstr.plotindex[0]
            ihi = initstr.plotindex[1]
            ifsf_pltdoublet,galshort,(doubletcube.wave)[ilo:ihi],$
                             (doubletcube.dat)[ilo:ihi],$
                             (doubletcube.cont)[ilo:ihi],$
                             (doubletcube.flux)[ilo:ihi],$
                             parinit.value,doublet,dir,$
                             outfile+doublet+'_init',zuse,/init
            goto,fullstop
         endif

         if (~ keyword_set(weights)) then begin
           param = Mpfitfun(initstr.fcnfitdoublet,$
                            (doubletcube.wave)[initstr.fitindex[0]:initstr.fitindex[1]],$
                            (doubletcube.dat)[initstr.fitindex[0]:initstr.fitindex[1]],$
                            (doubletcube.err)[initstr.fitindex[0]:initstr.fitindex[1]],$
                            parinfo=parinit,perror=perror,maxiter=500,$
                            bestnorm=chisq,covar=covar,yfit=specfit,dof=dof,$
                            nfev=nfev,niter=niter,status=status,quiet=quiet,$
                            npegged=npegged,ftol=1D-6,errmsg=errmsg,$
                            functargs={doubletname:doublet})
         endif else begin
           param = Mpfitfun(initstr.fcnfitdoublet,$
                            (doubletcube.wave)[initstr.fitindex[0]:initstr.fitindex[1]],$
                            (doubletcube.dat)[initstr.fitindex[0]:initstr.fitindex[1]],$
                            (doubletcube.err)[initstr.fitindex[0]:initstr.fitindex[1]],$
                            parinfo=parinit,perror=perror,maxiter=500,$
                            bestnorm=chisq,covar=covar,yfit=specfit,dof=dof,$
                            nfev=nfev,niter=niter,status=status,quiet=quiet,$
                            npegged=npegged,ftol=1D-6,errmsg=errmsg,$
                            functargs={doubletname:doublet},/WEIGHTS)          
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
;         if ~ tag_exist(initstr,'noemlinfit') then zuse = struct.zstar $
;         else 
         zuse = initstr.zsys_gas
         if ~ noplot then begin
          print, 'Plotting...'
            if tag_exist(initstr,'argspltfitdoublet') then $
               ifsf_pltdoublet,galshort,(doubletcube.wave)[initstr.plotindex[0]:initstr.plotindex[1]],$
                                (doubletcube.dat)[initstr.plotindex[0]:initstr.plotindex[1]],$
                                (doubletcube.cont)[initstr.plotindex[0]:initstr.plotindex[1]],$
                                (doubletcube.flux)[initstr.plotindex[0]:initstr.plotindex[1]],$
                                param,doublet,dir,outfile+doublet+'_fit',zuse,$
                                _extra=initstr.argspltfitdoublet $
            else $
               ifsf_pltdoublet,galshort,(doubletcube.wave)[initstr.plotindex[0]:initstr.plotindex[1]],$
                                (doubletcube.dat)[initstr.plotindex[0]:initstr.plotindex[1]],$
                                (doubletcube.cont)[initstr.plotindex[0]:initstr.plotindex[1]],$
                                (doubletcube.flux)[initstr.plotindex[0]:initstr.plotindex[1]],$
                                param,doublet,dir,outfile+doublet+'_fit',zuse
         endif

;        Compute model equivalent widths
         weq=1
         doubletemflux=1
         vwtabs=1
;        Velocity array
         delz = doubletcube.wave/(linelist[linename]*(1d + redshift)) - 1d
         veltmp = c*((delz+1d)^2d -1d)/((delz+1d)^2d +1d)
         modspec = ifsf_doubletfcn(doubletcube.wave,param,doubletname=doublet,$
                                   weq=weq,doubletemflux=doubletemflux,$
                                   vels=veltmp,vwtabs=vwtabs)

;;        Compute errors in fit
;         if ~ keyword_set(noerr) then begin
;            if dofirstemfit then doubletemfix_use = first_doubletemfix $
;            else doubletemfix_use = doubletemfix
;            if keyword_set(nomc) then plotonly=1b else plotonly=0b
;            errors = ifsf_fitdoubleterr([nhei,ndoubletabs,ndoubletem],(doubletcube.wave)[*],$
;                                    modspec,(doubletcube.err)[*],$
;                                    parinit,$
;                                    outfile+'_doublet_errs.ps',outfile+'_doublet_mc.xdr',$
;                                    dofirstemfit=dofirstemfit,$
;                                    first_parinit=first_parinit,$
;                                    first_modflux=first_modflux,$
;                                    heifix=heifix,doubletabsfix=doubletabsfix,$
;                                    doubletemfix=doubletemfix_use,$
;                                    niter=initstr.mcniter,$
;                                    nsplit=nsplit,quiet=quiet,weqerr=weqerr,$
;                                    doubletemfluxerr=doubletemfluxerr,noplot=noplot,$
;                                    plotonly=plotonly)
;         endif
;         
;         ifsf_printdoubletpar,doubletparlun,i+1,j+1,param

;        Initialize cubes to hold physical quantities
         if firstfit then begin
            doubletfit = $
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
;               doublet absorption line parameters
                cf: dblarr(ncols,nrows,maxncomp)+bad,$
                cferr: dblarr(ncols,nrows,maxncomp,2)+bad,$
                tau: dblarr(ncols,nrows,maxncomp)+bad,$
                tauerr: dblarr(ncols,nrows,maxncomp,2)+bad,$
                waveabs: dblarr(ncols,nrows,maxncomp)+bad,$
                waveabserr: dblarr(ncols,nrows,maxncomp,2)+bad,$
                sigmaabs: dblarr(ncols,nrows,maxncomp)+bad,$
                sigmaabserr: dblarr(ncols,nrows,maxncomp,2)+bad,$
;               doublet emission line parameters
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
         doubletfit.chisq[i,j]=chisq
         if dofirstemfit then doubletfit.chisq_emonly[i,j]=chisq_emonly
         doubletfit.dof[i,j]=dof
         doubletfit.niter[i,j]=niter
         doubletfit.weqabs[i,j,0:ndoubletabs]=weq.abs
         doubletfit.weqem[i,j,0:ndoubletem]=weq.em
         doubletfit.totfluxem[i,j,0:ndoubletem]=doubletemflux
;         if ~ keyword_set(noerr) then begin
;            doubletfit.weqabserr[i,j,*]=weqerr[0,*]
;            doubletfit.weqemerr[i,j,*]=weqerr[1,*]
;            doubletfit.totfluxemerr[i,j,*]=doubletemfluxerr
;         endif
         if ndoubletabs gt 0 then begin
            iarr = 3+ dindgen(ndoubletabs)*4
            doubletfit.cf[i,j,0:ndoubletabs-1]=param[iarr]
            doubletfit.tau[i,j,0:ndoubletabs-1]=param[iarr+1]
            doubletfit.waveabs[i,j,0:ndoubletabs-1]=param[iarr+2]
            doubletfit.sigmaabs[i,j,0:ndoubletabs-1]=param[iarr+3]
;            if ~ keyword_set(noerr) then begin
;               doubletfit.cferr[i,j,0:ndoubletabs-1,*]=errors[iarr-3,*]
;               doubletfit.tauerr[i,j,0:ndoubletabs-1,*]=errors[iarr-3+1,*]
;               doubletfit.waveabserr[i,j,0:ndoubletabs-1,*]=errors[iarr-3+2,*]
;               doubletfit.sigmaabserr[i,j,0:ndoubletabs-1,*]=errors[iarr-3+3,*]
;            endif
         endif
         if ndoubletem gt 0 then begin
            iarr = 3+ndoubletabs*4 + dindgen(ndoubletem)*4
            doubletfit.waveem[i,j,0:ndoubletem-1]=param[iarr]
            doubletfit.sigmaem[i,j,0:ndoubletem-1]=param[iarr+1]
            doubletfit.flux[i,j,0:ndoubletem-1]=param[iarr+2]
            doubletfit.frat[i,j,0:ndoubletem-1]=param[iarr+3]
;            if ~ keyword_set(noerr) then begin
;               doubletfit.waveemerr[i,j,0:ndoubletem-1,*]=errors[iarr-3,*]
;               doubletfit.sigmaemerr[i,j,0:ndoubletem-1,*]=errors[iarr-3+1,*]
;               doubletfit.fluxerr[i,j,0:ndoubletem-1,*]=errors[iarr-3+2,*]
;               doubletfit.fraterr[i,j,0:ndoubletem-1,*]=errors[iarr-3+3,*]               
;            endif
         endif

nofit:

      endfor
      
   endfor

;  Velocity calculations
   comps=maxncomp
   delz=MAKE_ARRAY(comps)
   velocity=MAKE_ARRAY(comps)
   FOR M = 0, comps-1 DO BEGIN
      delz[M] = param[4+4*M]/(linelist[linename]*(1d + redshift)) - 1d
      velocity[M] = c*((delz[M]+1d)^2d -1d)/((delz[M]+1d)^2d +1d)
   ENDFOR


;  Velocity calculations
   if param[1] gt 0 then begin
      comps_em=param[1]
      delz=MAKE_ARRAY(comps)
      velocity_em=MAKE_ARRAY(comps)
      FOR M = 0, comps_em-1 DO BEGIN
        delz[M] = param[2+param[0]*4+4*M]/(linelist[linename]*(1d + redshift)) - 1d
        velocity_em[M] = c*((delz[M]+1d)^2d -1d)/((delz[M]+1d)^2d +1d)
      ENDFOR
   endif
    

finish:
  
   lineofdashes = strjoin(replicate('-',62))

   if ~ keyword_set(noxdr) then begin
      output=initstr.outdir+initstr.galaxy+doublet+'par_best.txt'
      openw, lun, output, /GET_LUN
      printf, lun, lineofdashes
      printf, lun, 'Absorption',FORMAT='(A0)'
      printf, lun, lineofdashes
      printf, lun, param[0], '[Number of Components]',FORMAT='(I-3,A0)'
      printf, lun, weq.abs[0], ' [Total equivalent width in A]',FORMAT='(D0.4,A0)'
      printf, lun, vwtabs[0], ' [Weighted avg. vel. in km/s]',FORMAT='(D0.2,A0)'
      printf, lun, vwtabs[1], ' [Weighted RMS vel. in km/s]',FORMAT='(D0.2,A0)'
      printf, lun, lineofdashes
      printf, lun,'Cov. Factor','Tau',$
              'Wave(A)','Sigma(km/s)','Vel(km/s)',$
              FORMAT='(A-12,A-12,A-12,A-12,A-12)'
      printf, lun, lineofdashes
      FOR M = 0, comps-1 DO BEGIN
         printf,lun, param[2+4*M:2+4*M+3],velocity[M],$
                FORMAT='(F-12.4,F-12.4,F-12.4,F-12.4,F-12.4)'
      ENDFOR
      if param[1] gt 0 then begin
         printf, lun, lineofdashes
         printf, lun, 'Emission',FORMAT='(A0)'
         printf, lun, lineofdashes
         printf, lun, param[1], '[Number of Components]',FORMAT='(I-3,A0)'
         printf, lun, lineofdashes
         printf, lun,'Flux','2796/2803',$
                'Wave(A)','Sigma(km/s)','Vel(km/s)',$
                FORMAT='(A-12,A-12,A-12,A-12,A-12)'
         printf, lun, lineofdashes
         FOR M = 0, comps_em-1 DO BEGIN
            printf,lun, param[2+4*param[0]+4*M+2:2+4*param[0]+4*M+3],$
                   param[2+4*param[0]+4*M:2+4*param[0]+4*M+1],$
                   velocity_em[M],$
                   FORMAT='(F-12.4,F-12.4,F-12.4,F-12.4,F-12.4)'
         ENDFOR
      endif
      close, lun
      FREE_LUN, lun
   endif

fullstop:

;   free_lun,doubletparlun
   print,'Total runtime: ',systime(1)-starttime,' s.',$
         format='(/,A0,I0,A0,/)'

end
