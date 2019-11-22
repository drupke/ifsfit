; docformat = 'rst'
;
;+
;
; Fit multiplets. Plot fit and write fit parameters to a file.
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
;    dir: in, required, type=string
;      Directory path where spectra data is located and where results 
;      are to be output.
;    galshort: in, required, type=string
;      Shorthand for galaxy to fit.
;    multiplet: in, required, type=string
;      Multiplet profiles to fit. Captilization is required.
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
;      2019nov20, DSNR, copied from IFSF_FITDOUBLET
;    
; :Copyright:
;    Copyright (C) 2019 David S. N. Rupke
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
pro ifsf_fitmultiplet,dir,galshort,multiplet,fcngalinfo,$
                      cols=cols,rows=rows,verbose=verbose,$
                      noxdr=noxdr,noplot=noplot,weights=weights,$
                      noerr=noerr,nomc=nomc,init=init,$
                      argsgalinfo=argsgalinfo

   bad = 1d99
   c = 299792.458d

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
   linelist = ifsf_linelist([initstr.reflinename,initstr.linenames])

   if tag_exist(initstr,'taumax') then taumax=initstr.taumax $
   else taumax = 5d

   ext='.txt'

   ncols = 1
   nrows = 1
   oned=1
   nz = 1
   cube = {err: initstr.error, $
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
            if oned then lab = string(i+1,format='(I04)') $
            else lab = string(i+1,'_',j+1,format='(I04,A,I04)')
            outfile = initstr.outdir+initstr.galaxy
         endif

;        Get UVAbs absorption parameters
         nabs = initstr.nabs[i,j]
         if nabs gt 0 then begin
            if tag_exist(initstr,'abs_cfinit') then BEGIN
               cfinit=initstr.abs_cfinit
             ENDIF else cfinit = dblarr(nabs)+0.5d
             if tag_exist(initstr,'abs_tauinit') then BEGIN
               tauinit=initstr.abs_tauinit
             ENDIF else tauinit = dblarr(nabs)+0.5d
             winit = reform(((initstr.abs_zinit)[i,j,0:nabs-1]$
               +1d)*linelist[initstr.reflinename],nabs)
             siginit = reform((initstr.abs_siginit)$
               [i,j,0:nabs-1],nabs)
             initabs = [[cfinit],[tauinit],[winit],[siginit]]
         endif else initabs=0

;        Fill out parameter structure with initial guesses and constraints

         if nabs gt 0 AND tag_exist(initstr,'abs_fix') then $
            absfix=reform((initstr.abs_fix)[i,j,0:nabs-1,*],nabs,4) $
         else absfix=0b

         parinit = $
            call_function(initstr.fcninitpar,initabs,$
                          initstr.abs_siglim,initstr.reflinename,$
                          absfix=absfix,taumax=taumax)

;        Plot initial guess if requested
         if keyword_set(init) then begin
            zuse = initstr.zsys_gas
            ilo = initstr.plotindex[0]
            ihi = initstr.plotindex[1]
            if tag_exist(initstr,'argspltmultiplet') then $
               ifsf_pltmultiplet,galshort,(cube.wave)[ilo:ihi],$
                                 (cube.dat)[ilo:ihi],$
                                 (cube.cont)[ilo:ihi],$
                                 (cube.flux)[ilo:ihi],$
                                 parinit.value,linelist,initstr,multiplet,dir,$
                                 outfile+multiplet+'_init',zuse,/init,$
                                 _extra=initstr.argspltmultiplet $
            else $
               ifsf_pltmultiplet,galshort,(cube.wave)[ilo:ihi],$
                                 (cube.dat)[ilo:ihi],$
                                 (cube.cont)[ilo:ihi],$
                                 (cube.flux)[ilo:ihi],$
                                 parinit.value,linelist,initstr,multiplet,dir,$
                                 outfile+multiplet+'_init',zuse,/init
                       
            goto,fullstop
         endif

         if (~ keyword_set(weights)) then begin
           param = Mpfitfun(initstr.fcnfitmultiplet,$
                            (cube.wave)[initstr.fitindex[0]:initstr.fitindex[1]],$
                            (cube.dat)[initstr.fitindex[0]:initstr.fitindex[1]],$
                            (cube.err)[initstr.fitindex[0]:initstr.fitindex[1]],$
                            parinfo=parinit,perror=perror,maxiter=500,$
                            bestnorm=chisq,covar=covar,yfit=specfit,dof=dof,$
                            nfev=nfev,niter=niter,status=status,quiet=quiet,$
                            npegged=npegged,ftol=1D-6,errmsg=errmsg,$
                            functargs={refloglf:initstr.refloglf,$
                                       refmultwave:initstr.refmultwave,$
                                       loglf:initstr.loglf,$
                                       multwave:initstr.multwave})
         endif else begin
           param = Mpfitfun(initstr.fcnfitmultiplet,$
                            (cube.wave)[initstr.fitindex[0]:initstr.fitindex[1]],$
                            (cube.dat)[initstr.fitindex[0]:initstr.fitindex[1]],$
                            (cube.err)[initstr.fitindex[0]:initstr.fitindex[1]],$
                            parinfo=parinit,perror=perror,maxiter=500,$
                            bestnorm=chisq,covar=covar,yfit=specfit,dof=dof,$
                            nfev=nfev,niter=niter,status=status,quiet=quiet,$
                            npegged=npegged,ftol=1D-6,errmsg=errmsg,$
                            functargs={refloglf:initstr.refloglf,$
                                       refmultwave:initstr.refmultwave,$
                                       loglf:initstr.loglf,$
                                       multwave:initstr.multwave},/WEIGHTS)          
         endelse
         if status eq 5 then print,'IFSF_FITMULTIPLET: Max. iterations reached.'
         if status eq 0 OR status eq -16 then begin
            print,'IFSF_FITMULTIPLET: Error in MPFIT. Aborting.'
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
            if tag_exist(initstr,'argspltmultiplet') then $
               ifsf_pltmultiplet,galshort,(cube.wave)[initstr.plotindex[0]:initstr.plotindex[1]],$
                                (cube.dat)[initstr.plotindex[0]:initstr.plotindex[1]],$
                                (cube.cont)[initstr.plotindex[0]:initstr.plotindex[1]],$
                                (cube.flux)[initstr.plotindex[0]:initstr.plotindex[1]],$
                                param,linelist,initstr,multiplet,dir,outfile+multiplet+'_fit',zuse,$
                                _extra=initstr.argspltmultiplet $
            else $
               ifsf_pltmultiplet,galshort,(cube.wave)[initstr.plotindex[0]:initstr.plotindex[1]],$
                                (cube.dat)[initstr.plotindex[0]:initstr.plotindex[1]],$
                                (cube.cont)[initstr.plotindex[0]:initstr.plotindex[1]],$
                                (cube.flux)[initstr.plotindex[0]:initstr.plotindex[1]],$
                                param,linelist,initstr,multiplet,dir,outfile+multiplet+'_fit',zuse
         endif

;        Compute model equivalent widths
         weq=1d
         modspec = ifsf_multipletfcn(cube.wave,param,refloglf=initstr.refloglf,$
                                     refmultwave=initstr.refmultwave,$
                                     loglf=initstr.loglf,multwave=initstr.multwave,$
                                     weq=weq)


;        Initialize cubes to hold physical quantities
         if firstfit then begin
            fit = $
;               Fit properties
               {chisq: dblarr(ncols,nrows)+bad,$
                dof: dblarr(ncols,nrows)+bad,$
                niter: dblarr(ncols,nrows)+bad,$
;               Equivalent widths and fluxes
                weqabs: dblarr(ncols,nrows,1+maxncomp)+bad,$
                weqabserr: dblarr(ncols,nrows,2)+bad,$
;               multiplet absorption line parameters
                cf: dblarr(ncols,nrows,maxncomp)+bad,$
                cferr: dblarr(ncols,nrows,maxncomp,2)+bad,$
                tau: dblarr(ncols,nrows,maxncomp)+bad,$
                tauerr: dblarr(ncols,nrows,maxncomp,2)+bad,$
                waveabs: dblarr(ncols,nrows,maxncomp)+bad,$
                waveabserr: dblarr(ncols,nrows,maxncomp,2)+bad,$
                sigmaabs: dblarr(ncols,nrows,maxncomp)+bad,$
                sigmaabserr: dblarr(ncols,nrows,maxncomp,2)+bad}
            firstfit = 0
         endif
;        Populate cubes of physical properties
         fit.chisq[i,j]=chisq
         fit.dof[i,j]=dof
         fit.niter[i,j]=niter
         fit.weqabs[i,j,0:nabs]=weq.abs
         if nabs gt 0 then begin
            iarr = 1 + dindgen(nabs)*4
            fit.cf[i,j,0:nabs-1]=param[iarr]
            fit.tau[i,j,0:nabs-1]=param[iarr+1]
            fit.waveabs[i,j,0:nabs-1]=param[iarr+2]
            fit.sigmaabs[i,j,0:nabs-1]=param[iarr+3]
         endif

nofit:

      endfor
      
   endfor

;  Velocity calculations
   comps=maxncomp
   delz=MAKE_ARRAY(comps)
   velocity=MAKE_ARRAY(comps)
   FOR M = 0, comps-1 DO BEGIN
      delz[M] = param[3+4*M]/(linelist[initstr.reflinename]*(1d + redshift)) - 1d
      velocity[M] = c*((delz[M]+1d)^2d -1d)/((delz[M]+1d)^2d +1d)
   ENDFOR

finish:
  
   if ~ keyword_set(noxdr) then begin
      output=initstr.outdir+initstr.galaxy+multiplet+'par_best.txt'
      openw, lun, output, /GET_LUN
      printf, lun, param[0], '[Number of Absorption Components]',FORMAT='(I-3,A0)'
      printf, lun, weq.abs[0], ' [Total equivalent width in A]',FORMAT='(D0.4,A0)'
      printf, lun,'Covering Factor','Optical Depth',$
                  'Wavelength(A)','Sigma(km/s)','Velocity(km/s)',$
                  FORMAT='(5A-20)'
      FOR M = 0, maxncomp-1 DO $
         printf,lun, param[1+4*M:1+4*M+3],velocity[M],FORMAT='(5F-20.4)'
      FREE_LUN, lun
   endif

fullstop:

   print,'Total runtime: ',systime(1)-starttime,' s.',$
         format='(/,A0,I0,A0,/)'

end
