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
;      2015jun03, DSNR, adjusted treatment of emission line sigma limits
;      2016nov03, DSNR, added convolution with spectral resolution
;    
; :Copyright:
;    Copyright (C) 2013--2016 David S. N. Rupke
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
   linelist = ifsf_linelist(['NaD1','NaD2','HeI5876'],/quiet)

   if tag_exist(initnad,'taumax') then taumax=initnad.taumax $
   else taumax = 5d

   if tag_exist(initnad,'nadem_siglim') then $
      nadem_siglim = initnad.nadem_siglim $
   else nadem_siglim = 0d

;  NADFIT tags
   tags2d = ['CHISQ','CHISQ_EMONLY','DOF','NITER']
   tags3d = ['WEQABS','WEQABSERR','WEQEM','WEQEMERR','TOTFLUXEM','TOTFLUXEMERR',$
             'WAVEHEI','SIGMAHEI','FLUXHEI','CF','TAU','WAVEABS','SIGMAABS',$
             'WAVEEM','SIGMAEM','FLUX','FRAT']
   tags4d = ['CFERR','TAUERR','WAVEABSERR','SIGMAABSERR',$
             'WAVEEMERR','SIGMAEMERR','FLUXERR','FRATERR']

;  Estimated spectral resolution for GMOS, B600 grating based on measurements.
;  Website says R = 1688 at 4610 A for 0.5" slit, with IFU 0.31" eff. slit.
;  This gives 1.69 A FWHM. I measure sometimes closer to 1.5-1.6.
;  Sigma is then in the range 0.64 -- 0.72. Use the former for flexibility.
   if ~ tag_exist(initnad,'specres') then specres = 0.64d $
   else specres = initnad.specres
   argsfitnad = {specres: specres}

   ifsf_printnadpar,nadparlun,outfile=initdat.outdir+initdat.label+'.nad.dat'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if ~ tag_exist(initdat,'noemlinfit') then $
      restore,file=initdat.outdir+initdat.label+'.lin.xdr'
   restore,file=initdat.outdir+initdat.label+'.nadspec.xdr'

   nadsize = size(nadcube.wave)
   ncols = nadsize[1]
   nrows = nadsize[2]
   if nrows eq 1 then oned=1b else oned=0b
   nz = nadsize[3]
   
   if tag_exist(initdat,'vormap') then begin
     vormap=initdat.vormap
     nvorcols = max(vormap)
     vordone = bytarr(nvorcols)
     vorrefnad = intarr(nvorcols,2)
     vorreflines = intarr(nvorcols,2)
     for i=1,nvorcols do begin
        ivor = where(vormap eq i,ctivor)
        xyvor = array_indices(vormap,ivor[0])
        vorreflines[i-1,*] = xyvor
     endfor

   endif

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

;        Check data quality
         igd = where(nadcube.weq[i,j,*] ne bad,ctgd)
         if ctgd eq 0 then begin
            message,'Skipping this spectrum b/c <SNR> < threshold (IFSF_NORMNAD).',/cont
            goto,nofit
         endif

;;        Check that line was empirically detected
;         igd = where(nadcube.iweq[i,j,*] ne -1,ctgd)
;         if ctgd eq 0 then begin
;            message,'Skipping this spectrum b/c no empirical detection.',/cont
;            goto,nofit
;         endif

         didvor = 0b
         if tag_exist(initdat,'vormap') then begin
           if finite(initdat.vormap[i,j]) AND $
              initdat.vormap[i,j] ne bad then begin
              if vordone[vormap[i,j]-1] then didvor = 1b
              ireflinefit = vorreflines[vormap[i,j]-1,0]
              jreflinefit = vorreflines[vormap[i,j]-1,1]
           endif else begin
              goto,nofit
           endelse
         endif else begin
            ireflinefit = i
            jreflinefit = j
         endelse
         if didvor then goto,copyvor

         if ~ tag_exist(initnad,'noemlinfit') then begin
;           Restore continuum + emission-line fit
            if oned then begin
               inlab = string(i+1,format='(I04)')
               outlab = inlab
            endif else begin
               inlab = string(ireflinefit+1,'_',jreflinefit+1,format='(I04,A,I04)')
               outlab = string(i+1,'_',j+1,format='(I04,A,I04)')
            endelse
            infile = initdat.outdir+initdat.label+'_'+inlab+'.xdr'
            outfile = initdat.outdir+initdat.label+'_'+outlab
            if ~ file_test(infile) then begin
               message,'No XDR file for '+string(ireflinefit+1,format='(I0)')+', '+$
                       string(jreflinefit+1,format='(I0)')+'.',/cont
               goto,nofit
            endif
            restore,file=infile
         endif

;        Get fit range
         if tag_exist(initnad,'nadfitran') then begin
            ilamlo = value_locate((nadcube.wave)[i,j,*],initnad.nadfitran[0])
            ilamhi = value_locate((nadcube.wave)[i,j,*],initnad.nadfitran[1])
         endif else begin
            ilamlo = 0
            ilamhi = n_elements((nadcube.wave)[i,j,*])-1
         endelse
         passwav = reform((nadcube.wave)[i,j,ilamlo:ilamhi],ilamhi-ilamlo+1)
         passdat = reform((nadcube.dat)[i,j,ilamlo:ilamhi],ilamhi-ilamlo+1)
         passcont = reform((nadcube.cont)[i,j,ilamlo:ilamhi],ilamhi-ilamlo+1)
         passerr = reform((nadcube.err)[i,j,ilamlo:ilamhi],ilamhi-ilamlo+1)

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
                     message,'No XDR file for ',i+1,', ',j+1,'.',$
                             format='(A0,I4,A0,I4,A0)',/cont
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
               if nhei gt 0 then begin
                  inithei = [[(linepars.wave)[refline,0:nhei-1]/$
                              reflinelist[refline]*linelist['HeI5876']],$
                             [(linepars.sigma)[refline,0:nhei-1]],$
                             [dblarr(nhei)+0.1d]]
                  heifix = bytarr(nhei,3)
                  heifix[*,0:1] = 1b              
               endif else begin
                  inithei=0
                  heifix=0
               endelse
;           Case of allowing the line to vary freely
            endif else if tag_exist(initnad,'hei_zinit') AND $
                          tag_exist(initnad,'hei_siginit') AND $
                          initnad.nhei[i,j] gt 0 then begin
               nhei=initnad.nhei[i,j]
; Bugfix -- 7/13/15 -- thanks Elise Hampton
               inithei = [[(((initnad.hei_zinit)[i,j,0:nhei-1])[0:nhei-1]+1d)*$
                           linelist['HeI5876']],$
                          [((initnad.hei_siginit)[i,j,0:nhei-1])[0:nhei-1]],$
                          [dblarr(nhei)+0.1d]]
;               inithei = [[((initnad.hei_zinit)[i,j,0:nhei-1]+1d)*$
;                           linelist['HeI5876']],$
;                          [(initnad.hei_siginit)[i,j,0:nhei-1]],$
;                          [dblarr(nhei)+0.1d]]
               if tag_exist(initnad,'hei_fix') then $
                  heifix=reform((initnad.hei_fix)[i,j,0:nhei-1,*],nhei,3) $
               else heifix=bytarr(nhei,3)
            endif else begin
               message,'HeI5876 initialization parameters'+$
                       ' not properly specified.',/cont
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
               
            if nadem_siglim[0] eq 0 then message,'Emission '+$
               'line sigma limits not set (NADEM_SIGLIM) in INITNAD '+$
               'structure, but emission lines need to be fit.',/cont
               
            winit = reform(((initnad.nadem_zinit)[i,j,0:nnadem-1]+1d)$
                           *linelist['NaD1'],nnadem)
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
            
;           Fit the emission line region only if requested, and if NaD emission 
;           and absorption (or NaD emission and HeI emission) are to be fit.
            if tag_exist(initnad,'nadem_fitinit') AND $
               (nnadabs gt 0 OR nhei gt 0) then begin

               dofirstemfit=1b

               refit=1b
               while refit ne 0b do begin
               
;              Fill out parameter structure with initial guesses and constraints
               if tag_exist(initnad,'argsinitpar') then parinit = $
                  call_function(initnad.fcninitpar,0,0,initnadem,$
                                initnad.nadabs_siglim,nadem_siglim,$
                                heifix=0,nademfix=nademfix,$
                                _extra=initnad.argsinitpar) $
               else parinit = $
                  call_function(initnad.fcninitpar,0,0,initnadem,$
                                initnad.nadabs_siglim,nadem_siglim,$
                                heifix=0,nademfix=nademfix)

;               tmpdat = (nadcube.dat)[i,j,*]
;;              This block looks for automatically-determined absorption and 
;;              emission line indices (from IFSF_CMPNADWEQ, invoked by IFSFA)
;;              and uses these to only fit the emission line region by setting
;;              anything blueward to 1. Note that if only an emission line was 
;;              found, then the index used is shifted blueward slightly to make
;;              sure the entire line is included.   
;
;               if (nadcube.iweq)[i,j,1] ne -1 then $
;                  tmpind = (nadcube.iweq)[i,j,1] $
;               else if (nadcube.iweq)[i,j,2] ne -1 then $
;                  tmpind = (nadcube.iweq)[i,j,2]-1 $
;               else $
;                  message,'No absorption or emission indices. Aborting.'
;               tmpdat[0:tmpind] = 1d
                                
               param = Mpfitfun(initnad.fcnfitnad,$
                                passwav,$
                                passdat,$
                                passerr,$
                                parinfo=parinit,perror=perror,maxiter=100,$
                                bestnorm=chisq_emonly,covar=covar,$
                                yfit=specfit,dof=dof_emonly,$
                                nfev=nfev,niter=niter_emonly,status=status,$
                                quiet=quiet,$
                                npegged=npegged,ftol=1D-6,errmsg=errmsg,$
                                functargs=argsfitnad)
               if status eq 5 then message,'Max. iterations reached.',/cont
               if status eq 0 OR status eq -16 then $
                  message,'Error in MPFIT. Aborting.'

;              Reject insignificant components automatically
               if ~ tag_exist(initnad,'noautoreject') then begin
                  
               igdindices = where((nadcube.iweq)[i,j,*] ne -1,ctgdindices)
               if ctgdindices gt 0 then begin
                  ilinelo = min((nadcube.iweq)[i,j,igdindices])
                  ilinehi = max((nadcube.iweq)[i,j,igdindices])
                  niline = n_elements((nadcube.wave)[i,j,*])
                  inotlinelo = dindgen(ilinelo-1)
                  inotlinehi = dindgen(niline-ilinehi-1) + ilinehi + 1
                  rms = sqrt(median(((nadcube.dat)[i,j,[inotlinelo,inotlinehi]]-1d)^2d))
               endif else begin
                  rms = sqrt(median(((nadcube.dat)[i,j,*]-1d)^2d))
               endelse
               sigthresh = 1d
               ctbad = 0
               ctgd = nnadem
               ibad = !NULL
               igd = !NULL
               for k=0,nnadem-1 do begin
                  karr = 3 + k*4 + dindgen(4)
                  modcomp = ifsf_nadfcn((nadcube.wave)[i,j,*],[0,0,1,param[karr]],$
                                        specres=specres)
                  if max(modcomp)-1d le rms*sigthresh then begin
                     ctbad++
                     ctgd--
                     ibad = [ibad,k]
                  endif else begin
                     igd = [igd,k]
                  endelse
               endfor
               if ctbad gt 0 then begin
                  nnadem -= ctbad
                  if ctgd gt 0 then begin
                     winit = reform(((initnad.nadem_zinit)[i,j,igd]+1d)$
                                    *linelist['NaD1'],nnadem)
                     siginit = reform((initnad.nadem_siginit)[i,j,igd],nnadem)
                     if tag_exist(initnad,'nadem_finit') then $
                        finit = reform((initnad.nadem_finit)[i,j,igd],nnadem) $
                     else finit = dblarr(nnadem)+0.1d
                     if tag_exist(initnad,'nadem_rinit') then $
                        rinit = reform((initnad.nadem_rinit)[i,j,igd],nnadem) $
                     else rinit = dblarr(nnadem)+nad_emrat_init
                     initnadem = [[winit],[siginit],[finit],[rinit]]
                     if tag_exist(initnad,'nadem_fix') then $
                        nademfix=reform((initnad.nadem_fix)[i,j,igd,*],nnadem,4) $
                     else nademfix=0b
                  endif else begin
                     refit = 0b
                     initnadem=0
                     nademfix=0b
                     dofirstemfit=0b
                  endelse
               endif else begin
                 refit = 0b
                 first_param = param
                 first_nademfix = nademfix
                 first_parinit = parinit
                 first_modflux = passdat
                 initnadem = transpose(reform(param[3:2+nnadem*4],4,nnadem))
                 nademfix=rebin(transpose([1b,1b,1b,1b]),nnadem,4)                
               endelse
               
               endif else begin
                  refit = 0b
                  first_param = param
                  first_nademfix = nademfix
                  first_parinit = parinit
                  first_modflux = passdat
                  initnadem = transpose(reform(param[3:2+nnadem*4],4,nnadem))
                  nademfix=rebin(transpose([1b,1b,1b,1b]),nnadem,4)
               endelse

               endwhile
                           
            endif
         endif else begin
            initnadem=0
            nademfix=0b
         endelse

;        Fill out parameter structure with initial guesses and constraints

         if (nnadem eq 0 AND nnadabs eq 0 AND nhei eq 0) then begin
            message,'No components specified. Skipping this spaxel.',/cont
            goto,nofit
         endif

         if nnadabs gt 0 AND tag_exist(initnad,'nadabs_fix') then $
            nadabsfix=reform((initnad.nadabs_fix)[i,j,0:nnadabs-1,*],nnadabs,4) $
         else nadabsfix=0b

         refit=1b
         while refit ne 0b do begin

         if tag_exist(initnad,'argsinitpar') then parinit = $
            call_function(initnad.fcninitpar,inithei,initnadabs,initnadem,$
                          initnad.nadabs_siglim,nadem_siglim,$
                          heifix=heifix,nadabsfix=nadabsfix,nademfix=nademfix,$
                          _extra=initnad.argsinitpar) $
         else parinit = $
            call_function(initnad.fcninitpar,inithei,initnadabs,initnadem,$
                          initnad.nadabs_siglim,nadem_siglim,$
                          heifix=heifix,nadabsfix=nadabsfix,nademfix=nademfix)

         passdat = reform((nadcube.dat)[i,j,ilamlo:ilamhi],ilamhi-ilamlo+1)
         param = Mpfitfun(initnad.fcnfitnad,$
                          passwav,$
                          passdat,$
                          passerr,$
                          parinfo=parinit,perror=perror,maxiter=100,$
                          bestnorm=chisq,covar=covar,yfit=specfit,dof=dof,$
                          nfev=nfev,niter=niter,status=status,quiet=quiet,$
                          npegged=npegged,ftol=1D-6,errmsg=errmsg,$
                          functargs=argsfitnad)
         if status eq 5 then message,'Max. iterations reached',/cont
         if status eq 0 OR status eq -16 then $
            message,'Error in MPFIT. Aborting.'

;         refit=0b
;        Reject insignificant components automatically
         if ~ tag_exist(initnad,'noautoreject') then begin

         if nnadabs gt 0 then begin
            igdindices = where((nadcube.iweq)[i,j,*] ne -1,ctgdindices)
            if ctgdindices gt 0 then begin
               ilinelo = min((nadcube.iweq)[i,j,igdindices])
               ilinehi = max((nadcube.iweq)[i,j,igdindices])
               niline = n_elements((nadcube.dat)[i,j,*])
               inotlinelo = dindgen(ilinelo-1)
               inotlinehi = dindgen(niline-ilinehi-1) + ilinehi + 1
               rms = sqrt(median(((nadcube.dat)[i,j,[inotlinelo,inotlinehi]]-1d)^2d))
            endif else begin
               rms = sqrt(median(((nadcube.dat)[i,j,*]-1d)^2d))
            endelse
            sigthresh = 1d
            ctbad = 0
            ctgd = nnadabs
            ibad = !NULL
            igd = !NULL
            for k=0,nnadabs-1 do begin
               karr = 3+nhei*3 + k*4 + dindgen(4)
               modcomp = ifsf_nadfcn((nadcube.wave)[i,j,*],[0,1,0,param[karr]],$
                                     specres=specres)
               if 1d - min(modcomp) le rms*sigthresh then begin
                  ctbad++
                  ctgd--
                  ibad = [ibad,k]
               endif else begin
                  igd = [igd,k]
               endelse
            endfor
            if ctbad gt 0 then begin
               nnadabs -= ctbad
               if nnadabs gt 0 then begin
                  if tag_exist(initnad,'nadabs_cfinit') then $
                     cfinit = reform((initnad.nadabs_cfinit)[i,j,igd],nnadabs) $
                  else cfinit = dblarr(nnadabs)+0.5d
                  if tag_exist(initnad,'nadabs_tauinit') then $
                  tauinit = reform((initnad.nadabs_tauinit)[i,j,igd],nnadabs) $
                  else tauinit = dblarr(nnadabs)+0.5d
                  winit = reform(((initnad.nadabs_zinit)[i,j,igd]$
                     +1d)*linelist['NaD1'],nnadabs)
                  siginit = reform((initnad.nadabs_siginit)$
                     [i,j,igd],nnadabs)
                  initnadabs = [[cfinit],[tauinit],[winit],[siginit]]
                  if nnadabs gt 0 AND tag_exist(initnad,'nadabs_fix') then $
                     nadabsfix=reform((initnad.nadabs_fix)[i,j,igd,*],nnadabs,4) $
                  else nadabsfix=0b
               endif else begin
                  if tag_exist(initnad,'nadem_fitinit') AND $
                     initnadem[0] ne 0 then begin
                     param = first_param
                     refit=0b
                  endif else begin
                     message,'Skipping this spectrum b/c all components rejected.',/cont
                     goto,nofit
                  endelse
               endelse
            endif else refit=0b
         endif else refit=0b
         
         endif else refit=0b
                  
         endwhile
  
;         endif else initnadabs=0

;        Plot fit
;        
;        If the data was not first processed with IFSF, then set the redshift
;        for plotting to be the systemic redshift. Otherwise, use the stellar 
;        redshift determined from the fit.
         if tag_exist(initnad,'zref') then zref = initnad.zref $
         else if ~ tag_exist(initnad,'noemlinfit') then zref = struct.zstar $
         else zref = initdat.zsys_gas
         if ~ tag_exist(initnad,'noemlinfit') then zstar=struct.zstar else $
         zstar = 0d
         if ~ noplot then begin
            if tag_exist(initnad,'argspltfitnad') then $
               ifsf_pltnadfit,(nadcube.wave)[i,j,*],$
                              (nadcube.dat)[i,j,*],$
                              (nadcube.err)[i,j,*],$
                              param,outfile+'_nad_fit',zref,$
                              specres=specres,zstar=zstar,$
                              _extra=initnad.argspltfitnad $
            else $
               ifsf_pltnadfit,(nadcube.wave)[i,j,*],$
                              (nadcube.dat)[i,j,*],$
                              (nadcube.err)[i,j,*],$
                              param,outfile+'_nad_fit',zref,$
                              specres=specres,zstar=zstar
         endif

;        Compute model equivalent widths
         weq=1
         nademflux=1
         modspec = ifsf_nadfcn(passwav,param,weq=weq,$
                               nademflux=nademflux,cont=(nadcube.cont)[i,j,*],$
                               specres=specres)

;        Compute errors in fit
         if ~ keyword_set(noerr) then begin
            if dofirstemfit then nademfix_use = first_nademfix $
            else nademfix_use = nademfix
            if keyword_set(nomc) then plotonly=1b else plotonly=0b
            errors = ifsf_fitnaderr([nhei,nnadabs,nnadem],passwav,$
                                    modspec,passerr,passcont,$
                                    parinit,$
                                    outfile+'_nad_errs.ps',outfile+'_nad_mc.xdr',$
                                    dofirstemfit=dofirstemfit,$
                                    first_parinit=first_parinit,$
                                    first_modflux=first_modflux,$
                                    heifix=heifix,nadabsfix=nadabsfix,$
                                    nademfix=nademfix_use,$
                                    niter=initnad.mcniter,$
                                    nsplit=nsplit,quiet=quiet,weqerr=weqerr,$
                                    nademfluxerr=nademfluxerr,noplot=noplot,$
                                    plotonly=plotonly,specres=specres)
         endif
         
         ifsf_printnadpar,nadparlun,i+1,j+1,param

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
;               NaD absorption line parameters
                cf: dblarr(ncols,nrows,maxncomp)+bad,$
                cferr: dblarr(ncols,nrows,maxncomp,2)+bad,$
                tau: dblarr(ncols,nrows,maxncomp)+bad,$
                tauerr: dblarr(ncols,nrows,maxncomp,2)+bad,$
                waveabs: dblarr(ncols,nrows,maxncomp)+bad,$
                waveabserr: dblarr(ncols,nrows,maxncomp,2)+bad,$
                sigmaabs: dblarr(ncols,nrows,maxncomp)+bad,$
                sigmaabserr: dblarr(ncols,nrows,maxncomp,2)+bad,$
;               NaD emission line parameters
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
         if ~ keyword_set(noerr) then begin
            nadfit.weqabserr[i,j,*]=weqerr[0,*]
            nadfit.weqemerr[i,j,*]=weqerr[1,*]
            nadfit.totfluxemerr[i,j,*]=nademfluxerr
         endif
         if nhei gt 0 then begin
            iarr = 3 + dindgen(nhei)*3
            nadfit.wavehei[i,j,0:nhei-1]=param[iarr]
            nadfit.sigmahei[i,j,0:nhei-1]=param[iarr+1]
            nadfit.fluxhei[i,j,0:nhei-1]=param[iarr+2]
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
         
copyvor: 

         if tag_exist(initdat,'vormap') then begin
            if vordone[vormap[i,j]-1] then begin
               tagnames = tag_names(nadfit)
               iref = vorrefnad[vormap[i,j]-1,0]
               jref = vorrefnad[vormap[i,j]-1,1]
               foreach tag,tags2d do begin
                  itag = where(tag eq tagnames)
                  nadfit.(itag)[i,j] = nadfit.(itag)[iref,jref]
               endforeach
               foreach tag,tags3d do begin
                  itag = where(tag eq tagnames)
                  nadfit.(itag)[i,j,*] = nadfit.(itag)[iref,jref,*]
               endforeach
               foreach tag,tags4d do begin
                  itag = where(tag eq tagnames)
                  nadfit.(itag)[i,j,*,*] = nadfit.(itag)[iref,jref,*,*]
               endforeach
               thislab = string(i+1,'_',j+1,format='(I04,A,I04)')
               reflab = string(iref+1,'_',jref+1,format='(I04,A,I04)')
               thisout = initdat.outdir+initdat.label+'_'+thislab
               refout = initdat.outdir+initdat.label+'_'+reflab
               if file_test(refout+'_nad_fit.jpg') then begin
                  file_copy,refout+'_nad_fit.jpg',thisout+'_nad_fit.jpg',/over
                  print,'Using reference coordinate: [col,row]=[',iref+1,',',jref+1,']',$
                        format='(A0,I0,A0,I0,A0)'
               endif else begin
                  print,'Reference coordinate ([col,row]=[',iref+1,',',jref+1,$
                        ']) has no plot file; not copying.',$
                        format='(A0,I0,A0,I0,A0)'
               endelse
            endif
         endif

         if tag_exist(initdat,'vormap') then begin
            if finite(initdat.vormap[i,j]) AND $
               initdat.vormap[i,j] ne bad then begin
               if not vordone[vormap[i,j]-1] then begin
                  vordone[vormap[i,j]-1] = 1b
                  vorrefnad[vormap[i,j]-1,*] = [i,j]
               endif
            endif
         endif

nofit:

      endfor
      
   endfor

   if ~ keyword_set(noxdr) then begin
      if tag_exist(initnad,'outxdr') then outxdr=initnad.outxdr $
      else outxdr=initdat.outdir+initdat.label+'.nadfit.xdr'
      save,nadfit,file=outxdr
   endif

   free_lun,nadparlun

   print,'Total runtime: ',systime(1)-starttime,' s.',$
         format='(/,A0,I0,A0,/)'

end
