; docformat = 'rst'
;
;+
;
; Fit Na D absorption line. Plot fit and write fit parameters to a
; file.
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
;    
; :Copyright:
;    Copyright (C) 2013 David S. N. Rupke
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
pro ifsf_fitnad,initproc,cols=cols,rows=rows,verbose=verbose

;  NaD optical depth ratio (blue to red)
   tratio = 2.0093d
   nad_emrat_init = 1.5d

   starttime = systime(1)
   time = 0
   if keyword_set(verbose) then quiet=0 else quiet=1

   ; Get fit initialization
   initnad={dumy: 1}
   initdat = call_function(initproc,initnad=initnad)

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
         if refline ne '' then begin
            if refline ne 'HeI5876' AND $
               ~ tag_exist(initnad,'noemlinfit') then begin
;              Get reference line parameters
               reflinelist = ifsf_linelist([refline])
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
         endif else inithei=0

;        Get NaD absorption parameters
         nnadabs = initnad.nnadabs[i,j]
         if nnadabs gt 0 then begin
            if tag_exist(initnad,'nadabs_cfinit') then $
               cfinit = (initnad.nadabs_cfinit)[i,j,0:nnadabs-1] $
            else cfinit = dblarr(nnadabs)+0.5d
            if tag_exist(initnad,'nadabs_tauinit') then $
               tauinit = (initnad.nadabs_tauinit)[i,j,0:nnadabs-1] $
            else tauinit = dblarr(nnadabs)+0.5d
            winit = reform(((initnad.nadabs_zinit)[i,j,0:nnadabs-1]$
                            +1d)*linelist['NaD1'],nnadabs)
            siginit = reform((initnad.nadabs_siginit)$
                              [i,j,0:nnadabs-1],nnadabs)
            initnadabs = [[cfinit],[tauinit],[winit],[siginit]]
         endif else initnadabs=0

;        Get NaD emission parameters
         nnadem = initnad.nnadem[i,j]
         if tag_exist(initnad,'nadem_fix') then nademfix=initnad.nadem_fix $
         else nademfix=0b
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
         endif else initnadem=0

;        Fill out parameter structure with initial guesses and constraints
         if tag_exist(initnad,'argsinitpar') then parinit = $
            call_function(initnad.fcninitpar,inithei,initnadabs,initnadem,$
                          initnad.nadabs_siglim,initnad.nadem_siglim,$
                          heifix=heifix,nademfix=nademfix,$
                          _extra=initnad.argsinitpar) $
         else parinit = $
            call_function(initnad.fcninitpar,inithei,initnadabs,initnadem,$
                          initnad.nadabs_siglim,initnad.nadem_siglim,$
                          heifix=heifix,nademfix=nademfix)

         param = Mpfitfun(initnad.fcnfitnad,$
                          (nadcube.wave)[i,j,*],$
                          (nadcube.dat)[i,j,*],$
                          (nadcube.err)[i,j,*],$
                          parinfo=parinit,perror=perror,maxiter=100,$
                          bestnorm=chisq,covar=covar,yfit=specfit,dof=dof,$
                          nfev=nfev,niter=niter,status=status,quiet=quiet,$
                          npegged=npegged,ftol=1D-6,functargs=argslinefit,$
                          errmsg=errmsg)
         if status eq 5 then print,'IFSF_FITNAD: Max. iterations reached.'
         if status eq 0 OR status eq -16 then begin
            print,'IFSF_FITNAD: Error in MPFIT. Aborting.'
            outstr = 0
            goto,finish
         endif

         if ~ tag_exist(initnad,'noemlinfit') then zuse = struct.zstar $
         else zuse = initdat.zsys_gas
         if tag_exist(initnad,'argspltfitnad') then $
            ifsf_pltnadfit,(nadcube.wave)[i,j,*],$
                           (nadcube.dat)[i,j,*],$
                           param,outfile+'_nad_fit',struct.zstar,$
                           _extra=initnad.argspltfitnad $
         else $
            ifsf_pltnadfit,(nadcube.wave)[i,j,*],$
                           (nadcube.dat)[i,j,*],$
                           param,outfile+'_nad_fit',struct.zstar
      
         ifsf_printnadpar,nadparlun,i+1,j+1,param                  

nofit:

      endfor
      
   endfor


finish:

   free_lun,nadparlun

end
