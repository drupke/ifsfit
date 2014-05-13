; docformat = 'rst'
;
;+
;
; Fit Na D absorption line.
;
; :Categories:
;    IFSFIT
;
; :Returns:    
;
; :Params:
;    
;    
; :Keywords:
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

   starttime = systime(1)
   time = 0
   if keyword_set(verbose) then quiet=0 else quiet=1

   ; Get fit initialization
   initdat = call_function(initproc)

   ; Get linelist
   linelist = ifsf_linelist(['NaD1','NaD2','HeI5876'])

   if tag_exist(initdat,'nad_taumax') then taumax=initdat.nad_taumax $
   else taumax = 5d
   if tag_exist(initdat,'nad_fitran') then fitran=initdat.nad_fitran

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   restore,file=initdat.outdir+initdat.label+'.lin.xdr'
   restore,file=initdat.outdir+initdat.label+'.nadspec.xdr'
   nadsize = size(nadcube.wave)

;  if keyword_set(sigmax) then bmax=sigmax*sqrt(2d) $
;  else bmax = 1000d/(2d*sqrt(alog(2d)))

   ncols = nadsize[1]
   nrows = nadsize[2]
   nz = nadsize[3]
   refcol = initdat.nad_refcoords[0]
   refrow = initdat.nad_refcoords[1]

;;  Arrays to hold flags
;   donefit = dblarr(ncols,nrows)
;   nabscomp = dblarr(ncols,nrows)+1
;   nemcomp = intarr(ncols,nrows)
;
;;  Array of distances from reference spectrum
;   cols = rebin(dindgen(ncols)+1,ncols,nrows)
;   rows = rebin(reform(dindgen(nrows)+1,1,nrows),ncols,nrows)
;   dref = sqrt((cols-refcol)^2d + (rows-refrow)^2d)
;;  ... properly sorted in order of increasing distance
;   dref_isrt = sort(dref)
;   nspec = n_elements(rows)
;   cols_srt = reform(cols(dref_isrt),nspec)
;   rows_srt = reform(rows(dref_isrt),nspec)
;   dref_srt = reform(dref(dref_isrt),nspec)

;  Fit reference spectrum

   print,'Fitting [',refcol,',',refrow,'].',format='(A,I03,A,I03,A)'

;  Get HeI parameters
   refline = initdat.nad_heitie[refcol-1,refrow-1]
   if refline ne '' then begin

;     Restore fit    
      lab = string(refcol+1,'_',refrow+1,format='(I04,A,I04)')
      infile = initdat.outdir+initdat.label+'_'+lab+'.xdr'
      outfile = initdat.outdir+initdat.label+'_'+lab
      if ~ file_test(infile) then begin
         print,'IFSFA: No XDR file for ',i+1,', ',j+1,'.',$
            format='(A0,I4,A0,I4,A0)'
         goto,nofit
      endif
      restore,file=infile
      reflinelist = ifsf_linelist([refline])
      linepars = ifsf_sepfitpars(reflinelist,struct.param,struct.perror,$
                                 struct.parinfo)

;     Get reference line parameters      
      iem = where(linepars.flux[refline,*] ne 0d,nhei)
      if nhei gt 0 then $
         inithei = [[(linepars.wave)[refline,0:nhei-1]/$
                     reflinelist[refline]*linelist['HeI5876']],$
                   [(linepars.sigma)[refline,0:nhei-1]],$
                   [dblarr(nhei)+0.1d]] $
      else inithei=0

   endif else inithei=0

;  Get NaD absorption parameters
   nnadabs = initdat.nad_nnadabs[refcol-1,refrow-1]
   if nnadabs gt 0 then begin
      if tag_exist(initdat,'nad_nadabs_cfinit') then $
         cfinit = (initdat.nad_nadabs_cfinit)[refcol-1,refrow-1,0:nnadabs-1] $
      else cfinit = dblarr(nnadabs)+0.5d
      if tag_exist(initdat,'nad_nadabs_tauinit') then $
         tauinit = (initdat.nad_nadabs_tauinit)[refcol-1,refrow-1,0:nnadabs-1] $
      else tauinit = dblarr(nnadabs)+0.5d
      winit = reform(((initdat.nad_nadabs_zinit)[refcol-1,refrow-1,0:nnadabs-1]$
                     +1d)*linelist['NaD2'],nnadabs)
      siginit = reform((initdat.nad_nadabs_siginit)$
                       [refcol-1,refrow-1,0:nnadabs-1],nnadabs)
      initnadabs = [[cfinit],[tauinit],[winit],[siginit]]
   endif else initnadabs=0

;  Get NaD emission parameters
   nnadem = initdat.nad_nnadem[refcol-1,refrow-1]
   if nnadem gt 0 then begin
      winit = reform(((initdat.nad_nadem_zinit)[refcol-1,refrow-1,0:nnadem-1]$
                     +1d)*linelist('NaD2'),nnadem)
      siginit = reform((initdat.nad_nadem_siginit)$
                       [refcol-1,refrow-1,0:nnadem-1],nnadem)
      if tag_exist(initdat,'nad_nadem_finit') then $
         finit = (initdat.nad_nadem_finit)[refcol-1,refrow-1,0:nnadem-1] $
      else finit = dblarr(nnadem)+0.1d
      if tag_exist(initdat,'nad_nadem_rinit') then $
         rinit = (initdat.nad_nadem_rinit)[refcol-1,refrow-1,0:nnadem-1] $
      else rinit = dblarr(nnadem)+1.5d
      initnadem = [[winit],[siginit],[finit],[rinit]]
   endif else initnadem=0

;  Fill out parameter structure with initial guesses and constraints
   if tag_exist(initdat,'nad_argsinitpar') then parinit = $
      call_function(initdat.nad_fcninitpar,inithei,initnadabs,initnadem,$
                    initdat.nad_nadabs_siglim,initdat.nad_nadem_siglim,$
                    _extra=initdat.nad_argsinitpar) $
   else parinit = $
      call_function(initdat.nad_fcninitpar,inithei,initnadabs,initnadem,$
                    initdat.nad_nadabs_siglim,initdat.nad_nadem_siglim)

   param = Mpfitfun(initdat.nad_fcnfitnad,$
                    (nadcube.wave)[refcol-1,refrow-1,*],$
                    (nadcube.dat)[refcol-1,refrow-1,*],$
                    (nadcube.err)[refcol-1,refrow-1,*],$
                    parinfo=parinit,perror=perror,maxiter=100,$
                    bestnorm=chisq,covar=covar,yfit=specfit,dof=dof,$
                    nfev=nfev,niter=niter,status=status,quiet=quiet,$
                    npegged=npegged,ftol=1D-6,functargs=argslinefit,$
                    errmsg=errmsg)
  if status eq 0 OR status eq -16 then begin
     print,'IFSF_FITSPEC: Error in MPFIT. Aborting.'
     outstr = 0
     goto,finish
  endif

nofit:

finish:

end
