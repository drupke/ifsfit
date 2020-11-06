; docformat = 'rst'
;
;+
;
; Take outputs from IFSF and perform fitting loop. If loop is split among multiple
; cores, then DRT_BRIDGELOOP parses this file to feed into a batch file.
;
; :Categories:
;    IFSF
;
; :Returns:
;    None.
;    
; :Params:
;    ispax: in, required, type=int
;      Value of index over which to loop
;    colarr: in, required,type=intarr(2)
;      Column # of spaxel (0-offset)
;    rowarr: in, required, type=intarr(2)
;      Row # of spaxel (0-offset)
;    cube: in, required, type=structure
;      Output from IFSF_READCUBE, containing data
;    initdat: in, required, type=structure
;      Output from initialization routine, containing fit parameters
;    linelist: in, required, type=hash
;      Output from IFSF_LINELIST.
;    oned: in, required, type=byte
;      Whether data is in a cube or in one dimension (longslit)
;    onefit: in, required, type=byte
;      If set, ignore second fit
;    quiet: in, required, type=byte
;      verbosity switch from IFSF
;    
; :Keywords:
;    logfile: in, optional, type=strarr
;      Names of log files; one per spaxel.
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
;      2016sep18, DSNR, copied from IFSF into standalone procedure
;      2016sep26, DSNR, small change in masking for new treatment of spec. res.
;      2016oct20, DSNR, fixed treatment of SIGINIT_GAS
;      2016nov17, DSNR, added flux calibration
;      2018jun25, DSNR, added MC error calculation on stellar parameters
;
; :Copyright:
;    Copyright (C) 2016--2018 David S. N. Rupke
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
pro ifsf_fitloop,ispax,colarr,rowarr,cube,initdat,linelist,oned,onefit,quiet,$
                 logfile=logfile

;  Open logfile; loglun of -1 means STDOUT
   loglun=-1
   if keyword_set(logfile) then $
      if logfile[ispax] ne '' then $
         openw,loglun,logfile[ispax],/get_lun,/append

   masksig_secondfit_def = 2d

   i = colarr[ispax]
   j = rowarr[ispax]
   printf,loglun,'[col,row]=[',i+1,',',j+1,'] out of [',$
          cube.ncols,',',cube.nrows,']',format='(A0,I0,A0,I0,A0,I0,A0,I0,A0)'

   if oned then begin
      flux = cube.dat[*,i]
;     absolute value takes care of a few deviant points
      err = sqrt(abs(cube.var[*,i]))
      dq = cube.dq[*,i]
   endif else begin
      flux = reform(cube.dat[i,j,*],cube.nz)
      err = reform(sqrt(abs(cube.var[i,j,*])),cube.nz)
      dq = reform(cube.dq[i,j,*],cube.nz)
   endelse
;  Error maximum, for use later
   errmax = max(err)

   if tag_exist(initdat,'vormap') then begin
     tmpi = cube.vorcoords[i,0]
     tmpj = cube.vorcoords[i,1]
     i = tmpi
     j = tmpj     
     printf,loglun,'Reference coordinate: [col,row]=[',i+1,',',j+1,']',$
            format='(A0,I0,A0,I0,A0)'
     
   endif

   if oned then begin
      outlab = string(initdat.outdir,initdat.label,'_',i+1,$
         format='(A,A,A,I04,A,I04)')
   endif else begin
      outlab = string(initdat.outdir,initdat.label,'_',i+1,'_',j+1,$
         format='(A,A,A,I04,A,I04)')
   endelse


;;  Apply flux calibration
;   if tag_exist(initdat,'fluxunits') then begin
;      flux*=initdat.fluxunits
;      err*=initdat.fluxunits
;   endif
;
;;  Normalize
;   fluxmed = median(flux)
;   flux/=fluxmed
;   err/=fluxmed

;  Apply DQ plane
   indx_bad = where(dq gt 0,ct)
   if ct gt 0 then begin
      flux[indx_bad] = 0d
      err[indx_bad] = errmax*100d
   endif

   nodata = where(flux ne 0d AND finite(flux),ct)
   if ct ne 0 then begin

      if ~ tag_exist(initdat,'noemlinfit') then begin

         ;          Extract # of components and initial redshift guesses
         ;          specific to this spaxel, and write as hashes.
         ncomp = hash(initdat.lines)
         foreach line,initdat.lines do $
            if oned then ncomp[line] = (initdat.ncomp)[line,i] $
            else ncomp[line] = (initdat.ncomp)[line,i,j]
      endif

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; First fit
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      dofit=1b
      abortfit=0b
      while(dofit) do begin

         if ~ tag_exist(initdat,'noemlinfit') then $
            nocomp_emlist = $
               ncomp.where(0,complement=comp_emlist,ncomp=ct_comp_emlist) $
         else ct_comp_emlist=0

         if tag_exist(initdat,'siglim_gas') then begin
            size_siglim = size(initdat.siglim_gas)
            if size_siglim[0] eq 1 then siglim_gas = initdat.siglim_gas $
            else begin
               if oned then siglim_gas = initdat.siglim_gas[i,*] $
               else siglim_gas = initdat.siglim_gas[i,j,*]
            endelse
         endif else siglim_gas = 0b
         if tag_exist(initdat,'siginit_gas') then begin
            size_siginit = size(initdat.siginit_gas[initdat.lines[0]])
            if size_siginit[0] eq 1 then siginit_gas = initdat.siginit_gas $
            else begin
               siginit_gas = hash()
               if oned then $
                  foreach key,initdat.lines do $
                     siginit_gas[key] = initdat.siginit_gas[key,i,*] $
               else $
                  foreach key,initdat.lines do $
                     siginit_gas[key] = initdat.siginit_gas[key,i,j,*]
            endelse
         endif else siginit_gas = 0b

;        Initialize stellar redshift for this spaxel
         if oned then zstar = initdat.zinit_stars[i] $
         else zstar = initdat.zinit_stars[i,j]
         zstar_init = zstar

;        Regions to ignore in fitting. Set to max(err).
         if tag_exist(initdat,'cutrange') then begin
            sizecutrange = size(initdat.cutrange)
            if sizecutrange[0] eq 1 then begin
               indx_cut = where(cube.wave ge initdat.cutrange[0] AND $
                                cube.wave le initdat.cutrange[1],ct)
               if ct gt 0 then begin
                  dq[indx_cut]=1b
                  err[indx_cut]=errmax*100d
               endif
            endif else if sizecutrange[0] eq 2 then begin
               for k=0,sizecutrange[2]-1 do begin
                  indx_cut = where(cube.wave ge initdat.cutrange[0,k] AND $
                     cube.wave le initdat.cutrange[1,k],ct)
                  if ct gt 0 then begin
                     dq[indx_cut]=1b
                     err[indx_cut]=errmax*100d
                  endif
               endfor
            endif else message,'CUTRANGE not properly specified.'
         endif

;        Option to tweak cont. fit
         if tag_exist(initdat,'tweakcntfit') then $
            tweakcntfit = reform(initdat.tweakcntfit[i,j,*,*],3,$
            n_elements(initdat.tweakcntfit[i,j,0,*])) $
         else tweakcntfit = 0

;        Initialize starting wavelengths
         linelistz = hash()
         if ~ tag_exist(initdat,'noemlinfit') AND ct_comp_emlist gt 0 then $
            foreach line,initdat.lines do $
               if oned then $
                  linelistz[line] = $
                     reform(linelist[line]*(1d + (initdat.zinit_gas)[line,i,*]),$
                            initdat.maxncomp) $
               else $
                  linelistz[line] = $
                     reform(linelist[line]*(1d + (initdat.zinit_gas)[line,i,j,*]),$
                            initdat.maxncomp)
         linelistz_init = linelistz

         if ~ quiet then print,'IFSF_FITLOOP: First call to IFSF_FITSPEC'
         structinit = ifsf_fitspec(cube.wave,flux,err,dq,zstar,linelist,$
                                   linelistz,ncomp,initdat,quiet=quiet,$
                                   siglim_gas=siglim_gas,$
                                   siginit_gas=siginit_gas,$
                                   tweakcntfit=tweakcntfit,col=i+1,$
                                   row=j+1)

         testsize = size(structinit)
         if testsize[0] eq 0 then begin
            printf,loglun,'IFSF: Aborting.'
            abortfit=1b
         endif
         if not quiet then print,'FIT STATUS: ',structinit.fitstatus
         if structinit.fitstatus eq -16 then begin
            printf,loglun,'IFSF: Aborting.'
            abortfit=1b
         endif

         ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
         ; Second fit
         ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

         if ~ onefit AND ~ abortfit then begin

            if ~ tag_exist(initdat,'noemlinfit') AND ct_comp_emlist gt 0 then begin
               ;             Set emission line mask parameters
               linepars = ifsf_sepfitpars(linelist,structinit.param,$
                                          structinit.perror,$
                                          structinit.parinfo)
               linelistz = linepars.wave
               if tag_exist(initdat,'masksig_secondfit') then $
                  masksig_secondfit = initdat.masksig_secondfit $
               else masksig_secondfit = masksig_secondfit_def
               maskwidths = hash(initdat.lines)
               foreach line,initdat.lines do $
                  maskwidths[line] = masksig_secondfit*linepars.sigma_obs[line]
               maskwidths_tmp = maskwidths
               peakinit_tmp = linepars.fluxpk_obs
               siginit_gas_tmp = linepars.sigma
            endif else begin
               maskwidths_tmp=0
               peakinit_tmp=0
               siginit_gas_tmp=0
            endelse

            zstar_init2 = structinit.zstar
            if ~ quiet then print,'IFSF_FITLOOP: Second call to IFSF_FITSPEC'
            struct = ifsf_fitspec(cube.wave,flux,err,dq,structinit.zstar,$
                                  linelist,$
                                  linelistz,ncomp,initdat,quiet=quiet,$
                                  maskwidths=maskwidths_tmp,$
                                  peakinit=peakinit_tmp,$
                                  siginit_gas=siginit_gas_tmp,$
                                  siglim_gas=siglim_gas,$
                                  tweakcntfit=tweakcntfit,col=i+1,$
                                  row=j+1)
            testsize = size(struct)
            if testsize[0] eq 0 then begin
               printf,loglun,'IFSF: Aborting.'
               abortfit=1b
            endif
            if not quiet then print,'FIT STATUS: ',struct.fitstatus
            if struct.fitstatus eq -16 then begin
               printf,loglun,'IFSF: Aborting.'
               abortfit=1b
            endif

         endif else struct = structinit

         ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
         ; Check components
         ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

         if tag_exist(initdat,'fcncheckcomp') AND $
            ~ tag_exist(initdat,'noemlinfit') AND $
            ~ onefit AND $
            ~ abortfit AND $
            ct_comp_emlist gt 0 then begin

            siglim_gas = struct.siglim

            linepars = ifsf_sepfitpars(linelist,struct.param,$
                                       struct.perror,struct.parinfo)
            if tag_exist(initdat,'argscheckcomp') then goodcomp = $
               call_function(initdat.fcncheckcomp,linepars,initdat.linetie,$
                             ncomp,newncomp,siglim_gas,$
                             _extra=initdat.argscheckcomp) $
            else goodcomp = $
               call_function(initdat.fcncheckcomp,linepars,initdat.linetie,$
                             ncomp,newncomp,siglim_gas)

            if newncomp.count() gt 0 then begin
               foreach nc,newncomp,line do $
                  printf,loglun,'IFSF: Repeating the fit of ',line,$
                         ' with ',string(nc,format='(I0)'),' components.',$
                         format='(5A0)'
            endif else begin
               dofit=0b
            endelse

         endif else dofit=0b


         if ~ dofit AND ~ abortfit then begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Monte Carlo the errors
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

            if tag_exist(initdat,'mcerrors') then begin

;              Parameter values to collect
               mcdist = hash()
               mcdist['ct_ebv'] = dblarr(initdat.mcerrors.niter)
               mcdist['ct_ppxf_sigma'] = dblarr(initdat.mcerrors.niter)
               mcdist['zstar'] = dblarr(initdat.mcerrors.niter)
               nwave = n_elements(struct.wave)
               rans = randomn(seed,nwave,initdat.mcerrors.niter,/double)
               model = struct.cont_fit + struct.emlin_fit

               for k=0,initdat.mcerrors.niter-1 do begin
;                 attempt to correct for over- and under-estimated errors
;                 by multiplying by reduced chi-squared
                  if struct.ct_rchisq ne 0d then $
                     erradj = struct.spec_err * struct.ct_rchisq $
                  else erradj = struct.spec_err
                  moderr = model + rans[*,k]*erradj
                  imaxerr = where(struct.spec_err eq max(err))
                  moderr[imaxerr] = 0d

;; Option 1: re-run entire continuum + emission-line fit, two iterations
;                  zstar = zstar_init
;                  linelistz = linelistz_init
;                  mcstructinit = ifsf_fitspec(struct.wave,moderr,struct.spec_err,$
;                                              dq[struct.fitran_indx],zstar,$
;                                              linelist,linelistz,$
;                                              ncomp,initdat,quiet=quiet,$
;                                              siglim_gas=siglim_gas,$
;                                              siginit_gas=siginit_gas,$
;                                              tweakcntfit=tweakcntfit,col=i+1,$
;                                              row=j+1)
;                  
;                  if ~ onefit then begin
;
;                     if ~ tag_exist(initdat,'noemlinfit') AND $
;                        ct_comp_emlist gt 0 then begin
;                        linepars = ifsf_sepfitpars(linelist,mcstructinit.param,$
;                                                   mcstructinit.perror,$
;                                                   mcstructinit.parinfo)
;                        linelistz = linepars.wave
;                        if tag_exist(initdat,'masksig_secondfit') then $
;                           masksig_secondfit = initdat.masksig_secondfit $
;                        else masksig_secondfit = masksig_secondfit_def
;                        maskwidths = hash(initdat.lines)
;                        foreach line,initdat.lines do $
;                           maskwidths[line] = $
;                              masksig_secondfit*linepars.sigma_obs[line]
;                        maskwidths_tmp = maskwidths
;                        peakinit_tmp = linepars.fluxpk_obs
;                        siginit_gas_tmp = linepars.sigma
;                     endif else begin
;                        maskwidths_tmp=0
;                        peakinit_tmp=0
;                        siginit_gas_tmp=0
;                     endelse
;
;                     mcstruct = ifsf_fitspec(struct.wave,moderr,struct.spec_err,$
;                                             dq[struct.gd_indx],zstar,$
;                                             linelist,linelistz,$
;                                             ncomp,initdat,quiet=quiet,$
;                                             siglim_gas=siglim_gas,$
;                                             maskwidths=maskwidths_tmp,$
;                                             peakinit=peakinit_tmp,$
;                                             siginit_gas=siginit_gas_tmp,$
;                                             tweakcntfit=tweakcntfit,col=i+1,$
;                                             row=j+1)
;
;
;                  endif else mcstruct = mcstructinit
            
; Option 2: re-run second stellar fit only
                  if ~ onefit then begin
                     zstar = zstar_init2
                  endif else begin
                     zstar = zstar_init
                     maskwidths_tmp = 0d
                     peakinit_tmp = 0d
                     siginit_gas_tmp = 0d
                  endelse
                  if ~ tag_exist(initdat,'doemlinmask') then $
                     initdat_tmp = create_struct(initdat,'doemlinmask',1b) $
                  else initdat_tmp = initdat
                  if ~ tag_exist(initdat,'noemlinfit') then $
                     initdat_tmp = create_struct(initdat_tmp,'noemlinfit',1b)

; For correctness, error here should be erradj, but for consistency with 
; first call we'll leave it alone. Difference is a constant scaling factor,
; so should in principle only affect reduced chi squared.
                  mcstruct = ifsf_fitspec(struct.wave,moderr,$
                                          struct.spec_err,$
                                          dq[struct.fitran_indx],zstar,$
                                          linelist,linelistz,$
                                          ncomp,initdat_tmp,quiet=quiet,$
                                          siglim_gas=siglim_gas,$
                                          maskwidths=maskwidths_tmp,$
                                          peakinit=peakinit_tmp,$
                                          siginit_gas=siginit_gas_tmp,$
                                          tweakcntfit=tweakcntfit,col=i+1,$
                                          row=j+1)



                  mcdist['zstar',k] = mcstruct.zstar
                  mcdist['ct_ppxf_sigma',k] = mcstruct.ct_ppxf_sigma
                  mcdist['ct_ebv',k] = mcstruct.ct_ebv
                  
               endfor

               ct_errors = ifsf_pltmcstel(mcdist,outlab+'_cterrs.pdf')
               struct = create_struct(struct,'ct_errors',ct_errors)
               
            endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Save result to a file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

            save,struct,file=outlab+'.xdr'

         endif

      endwhile

   endif else $
      print,'IFSF_FITLOOP: No good data. Aborting.'
   
   if keyword_set(logfile) then $
      if logfile[ispax] ne '' then $
         free_lun,loglun

end
