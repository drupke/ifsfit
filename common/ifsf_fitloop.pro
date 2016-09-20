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
;
; :Copyright:
;    Copyright (C) 2016 David S. N. Rupke
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

   printf,loglun,'[col,row]=[',colarr[ispax]+1,',',rowarr[ispax]+1,'] out of [',$
          cube.ncols,',',cube.nrows,']',format='(A0,I0,A0,I0,A0,I0,A0,I0,A0)'
   i = colarr[ispax]
   j = rowarr[ispax]

   if oned then begin
      flux = cube.dat[*,i]
      ;       absolute value takes care of a few deviant points
      err = sqrt(abs(cube.var[*,i]))
      bad = cube.dq[*,i]
   endif else begin
      flux = reform(cube.dat[i,j,*],cube.nz)
      err = reform(sqrt(abs(cube.var[i,j,*])),cube.nz)
      bad = reform(cube.dq[i,j,*],cube.nz)
   endelse

   ;    Apply DQ plane
   indx_bad = where(bad gt 0,ct)
   if ct gt 0 then begin
      flux[indx_bad] = 0d
      err[indx_bad] = max(err)*100d
   endif

   nodata = where(flux ne 0d,ct)
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

         ;       Initialize stellar redshift for this spaxel
         if oned then zstar = initdat.zinit_stars[i] $
         else zstar = initdat.zinit_stars[i,j]

         ;       Ignore NaI D line for purposes of continuum fit by maximizing
         ;       error. This doesn't actually remove these data points from the fit,
         ;       but it minimizes their weight.
         if not tag_exist(initdat,'keepnad') then begin
            if not tag_exist(initdat,'nad_contcutrange') then begin
               nadran_rest = [5850d,5900d]
               nadran = (1d + zstar) * nadran_rest
            endif else nadran = initdat.nad_contcutrange
            indx_nad = where(cube.wave ge nadran[0] AND $
               cube.wave le nadran[1],ct)
            if ct gt 0 then err[indx_nad]=max(err)
         endif

         ;       Option to tweak cont. fit
         if tag_exist(initdat,'tweakcntfit') then $
            tweakcntfit = reform(initdat.tweakcntfit[i,j,*,*],3,$
            n_elements(initdat.tweakcntfit[i,j,0,*])) $
         else tweakcntfit = 0

         ;       Initialize starting wavelengths
         linelistz = hash(initdat.lines)
         if ~ tag_exist(initdat,'noemlinfit') AND ct_comp_emlist gt 0 then $
            foreach line,initdat.lines do $
            ;              if oned then $
            ;                 linelistz[line] = $
            ;                 reform(linelist[line]*(1d + (initdat.zinit_gas)[line,i,*]),$
            ;                           initdat.maxncomp) $
            ;              else $
            linelistz[line] = $
            reform(linelist[line]*(1d + (initdat.zinit_gas)[line,i,j,*]),$
            initdat.maxncomp)

         structinit = ifsf_fitspec(cube.wave,flux,err,zstar,linelist,$
            linelistz,ncomp,initdat,quiet=quiet,$
            siglim_gas=siglim_gas,$
            tweakcntfit=tweakcntfit)

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
                  structinit.perror,structinit.parinfo)
               linelistz = linepars.wave
               if tag_exist(initdat,'masksig_secondfit') then $
                  masksig_secondfit = initdat.masksig_secondfit $
               else masksig_secondfit = masksig_secondfit_def
               maskwidths = hash(initdat.lines)
               foreach line,initdat.lines do $
                  maskwidths[line] = masksig_secondfit*linepars.sigma[line]
               maskwidths_tmp = maskwidths
               peakinit_tmp = linepars.fluxpk
               siginit_gas_tmp = linepars.sigma
            endif else begin
               maskwidths_tmp=0
               peakinit_tmp=0
               siginit_gas_tmp=0
            endelse

            struct = ifsf_fitspec(cube.wave,flux,err,structinit.zstar,$
               linelist,$
               linelistz,ncomp,initdat,quiet=quiet,$
               maskwidths=maskwidths_tmp,$
               peakinit=peakinit_tmp,$
               siginit_gas=siginit_gas_tmp,$
               siglim_gas=siglim_gas,$
               tweakcntfit=tweakcntfit)
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

         ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
         ; Save result to a file
         ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

         if ~ dofit AND ~ abortfit then begin
            if oned then $
               save,struct,file=$
               string(initdat.outdir,initdat.label,'_',i+1,'.xdr',$
               format='(A,A,A,I04,A,I04,A)') $
            else $
               save,struct,file=$
               string(initdat.outdir,initdat.label,'_',i+1,'_',j+1,$
               '.xdr',format='(A,A,A,I04,A,I04,A)')
         endif

      endwhile

   endif
   
   if keyword_set(logfile) then $
      if logfile[ispax] ne '' then $
         free_lun,loglun

end
