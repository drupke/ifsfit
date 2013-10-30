pro gmos_fit_nad,gal,bin,sigfix=sigfix,taumax=taumax,sigmax=sigmax

;
; History
;  10jul21  DSNR  created
;

  if ~keyword_set(taumax) then taumax=5d

  if keyword_set(sigmax) then bmax=sigmax*sqrt(2d) $
  else bmax = 1000d/(2d*sqrt(alog(2d)))

  he_rest = 5875.661d
  nad1_rest = 5895.92d
  srcdir = '/Users/drupke/winds/routines_fitting/gmos/'

  initdat = gmos_initfit_spectra(gal,bin=bin)
  ncompinit = initdat.ncomp
  infile = initdat.infile
  fitdir = initdat.outdir
  fitran = initdat.nadfitran
  dx = initdat.dx
  dy = initdat.dy
  refx = initdat.nadxref
  refy = initdat.nadyref
; Equivalent width sigma thresholds for fitting 1 or 2 components to NaD
  nad_weq_sig_thresh_1comp = initdat.nad1thresh
  nad_weq_sig_thresh_2comp = initdat.nad2thresh
  zinit = initdat.zinit

; Arrays to hold flags
  donefit = dblarr(dx,dy)
  nabscomp = dblarr(dx,dy)+1
  nemcomp = intarr(dx,dy)
  if gal eq 'f08572nw' OR $
     gal eq 'f08572se' OR $
     gal eq 'f10565' OR $
     gal eq 'mrk231' OR $
     gal eq 'f17207w' OR $
     gal eq 'f17207e' then begin
     print,'Using 0 emission components.'
  endif else if gal eq 'mrk273' OR $
     gal eq 'vv705nw' OR $
     gal eq 'vv705se' then begin
        nemcomp = intarr(dx,dy)+1
        print,'Using 1 emission components.'
     endif else begin
     print,'GMOS_FIT_NAD: Galaxy not recognized.  Select # of emission components.'
     stop
  endelse

; Array of distances from reference spectrum
  cols = rebin(dindgen(dx)+1,dx,dy)
  rows = rebin(reform(dindgen(dy)+1,1,dy),dx,dy)
  dref = sqrt((cols-refx)^2d + (rows-refy)^2d)
; ... properly sorted in order of increasing distance
  dref_isrt = sort(dref)
  nspec = n_elements(rows)
  cols_srt = reform(cols(dref_isrt),nspec)
  rows_srt = reform(rows(dref_isrt),nspec)
  dref_srt = reform(dref(dref_isrt),nspec)

; Pick out spectra to fit using Weq criterion
  dofit = dblarr(dx,dy)
  readcol300,fitdir+gal+'.lines.dat',$
             col_f,row_f,comp_f,o1,o1e,ha,hae,n2,n2e,$
             nad_weq,nad_weqe,$
             /skip,/silent,format='(I,I,I,D,D,D,D,D,D,D,D)'
  ndofit=0
  for i=0,n_elements(col_f)-1 do begin
     if comp_f[i] eq 1 then begin
        if nad_weq[i] ge nad_weq_sig_thresh_1comp*nad_weqe[i] then begin
           dofit[col_f[i]-1,row_f[i]-1] = 1
           ndofit++
           if nad_weq[i] ge nad_weq_sig_thresh_2comp*nad_weqe[i] then $
              nabscomp[col_f[i]-1,row_f[i]-1] = 2
           haflux = ha[i]
           if comp_f[i+1] eq 2 then haflux += ha[i+1]
;  Parameters for Mrk 231 fit
           if gal eq 'mrk231' then begin
              if (haflux gt 0 AND $
                  nad_weq[i] ge 6d*nad_weqe[i] AND $
                  (col_f[i] le 13 OR row_f[i] ge 5)) $
              then nemcomp[col_f[i]-1,row_f[i]-1] = 1
           endif
        endif
     endif
  endfor
  dofit = reform(dofit(dref_isrt),nspec)
  idofit=1

; Fit reference spectrum

  print,'Fitting [',refx,',',refy,'].  Spec. ',idofit,' of ',ndofit,'.',$
        format='(A,I03,A,I03,A,I0,A,I0,A)'
  
  spec = string(refx,'_',refy,format='(I04,A0,I04)')
  spectot= fitdir+gal+'_'+spec
  parin  = spectot+'_nad_parin.dat'
  parout = spectot+'_nad_parout.dat'
  specin = spectot+'_nad_spec.dat'
  specout= spectot+'_nad_fit'
  xdr = spectot+'.xdr'

  restore,file=xdr

; Get emission-line parameters
  nem = nemcomp[refx-1,refy-1]
  if nem gt 0 then begin
     if n_elements(struct.param) gt 1 then $
        nem_emfit = struct.param[1] $
     else $
        nem_emfit = 0
     if nem_emfit gt 0 then begin
        linepars = sepfitpars(struct.param,struct.perror)
        haind = where(struct.linelabel eq 'Halpha')
        i_hapk = 0
        if nem_emfit gt 1 then begin
           for i=1,nem_emfit-1 do begin
              if linepars.fluxpk[haind,i] gt $
                 linepars.fluxpk[haind,0] then $
                    i_hapk = i
           endfor
        endif
        ; if nem eq 2 then begin
        ;    if nem_emfit gt 2 then begin
        ;    for i=1,nem_emfit-1 do begin
        ;       if linepars.fluxpk[haind,i] gt $
        ;          linepars.fluxpk[haind,0] then $
        ;             i_hapk = i
        ;    endfor
        ; endif
        he_lam = he_rest * (1d + struct.z.gas[i_hapk])
        he_b = linepars.sigma[haind,i_hapk] * sqrt(2d) / $
               he_lam * 299792d
     endif
  endif

; First with 2 components ...

  nabs = 2

  nadr = (nad1_rest*(1d + zinit[refx-1,refy-1,0])) + [-10,0,5]
  nadb = (nad1_rest*(1d + zinit[refx-1,refy-1,0])) + [-18,-10,-2]

  openw,lun,parin,/get_lun
  printf,lun,nabs,nem,fitran[0],fitran[1],format='(I8,I8,D8.1,D8.1)'
  fparin = '(D8.2,D8.2,D8.2,I8)'
  printf,lun,0.01,0.5,1.0,1,format=fparin
  printf,lun,0.01,0.5,taumax,1,format=fparin
;  printf,lun,6125,6135,6140,1,format=fparin
  printf,lun,nadr[0],nadr[1],nadr[2],1,format=fparin
  printf,lun,60,400,bmax,1,format=fparin
  printf,lun,0.01,0.5,1.0,1,format=fparin
  printf,lun,0.01,0.5,taumax,1,format=fparin
;  printf,lun,6135,6143,6155,1,format=fparin
  printf,lun,nadb[0],nadb[1],nadb[2],1,format=fparin
  printf,lun,60,200,bmax,1,format=fparin
  if nem eq 1 then begin
     printf,lun,0,0.1,2,1,format=fparin
     printf,lun,0,he_lam,0,0,format=fparin
     printf,lun,0,he_b,0,0,format=fparin
  endif
  free_lun,lun
  
  file_copy,srcdir+'lmfit.h',srcdir+'lmfit.h.BAK',/over
  lmfith = strarr(15)
  dummy = ''
  openr,lun,srcdir+'lmfit.h',/get_lun
  for j=0,11 do begin
     readf,lun,dummy
     lmfith[j] = dummy        
  endfor
  lmfith[12] = '#define FILE_PARI "'+parin+'"'
  lmfith[13] = '#define FILE_OUT  "'+parout+'"'
  lmfith[14] = '#define FILE_DAT  "'+specin+'"'
  free_lun,lun
  
  openw,lun,srcdir+'lmfit.h',/get_lun
  for j=0,14 do printf,lun,lmfith[j]
  free_lun,lun
     
  spawn,'gcc -o '+fitdir+'lmfit '+srcdir+'lmfit.c -I'+$
        srcdir+'/lib -lm'
  spawn,fitdir+'lmfit'
  
  gmos_plotnadfit,specin,parout,specout,struct.z

  spawn,'mv '+spectot+'_nad_parout.dat'+' '+spectot+'_nad_parout.2comp.dat'
  spawn,'mv '+specout+'.jpg '+specout+'.2comp.jpg'

; ... then with 1.
  
  nabs = 1

  nadr = (nad1_rest*(1d + zinit[refx-1,refy-1,0])) + [-15,0,10]

  openw,lun,parin,/get_lun
  printf,lun,nabs,nem,fitran[0],fitran[1],format='(I8,I8,D8.1,D8.1)'
  fparin = '(D8.2,D8.2,D8.2,I8)'
  printf,lun,0.01,0.5,1.0,1,format=fparin
  printf,lun,0.01,0.5,taumax,1,format=fparin
;  printf,lun,6125,6140,6150,1,format=fparin
  printf,lun,nadr[0],nadr[1],nadr[2],1,format=fparin
  if keyword_set(sigfix) then $
     printf,lun,60,sigfix*sqrt(2d),bmax,0,format=fparin $
  else $
     printf,lun,60,400,bmax,1,format=fparin
  if nem eq 1 then begin
     printf,lun,0,0.1,2,1,format=fparin
     printf,lun,0,he_lam,0,0,format=fparin
     printf,lun,0,he_b,0,0,format=fparin
  endif
  free_lun,lun
  
  file_copy,srcdir+'lmfit.h',srcdir+'lmfit.h.BAK',/over
  lmfith = strarr(15)
  dummy = ''
  openr,lun,srcdir+'lmfit.h',/get_lun
  for j=0,11 do begin
     readf,lun,dummy
     lmfith[j] = dummy        
  endfor
  lmfith[12] = '#define FILE_PARI "'+parin+'"'
  lmfith[13] = '#define FILE_OUT  "'+parout+'"'
  lmfith[14] = '#define FILE_DAT  "'+specin+'"'
  free_lun,lun
  
  openw,lun,srcdir+'lmfit.h',/get_lun
  for j=0,14 do printf,lun,lmfith[j]
  free_lun,lun
     
  spawn,'gcc -o '+fitdir+'lmfit '+srcdir+'lmfit.c -I'+$
        srcdir+'/lib -lm'
  spawn,fitdir+'lmfit'

  gmos_plotnadfit,specin,parout,specout,struct.z

  donefit[refx-1,refy-1] = 1
  cols_done = [refx]
  rows_done = [refy]
  nabs_done = [2]
  idofit++

; Fit rest of spectra
  
  for i=1,nspec-1 do begin
     
     if dofit[i] then begin

        badcomp=0
        nabs = nabscomp[cols_srt[i]-1,rows_srt[i]-1]
        nem = nemcomp[cols_srt[i]-1,rows_srt[i]-1]
       
        print,'Fitting [',cols_srt[i],',',rows_srt[i],'].',$
              format='(A,I03,A,I03,A,$)'

        spec = string(cols_srt[i],'_',rows_srt[i],$
                      format='(I04,A0,I04)')
        spectot= fitdir+gal+'_'+spec
        parin  = spectot+'_nad_parin.dat'
        parout = spectot+'_nad_parout.dat'
        specin = spectot+'_nad_spec.dat'
        specout= spectot+'_nad_fit'
        xdr = spectot+'.xdr'

;       Get emission-line parameters
        restore,file=xdr
        if nem gt 0 then begin
           if n_elements(struct.param) gt 1 then $
              nem_emfit = struct.param[1] $
           else $
              nem_emfit = 0
           if nem_emfit gt 0 then begin
              linepars = sepfitpars(struct.param,struct.perror)
              haind = where(struct.linelabel eq 'Halpha')
              i_hapk = 0
              if nem_emfit eq 2 then $
                 if linepars.fluxpk[haind,1] gt $
                 linepars.fluxpk[haind,0] then $
                    i_hapk = 1
              he_lam = he_rest * (1d + struct.z.gas[i_hapk])
              he_b = linepars.sigma[haind,i_hapk] * sqrt(2d) / $
                     he_lam * 299792d
           endif
        endif

        if file_test(parout) then spawn,'rm '+parout

fit:


;       Find nearest spectrum that has already been fit for input
;       parameters.
        dref_new = sqrt((cols-cols_srt[i])^2d + (rows-rows_srt[i])^2d)
        dref_new_isrt   = sort(dref_new)
        dref_new_srt    = reform(dref_new[dref_new_isrt],nspec)
        donefit_new_srt = reform(donefit[dref_new_isrt],nspec)
        nabscomp_new_srt= reform(nabscomp[dref_new_isrt],nspec)
        cols_new_srt    = reform(cols[dref_new_isrt],nspec)
        rows_new_srt    = reform(rows[dref_new_isrt],nspec)
        foundnearest = 0
        j = 0
        while ~ foundnearest AND j lt nspec do begin
           if donefit_new_srt[j] eq 1 AND $
              nabscomp_new_srt[j] eq nabs then foundnearest++
           j++
        endwhile
        j--
        if ~ foundnearest then begin
           useref = 1
           j = where(cols_new_srt eq refx AND rows_new_srt eq refy)
        endif else useref = 0

        if (nabscomp[cols_srt[i]-1,rows_srt[i]-1] eq 2 OR $
            badcomp ne 0) AND $
           nabs eq 1 then print,'',format='(A18,$)'
        print,'  Ref: [',cols_new_srt[j],',',rows_new_srt[j],$
              '].  Spec. ',idofit,' of ',ndofit,'.  Nabs =',nabs,$
              ', Nem =',nem,'.',$
              format='(A,I03,A,I03,A,I0,A,I0,A,I0,A,I0,A)'

;       Get output parameters from nearest fit        
        parfile_nearest = fitdir+gal+'_'+$
                          string(cols_new_srt[j],'_',rows_new_srt[j],$
                                 format='(I04,A0,I04)')+$
                          '_nad_parout'
        if nabs eq 2 then parfile_nearest += '.2comp'
        parfile_nearest+='.dat'

        gmos_readnadpars,parfile_nearest,abspars,empars,opars
        
        openw,lun,parin,/get_lun
        printf,lun,nabs,nem,fitran[0],fitran[1],$
               format='(I8,I8,D8.1,D8.1)'
        fparin = '(D8.2,D8.2,D8.2,I8)'
        printf,lun,0.01,abspars[0,0],1.0,1,format=fparin
        if abspars[1,0] gt 4.9 then $
           printf,lun,0.01,1.0,taumax,1,format=fparin $
        else $
           printf,lun,0.01,abspars[1,0],taumax,1,format=fparin
        printf,lun,abspars[2,0]-10d,abspars[2,0],abspars[2,0]+10d,1,$
               format=fparin
        if keyword_set(sigfix) then $
           printf,lun,60,sigfix*sqrt(2d),bmax,0,format=fparin $
        else $
           printf,lun,60,abspars[3,0],bmax,1,format=fparin
        if nabs eq 2 then begin
           printf,lun,0.01,abspars[0,1],1.0,1,format=fparin
           if abspars[1,1] gt 4.9 then $
              printf,lun,0.01,1.0,taumax,1,format=fparin $
           else $
              printf,lun,0.01,abspars[1,1],taumax,1,format=fparin
           printf,lun,abspars[2,1]-10d,abspars[2,1],abspars[2,1]+10d,1,$
                  format=fparin
           printf,lun,60,abspars[3,1],bmax,1,format=fparin
        endif
        if nem eq 1 then begin
           printf,lun,0,0.1,2,1,format=fparin
           printf,lun,0,he_lam,0,0,format=fparin
           printf,lun,0,he_b,0,0,format=fparin
        endif
        free_lun,lun

        file_copy,srcdir+'lmfit.h',srcdir+'lmfit.h.BAK',/over
        lmfith = strarr(15)
        dummy = ''
        openr,lun,srcdir+'lmfit.h',/get_lun
        for j=0,11 do begin
           readf,lun,dummy
           lmfith[j] = dummy        
        endfor
        lmfith[12] = '#define FILE_PARI "'+parin+'"'
        lmfith[13] = '#define FILE_OUT  "'+parout+'"'
        lmfith[14] = '#define FILE_DAT  "'+specin+'"'
        free_lun,lun

        openw,lun,srcdir+'lmfit.h',/get_lun
        for j=0,14 do printf,lun,lmfith[j]
        free_lun,lun
     
        spawn,'gcc -o '+fitdir+'lmfit '+srcdir+'lmfit.c -I'+$
              srcdir+'/lib -lm'
        spawn,fitdir+'lmfit'

        if file_test(parout) then begin

           if nabs eq 2 then begin
              gmos_plotnadfit,specin,parout,specout,struct.z
              badcomp = gmos_checknad(specin,parout,fitran)
              if badcomp gt 0 then begin
                 nabscomp[cols_srt[i]-1,rows_srt[i]-1] = 1
                 print,'','Rejecting 2-component fit.',format='(A20,A)'
              endif
              spawn,'mv '+spectot+'_nad_parout.dat'+' '+spectot+$
                    '_nad_parout.2comp.dat'
              spawn,'mv '+specout+'.jpg '+specout+'.2comp.jpg'
              nabs = 1
              goto,fit
           endif else begin
              gmos_plotnadfit,specin,parout,specout,struct.z
              donefit[cols_srt[i]-1,rows_srt[i]-1] = 1
              cols_done = [cols_done,cols_srt[i]]
              rows_done = [rows_done,rows_srt[i]]
              nabs_done = [nabs_done,$
                           nabscomp[cols_srt[i]-1,rows_srt[i]-1]]
           endelse

        endif

        idofit++

     endif

  endfor

  gmos_printnadpars,gal,cols_done,rows_done,nabs_done,$
                    fitdir,gal+'.nad',struct.z

end
