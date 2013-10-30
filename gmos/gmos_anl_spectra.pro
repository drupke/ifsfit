pro gmos_anl_spectra,showplot=showplot
;
; History
;  09may13  DSNR  created
;  09jun07  DSNR  added multiple components
;

  if keyword_set(showplot) then zbuf=0 else zbuf=1

  fwhmtosig = 2d*sqrt(2d*alog(2d))

; For GMOS spectra
;;   rootdir = '/Users/drupke/winds/gmos/f10565/'
;;   outpre = 'f10565p2448'
;;   outfile = rootdir+'specfit/'+outpre
;;   startrow=1
;;   nrow=11
;;   startcol=1
;;   ncol=13
; For KPNO longslit spectra
  rootdir = '/Users/drupke/winds/gmos/f10565/kpnolong/'
  outpre = 'f10565p2448'
  outfile = rootdir+'specfit/'+outpre
  startrow=1
  nrow=9
  startcol=1
  ncol=1
;;   startrow=6
;;   nrow=1
;;   startcol=9
;;   ncol=1
; for corner spectra
;;   startrow=101
;;   nrow=2
;;   startcol=101
;;   ncol=2

; Get linelist
  linelist = gmos_initlinelist()
; which lines to output
  outlines=['[OI]6300','Halpha','[NII]6548']
  gmos_printlinepars,[''],[0],[0],outfile+'.lines.dat','',/init,$
                     whichlines=outlines

; Open output file for fit results
  openw,fitunit,outfile+'.fit.dat',/get_lun
  printf,fitunit,'# Col/Row','Rchi2','Niter','fwhm','z',$
         format='(A-10,A6,A6,A8,A9)'

  for i=startrow,startrow+nrow-1 do begin

     for j=startcol,startcol+ncol-1 do begin
      
        print,'Processing spectrum ',i,', ',j,format='(A,I0,A,I0)'

;       GMOS
        lab = string(i,'_',j,format='(I03,A,I03)')
;       KPNO
        lab = string(i,format='(I0)')
        infile = rootdir+'specfit/'+outpre+"_"+lab+'.genx'

        if (i eq 1 AND j eq 1) then append=0 else append=1

;       Restore fit from binary file
        restgen,struct=struct,file=infile

        if n_elements(struct.param) gt 1 then $
           ncomp = struct.param[1] $
        else $
           ncomp = 0

        if ncomp then begin
;       Separate best-fit parameters
           linepars = sepfitpars(struct.param,struct.perror)
           nlines = n_elements(linepars.flux)
        endif

;       Plot results
        if ~ keyword_set(noplot) then begin
           plotstelfit,struct,outfile+'_'+lab+'_stel',zbuf=zbuf
           gmos_plotstronglines,struct,outfile+'_'+lab+'_strg',zbuf=zbuf
        endif

;       Print fit parameters to a text file
        printf,fitunit,struct.obj_id,$
               struct.redchisq,$
               struct.niter,$
               format='(A-10,D6.2,I6.0,$)'        
        if ncomp gt 0 then $
           printf,fitunit,fwhmtosig*linepars.sigma[0,0],$
                  struct.z[0],format='(D8.2,D9.5,$)' $
        else printf,fitunit,-1,-1,format='(D8.2,D9.5,$)'
        if ncomp eq 2 then $
           printf,fitunit,fwhmtosig*linepars.sigma[0,1],$
                  struct.z[1],format='(D8.2,D9.5,$)' $
        else printf,fitunit,-1,-1,format='(D8.2,D9.5,$)'
        printf,fitunit,''

;       Print line fluxes to a text file
        if ncomp then $
           gmos_printlinepars,struct.linelabel,linepars.flux,linepars.fluxerr,$
                              outfile+'.lines.dat',lab,/append,$
                              whichlines=outlines

     endfor

  endfor

  free_lun,fitunit

end
