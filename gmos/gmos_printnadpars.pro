;
; History
;  10jul22  DSNR  created
;

pro gmos_printnadpars,gal,cols,rows,nabs,fitdir,outfile,galz

  nad1_rest = 5895.92d
  nad2_rest = 5889.95d

  nspec = n_elements(rows)

  openw,lun,fitdir+outfile+'.dat',/get_lun
  printf,lun,'#Col','Row','Cmp','z','FWHM','tau','C_f','W_eq',$
         format='(A-4,2A4,5A10)'
  openw,lun2,fitdir+outfile+'.2comp.dat',/get_lun
  printf,lun2,'#Col','Row','Cmp','z','FWHM','tau','C_f','W_eq',$
         format='(A-4,2A4,5A10)'

  for i=0,nspec-1 do begin

     spec = string(cols[i],'_',rows[i],format='(I04,A0,I04)')
     spectot= fitdir+gal+'_'+spec
     specin = spectot+'_nad_spec.dat'
     parout = spectot+'_nad_parout.dat'
     gmos_readnadpars,parout,abspars,empars,opars
     weq = gmos_compweq(specin,parout,galz)
     printf,lun,cols[i],rows[i],1,abspars[2,0]/nad1_rest-1d,$
            abspars[3,0]*2d*sqrt(alog(2d)),abspars[1,0],abspars[0,0],$
            weq,format='(3I4,D10.6,D10.2,D10.4,D10.4,D10.4)'

     if nabs[i] eq 2 then begin

        specin = spectot+'_nad_spec.dat'
        parout2 = spectot+'_nad_parout.2comp.dat'
        gmos_readnadpars,parout2,abspars,empars,opars
        weq = gmos_compweq(specin,parout2,galz)
        for j=0,opars.nabs-1 do $
           printf,lun2,cols[i],rows[i],j+1,abspars[2,j]/nad1_rest-1d,$
                  abspars[3,j]*2d*sqrt(alog(2d)),abspars[1,j],$
                  abspars[0,j],weq,$
                  format='(3I4,D10.6,D10.2,D10.4,D10.4,D10.4)'
        
     endif

  endfor
  
  free_lun,lun
  free_lun,lun2

end
