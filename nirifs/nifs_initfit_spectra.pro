function nifs_initfit_spectra,gal,bin,cpass=cpass

  if bin eq 3 then binlab = '6x6' else binlab = '2x2'
  ;; startempfile='/Users/drupke/src/idl/uhspecfit/stellar_spectra/'+$
  ;;              'wallace_hinkle/wallacehinkle.xdr'
  startempfile=''
  
  fcnlinefit='manygauss'
  fcninitpar='nifs_initparinfo'
  fcnstrlines='nifs_plotstronglines'

  vdisp=60d ; Corresponds to 4 A at 2um
  sigguess=0
  argstrlines=0
     
  if gal eq 'mrk231' then begin
     outdir = '/Users/drupke/winds/ao/specfits/'+gal+'/'+binlab+'/'
     rootindir = '/Users/drupke/winds/ao/red/'+gal+'/'
     rootindir += 'science/merged/'
     if bin eq 3 then begin
        dx = 15
        dy = 15
        cx = 8d
        cy = 8d
        pscale = 0.1035d*3d        ; in arcseconds / pixel
     endif else if bin eq 99 then begin
        dx = 1
        dy = 1
        cx = 1d
        cy = 1d
        pscale = 0.1035d        ; in arcseconds / pixel
     endif else begin
        dx = 45
        dy = 45
        cx = 23d
        cy = 23d
        pscale = 0.1035d        ; in arcseconds / pixel
     endelse
     doubleline = strarr(dx,dy)
     doubleline[*,*] = 'h2'
     if ~keyword_set(cpass) then begin
        ncomp = dblarr(dx,dy)+2
;; ;       Outer radius for fitting two H_2 components, in arcseconds
;;         cols = rebin(dindgen(dx)+1,dx,dy)
;;         rows = rebin(reform(dindgen(dy)+1,1,dy),dx,dy)
;;         windthresh = 1d
;;         dnuc = sqrt((cols-cx)^2d + (rows-cy)^2d)
;;         iwind = where(dnuc le windthresh/pscale)
;;         ncomp[iwind] = 3
;; ;       Region in which to fit two Paa components
;;         if bin eq 3 then begin
;; ;          Have not checked these #s ...
;;            col_hii = [4]
;;            row_hii = [6]
;;            ncomp[6,6] = -2
;;         endif else begin
;;            col_hii = [11,11,11,12,12,12,13,13,13]
;;            row_hii = [20,21,22,20,21,22,20,21,22]
;;         endelse
;;         col_hii--
;;         row_hii--
;;         for i=0,n_elements(col_hii)-1 do begin
;;            ncomp[col_hii[i],row_hii[i]] = 3
;;            doubleline[col_hii[i],row_hii[i]] = 'paa'
;;         endfor
        zinit = dblarr(dx,dy,3)+0.0422d
     endif else begin
        restore,file='/Users/drupke/winds/ao/red/mrk231/'+$
                'science/merged/mrk231_initcpass.xdr'
        ncomp = initcpass.ncomp
        zinit = initcpass.zinit
     endelse
     if bin eq 999 then begin
        infile=rootindir+gal+'all'+binlab+'k_h2sum_noct_2comp_ord1cfit_c1.fits'
        ;; infile=rootindir+gal+'all'+binlab+'k_h2sum_noct_2comp_ord2cfit_c1.fits'
        ctoutfile=rootindir+gal+'all'+binlab+'k_h2sum_c1_ct.fits'
        qsotmpfile = 0
        fcncontfit='fit_no_continuum'
        plotstelfit='plotstelfit'
        ncomp = -3
        argscontfit = {addnorm: 2463.38}
        ;; argscontfit = {addnorm: 2289d}
        qsocntargs = 0
     endif else begin
        infile=rootindir+gal+'all'+binlab+'k.fits'
        ctoutfile=rootindir+gal+'all'+binlab+'k_ct.fits' 
        qsotmpfile = rootindir+gal+'nuck_corr.fits'
        fcncontfit='nifs_fit_qso_continuum'
        plotstelfit='nifs_plotqsofit'
        fitord=0
        qsoonly=1
        qsoord=3
        expterms=4
        argscontfit = {$
                      fitord:fitord,$
                      qsoord:qsoord,$
                      expterms:expterms,$
                      qsotmp:qsotmpfile,$
                      qsoonly:qsoonly $
                      }
        qsocntargs = {$
                     fitord:fitord,$
                     qsoord:qsoord,$
                     expterms:expterms,$
                     qsotmp:qsotmpfile, $
                     qsoonly:qsoonly $
                     }
     endelse
     fitran_rest_lo=dblarr(dx,dy)+18500
     fitran_rest_hi=dblarr(dx,dy)+21900
     fitblr=dblarr(dx,dy)
;    Sigma = 1.6 A corresponds to R = 5290 (from NIFS website) at 2um
;    Sigma = 30 A corresponds to 1060 km/s FWHM at 2um
;    Sigma = 85 A corresponds to 3000 km/s FWHM at 2um
     siglim=[2d,50d]
     argsinitpar = {siglim:siglim}
;    Switch to negative values for 'NIFS_CHECKCOMP' so that pegged
;    values are not allowed.
     siglim=-siglim
     argslinelist = replicate({hei206:1},dx,dy)
     argscheckcomp = {paasig: 3d,$
                      dblsigh2: 3d,$
                      singsigh2: 100d,$
                      siglim: siglim,$
                      use_h2_10_s2: [0,0]}
     outlines = ['H2_10_S3','H2_10_S2','H2_10_S1',$
                 'HeI206',$
                 'Paa','Brd','Brg']

  endif
     
  if gal eq 'f08572nw' then begin

     pipever = 3.2
     outdir = '/Users/drupke/winds/ao/specfits/'+gal+'/'+$
              binlab+'/v'+string(pipever,format='(D0.1)')+'/'
     rootindir = '/Users/drupke/winds/ao/red/'+gal+$
                 '/v'+string(pipever,format='(D0.1)')+'/'
     if pipever eq 2.3 then begin
        ;; rootindir = '/Users/drupke/winds/ao/red/'+gal+'/noctcl/'
        rootindir += 'meanclip2p35/'
        dx = 42
        dy = 14
        cx = 20d
        cy = 7d
     endif else begin
        dx = 41
        dy = 14
        cx = 23d
        cy = 7d
     endelse
;    for testing nuclear spaxel fit
     ;; rootindir = '/Users/drupke/winds/ao/red/'+gal+'/nucleus/'
     infile = rootindir+'s130131_a016004_tlc_Kbb_035_bin2.fits'
     ;; infile = rootindir+'s130131_a016004_ntlc_Kbb_035_bin2.fits' 
     ctoutfile=rootindir+gal+'_ct.fits'
     fcncontfit='nifs_fit_flat_continuum'
     ;; fcncontfit='fit_continuum'
     plotstelfit='plotstelfit'
     pscale = 0.07d           ; in arcseconds / pixel
     ncomp = dblarr(dx,dy)+2
     fitran_rest_lo=dblarr(dx,dy)+18500
     fitran_rest_hi=dblarr(dx,dy)+22000
     fitord=3
     fitblr = dblarr(dx,dy)
     if pipever eq 2.3 then begin
        ncomp[*,0] = 0
        ncomp[0:12,*] = 1
        ncomp[25:dx-1,*] = 1
; Components adjusted to add or remove systemic H_2
        ncomp[13,2:4] = -1
        ncomp[14,7] = -2
; Components adjusted to include or remove OF
        ncomp[14:20,2:7] = 3
        ncomp[14,3] = -1
        ncomp[15,2] = -2
        ;; ncomp[15,3] = -2 ; a bit iffy, but let's include it
        ncomp[15,5] = -3        ; a bit iffy, but let's include it
        ncomp[16,6] = -3        ; a bit iffy, but let's include it
        ncomp[17,5] = -3
        ncomp[17,6] = -3
        ;; ncomp[18,2] = -2
        ;; ncomp[18,3] = -3 ; tried this but just fits a narrow line
        ncomp[18,5] = -3
        ncomp[18,6] = -3
        ncomp[19,4] = -3
        ncomp[19,5] = -3
        ncomp[19,2] = -2
        ncomp[20,2] = -2
        fitran_rest_hi[19,2] = 22600/1.059
        fitran_rest_hi[19:20,1] = 22400/1.059
        fitran_rest_hi[19:20,0] = 22200/1.059
        fitblr[19,6] = 1
        zsys = 0.059d
     endif else begin
        ncomp[*,0] = 0
; Components adjusted to add or remove systemic H_2
        ncomp[*,1] = 1
        ncomp[*,12] = 1 
        ncomp[0:15,*] = 1
        ncomp[16:17,*] = 1
        ncomp[20,3] = -2
        ncomp[26:27,11] = 1
        ncomp[27,7] = -2
        ncomp[28,3:4] = 1
        ncomp[29:dx-1,*] = 1
; Components adjusted to include or remove OF
        ncomp[22:27,2:7] = 3
        ncomp[22,2] = 2
        ncomp[22,5] = -3
        ncomp[23,6] = -3
        ncomp[25,6] = -3        ; a bit iffy, but let's include it
        ncomp[26,2] = 1
        ;; ncomp[26,3] = 2      ; a bit iffy, but let's include it
        ncomp[26,5] = -3        ; a bit iffy, but let's include it
        ncomp[27,3] = 1
        fitblr[22,6] = 1
        zsys = 0.058d
     endelse
     zinit = dblarr(dx,dy,3)+zsys
     zinit[*,*,2] = zsys-0.003d
     doubleline = strarr(dx,dy)
     doubleline[*,*] = 'h2'
; Region in which to fit two Paa components
     ;; col_hii = [20]
     ;; row_hii = [7]
     ;; col_hii--
     ;; row_hii--
     ;; for i=0,n_elements(col_hii)-1 do begin
     ;;    ncomp[col_hii[i],row_hii[i]] = 3
     ;;    doubleline[col_hii[i],row_hii[i]] = 'paa'
     ;; endfor

     argscontfit = {fitord:fitord}
;    for testing nuclear spaxel fit
     ;; argscontfit = {fitord:fitord,norefit:1}
;    Sigma = 1.6 A corresponds to R = 5290 (from NIFS website) at 2um
;    Sigma = 30 A corresponds to 1060 km/s FWHM at 2um
;    Sigma = 85 A corresponds to 3000 km/s FWHM at 2um
     siglim=[3d,30d]
     argsinitpar = {siglim:siglim}
;    Switch to negative values for 'NIFS_CHECKCOMP' so that pegged
;    values are not allowed.
     siglim=-siglim
     ;; argslinelist = replicate({hei187:0,hei206:0,h2s4:0},dx,dy) ;
     argslinelist = replicate({hei187:0,hei206:0},dx,dy)
; Outer radius for fitting two H_2 components, in arcseconds
     cols = rebin(dindgen(dx)+1,dx,dy)
     rows = rebin(reform(dindgen(dy)+1,1,dy),dx,dy)
     heithresh = 0.3d
     dnuc = sqrt((rows-cy)^2d + (cols-cx)^2d)
     ihei = where(dnuc le heithresh/pscale)
     ;; argslinelist[ihei] = {hei187:1,hei206:1,h2s4:1}
     ;; outlines = ['H2_10_S3','H2_10_S2','H2_10_S1',$
     ;;             'Paa','Brd','Brg','HeI187','HeI206','H2_10_S4']
     argslinelist[ihei] = {hei187:1,hei206:1}
     argscheckcomp = {paasig: 3d,$
                      dblsigh2: 3d,$
                      singsigh2: 5d,$
                      siglim: siglim,$
                      use_h2_10_s2: [1,0]}
     outlines = ['H2_10_S3','H2_10_S2','H2_10_S1',$
                 'Paa','Brd','Brg','HeI187','HeI206']
     qsocntargs=0

  endif

; Make sure data exists  
  if ~ file_test(infile) then begin
     print,"NIFS_INITFIT_SPECTRA: Data cube not found."
     stop
  endif

  init = {zinit            : zinit, $
          ncomp            : ncomp, $
          infile           : infile, $
          outdir           : outdir, $
          ctoutfile        : ctoutfile, $
          fcnlinefit       : fcnlinefit, $
          fcninitpar       : fcninitpar, $
          fcncontfit       : fcncontfit, $
          fcnstrlines      : fcnstrlines, $
          plotstelfit      : plotstelfit, $
          argsinitpar      : argsinitpar, $
          argscontfit      : argscontfit, $
          argstrlines      : argstrlines, $
          argslinelist     : argslinelist, $
          argscheckcomp     : argscheckcomp, $
          qsocntargs       : qsocntargs, $
          fitran_rest_lo   : fitran_rest_lo, $
          fitran_rest_hi   : fitran_rest_hi, $
          vdisp            : vdisp, $
          startempfile     : startempfile,$
          outlines         : outlines, $
          dx               : dx,$
          dy               : dy, $
          cx               : cx,$
          cy               : cy, $
          siglim           : siglim, $
          sigguess         : sigguess, $
          doubleline       : doubleline, $
          fitblr           : fitblr $
         }

  return,init

end
