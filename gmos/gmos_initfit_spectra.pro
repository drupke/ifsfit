function gmos_initfit_spectra,gal,bin=bin

  if ~keyword_set(bin) then bin=9
  bin = double(bin)

; Default infile
  infile = ''
; Directory for input files
  rootindir = '/Users/drupke/winds/gmos/'
; IDL save file with stellar template structure
  startempfile = '/Users/drupke/src/idl/uhspecfit/stellar_models/'+$
                 'gonzalezdelgado/SSPGeneva_z020.sav'




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Sky line measurements
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  zinit = dblarr(1496)
  ncomp = intarr(1496)+1
  fcninitpar='gmos_initparinfo'
  fcnlinefit='gmos_manygauss'
  fcncontfit='gmos_fit_flat_continuum'
  plotstelfit='plotstelfit'
  fcnstrlines='gmos_plotskylines'
  plotstelfit_args = {zbuf:1}
  argscontfit = {fitord:2}
  argstrlines = 0
  argslinelist = 0
  vdisp=75d
  fluxsigthresh = 3
  ctoutfile=0
  qsoutfile=0
  qsoutfile_ha=0
  qsocntargs=0
  dx = 0
  dy = 0
  cx = 0
  cy = 0
  nad1thresh = 0
  nad2thresh = 0
  nadxref = 0
  nadyref = 0
  nadfitran = 0
  nadfitord = 0
  nadrestcomp = 0
  fitran_rest=[6285d,6380d]
  siglim = [0.699d,0.701d]
  sigguess = dblarr(7)+0.7d
  argsinitpar = {siglim:siglim}
  outlines = ['[OI]6300']

; Mrk 231
  if gal eq 'mrk231a' then begin
     infile=rootindir+'mrk231/red/ctexrdat_3a.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/mrk231/o1/expa/'
  endif
  if gal eq 'mrk231b' then begin
     infile=rootindir+'mrk231/red/ctexrdat_3b.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/mrk231/o1/expb/'
  endif
  if gal eq 'mrk231c' then begin
     infile=rootindir+'mrk231/red/ctexrdat_3c.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/mrk231/o1/expc/'
  endif
  if gal eq 'mrk231d' then begin
     infile=rootindir+'mrk231/red/ctexrdat_3d.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/mrk231/o1/expd/'
  endif
  if gal eq 'mrk231e' then begin
     infile=rootindir+'mrk231/red/ctexrdat_3e.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/mrk231/o1/expe/'
  endif
; Mrk 273
  if gal eq 'mrk273a' then begin
     infile=rootindir+'mrk273/red/ctexrdat_2a.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/mrk273/o1/expa/'
  endif
  if gal eq 'mrk273b' then begin
     infile=rootindir+'mrk273/red/ctexrdat_2b.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/mrk273/o1/expb/'
  endif
  if gal eq 'mrk273c' then begin
     infile=rootindir+'mrk273/red/ctexrdat_2c.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/mrk273/o1/expc/'
  endif
  if gal eq 'mrk273d' then begin
     infile=rootindir+'mrk273/red/ctexrdat_2d.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/mrk273/o1/expd/'
  endif
; F08572+3915
  if gal eq 'f08572a' then begin
     infile=rootindir+'f08572/red/ctexrdat_1a.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f08572/o1/expa/'
  endif
  if gal eq 'f08572b' then begin
     infile=rootindir+'f08572/red/ctexrdat_1b.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f08572/o1/expb/'
  endif
  if gal eq 'f08572c' then begin
     infile=rootindir+'f08572/red/ctexrdat_1c.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f08572/o1/expc/'
  endif
  if gal eq 'f08572d' then begin
     infile=rootindir+'f08572/red/ctexrdat_2a.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f08572/o1/expd/'
  endif
  if gal eq 'f08572e' then begin
     infile=rootindir+'f08572/red/ctexrdat_2b.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f08572/o1/expe/'
  endif
  if gal eq 'f08572f' then begin
     infile=rootindir+'f08572/red/ctexrdat_3a.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f08572/o1/expf/'
  endif
  if gal eq 'f08572g' then begin
     infile=rootindir+'f08572/red/ctexrdat_4a.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f08572/o1/expg/'
  endif
  if gal eq 'f08572h' then begin
     infile=rootindir+'f08572/red/ctexrdat_4b.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f08572/o1/exph/'
  endif
  if gal eq 'f08572i' then begin
     infile=rootindir+'f08572/red/ctexrdat_4c.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f08572/o1/expi/'
  endif
; VV705
  if gal eq 'vv705a' then begin
     infile=rootindir+'vv705/red/ctexrdat_1a.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/vv705/o1/expa/'
  endif
  if gal eq 'vv705b' then begin
     infile=rootindir+'vv705/red/ctexrdat_1b.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/vv705/o1/expb/'
  endif
  if gal eq 'vv705c' then begin
     infile=rootindir+'vv705/red/ctexrdat_1c.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/vv705/o1/expc/'
  endif
  if gal eq 'vv705d' then begin
     infile=rootindir+'vv705/red/ctexrdat_2a.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/vv705/o1/expd/'
  endif
  if gal eq 'vv705e' then begin
     infile=rootindir+'vv705/red/ctexrdat_2b.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/vv705/o1/expe/'
  endif
  if gal eq 'vv705f' then begin
     infile=rootindir+'vv705/red/ctexrdat_2c.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/vv705/o1/expf/'
  endif
  if gal eq 'vv705g' then begin
     infile=rootindir+'vv705/red/ctexrdat_2d.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/vv705/o1/expg/'
  endif
  if gal eq 'vv705h' then begin
     infile=rootindir+'vv705/red/ctexrdat_3a.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/vv705/o1/exph/'
  endif
  if gal eq 'vv705i' then begin
     infile=rootindir+'vv705/red/ctexrdat_3b.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/vv705/o1/expi/'
  endif
; F17207
  if gal eq 'f17207a' then begin
     zinit = dblarr(1496)+0.0005
     infile=rootindir+'f17207/red/ctexrdat_1a.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f17207/nad/expa/'
     argslinelist = {addnad:1}
     fitran_rest=[5870d,5920d]
     siglim = [0.699d,0.701d]
     sigguess = dblarr(9)+0.7d
     argsinitpar = {siglim:siglim}
     outlines = ['NaD2','NaD1']
  endif
  if gal eq 'f17207b' then begin
     infile=rootindir+'f17207/red/ctexrdat_2a.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f17207/nad/expb/'
     argslinelist = {addnad:1}
     fitran_rest=[5870d,5920d]
     siglim = [0.699d,0.701d]
     sigguess = dblarr(9)+0.7d
     argsinitpar = {siglim:siglim}
     outlines = ['NaD2','NaD1']
  endif
  if gal eq 'f17207c' then begin
     infile=rootindir+'f17207/red/ctexrdat_3a.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f17207/nad/expc/'
     argslinelist = {addnad:1}
     fitran_rest=[5870d,5920d]
     siglim = [0.699d,0.701d]
     sigguess = dblarr(9)+0.7d
     argsinitpar = {siglim:siglim}
     outlines = ['NaD2','NaD1']
  endif
  if gal eq 'f17207d' then begin
     infile=rootindir+'f17207/red/ctexrdat_3b.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f17207/nad/expd/'
     argslinelist = {addnad:1}
     fitran_rest=[5870d,5920d]
     siglim = [0.699d,0.701d]
     sigguess = dblarr(9)+0.7d
     argsinitpar = {siglim:siglim}
     outlines = ['NaD2','NaD1']
  endif
  if gal eq 'f17207e' then begin
     infile=rootindir+'f17207/red/ctexrdat_3c.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f17207/nad/expe/'
     argslinelist = {addnad:1}
     fitran_rest=[5870d,5920d]
     siglim = [0.699d,0.701d]
     sigguess = dblarr(9)+0.7d
     argsinitpar = {siglim:siglim}
     outlines = ['NaD2','NaD1']
  endif
  if gal eq 'f17207f' then begin
     infile=rootindir+'f17207/red/ctexrdat_3d.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f17207/nad/expf/'
     argslinelist = {addnad:1}
     fitran_rest=[5870d,5920d]
     siglim = [0.699d,0.701d]
     sigguess = dblarr(9)+0.7d
     argsinitpar = {siglim:siglim}
     outlines = ['NaD2','NaD1']
  endif
  if gal eq 'f17207g' then begin
     zinit = dblarr(1496)+0.00035
     infile=rootindir+'f17207/red/ctexrdat_4a.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f17207/nad/expg/'
     argslinelist = {addnad:1}
     fitran_rest=[5870d,5920d]
     siglim = [0.699d,0.701d]
     sigguess = dblarr(9)+0.7d
     argsinitpar = {siglim:siglim}
     outlines = ['NaD2','NaD1']
  endif
  if gal eq 'f17207h' then begin
     zinit = dblarr(1496)+0.00035
     infile=rootindir+'f17207/red/ctexrdat_4b.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f17207/nad/exph/'
     argslinelist = {addnad:1}
     fitran_rest=[5870d,5920d]
     siglim = [0.699d,0.701d]
     sigguess = dblarr(9)+0.7d
     argsinitpar = {siglim:siglim}
     outlines = ['NaD2','NaD1']
  endif
  if gal eq 'f17207i' then begin
     zinit = dblarr(1496)+0.00035
     infile=rootindir+'f17207/red/ctexrdat_4c.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f17207/nad/expi/'
     argslinelist = {addnad:1}
     fitran_rest=[5870d,5920d]
     siglim = [0.699d,0.701d]
     sigguess = dblarr(9)+0.7d
     argsinitpar = {siglim:siglim}
     outlines = ['NaD2','NaD1']
  endif
; F10565
  if gal eq 'f10565a' then begin
     infile=rootindir+'f10565/red/ctexrdat_1a.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f10565/nad/expa/'
     argslinelist = {addnad:1}
     fitran_rest=[5870d,5920d]
     siglim = [0.699d,0.701d]
     sigguess = dblarr(9)+0.7d
     argsinitpar = {siglim:siglim}
     outlines = ['NaD2','NaD1']
  endif
  if gal eq 'f10565b' then begin
     infile=rootindir+'f10565/red/ctexrdat_1b.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f10565/nad/expb/'
     argslinelist = {addnad:1}
     fitran_rest=[5870d,5920d]
     siglim = [0.699d,0.701d]
     sigguess = dblarr(9)+0.7d
     argsinitpar = {siglim:siglim}
     outlines = ['NaD2','NaD1']
  endif
  if gal eq 'f10565c' then begin
     infile=rootindir+'f10565/red/ctexrdat_1c.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f10565/nad/expc/'
     argslinelist = {addnad:1}
     fitran_rest=[5870d,5920d]
     siglim = [0.699d,0.701d]
     sigguess = dblarr(9)+0.7d
     argsinitpar = {siglim:siglim}
     outlines = ['NaD2','NaD1']
  endif
  if gal eq 'f10565d' then begin
     infile=rootindir+'f10565/red/ctexrdat_1d.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f10565/nad/expd/'
     argslinelist = {addnad:1}
     fitran_rest=[5870d,5920d]
     siglim = [0.699d,0.701d]
     sigguess = dblarr(9)+0.7d
     argsinitpar = {siglim:siglim}
     outlines = ['NaD2','NaD1']
  endif
  if gal eq 'f10565e' then begin
     infile=rootindir+'f10565/red/ctexrdat_1e.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f10565/nad/expe/'
     argslinelist = {addnad:1}
     fitran_rest=[5870d,5920d]
     siglim = [0.699d,0.701d]
     sigguess = dblarr(9)+0.7d
     argsinitpar = {siglim:siglim}
     outlines = ['NaD2','NaD1']
  endif
  if gal eq 'f10565f' then begin
     infile=rootindir+'f10565/red/ctexrdat_1f.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f10565/nad/expf/'
     argslinelist = {addnad:1}
     fitran_rest=[5870d,5920d]
     siglim = [0.699d,0.701d]
     sigguess = dblarr(9)+0.7d
     argsinitpar = {siglim:siglim}
     outlines = ['NaD2','NaD1']
  endif
  if gal eq 'f10565g' then begin
     infile=rootindir+'f10565/red/ctexrdat_1g.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/f10565/nad/expg/'
     argslinelist = {addnad:1}
     fitran_rest=[5870d,5920d]
     siglim = [0.699d,0.701d]
     sigguess = dblarr(9)+0.7d
     argsinitpar = {siglim:siglim}
     outlines = ['NaD2','NaD1']
  endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Mrk 231
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



  if gal eq 'mrk231' then begin
;    Binning = 1x1
     if bin eq 1 then begin
        dx = 71
        dy = 82
        cx = 36d
        cy = 41d
        dthresh=0.5d            ; in arcseconds
        pscale = bin/10d           ; in arcseconds / pixel
        outstr = 'rb'+string(bin,format='(I0)')
        ncomp = dblarr(dx,dy)+2
        rows = rebin(dindgen(dx)+1,dx,dy)
        cols = rebin(reform(dindgen(dy)+1,1,dy),dx,dy)
        dnuc = sqrt((rows-cx)^2d + (cols-cy)^2d)
        inuc = where(dnuc le dthresh/pscale)
        ncomp[inuc] = 0
        fluxsigthresh = 7d
        nad1thresh = 5d
        nad2thresh = 10d
        nadxref = 7
        nadyref = 13
        nadfitran = [6100d,6150d]
     endif else if bin eq 3 then begin
;    Binning = 3x3
        dx = 21
        dy = 25
        cx = 11d
        cy = 13d
        dthresh=0.5d            ; in arcseconds
        pscale = bin/10d           ; in arcseconds / pixel
        outstr = 'rb'+string(bin,format='(I0)')
        ncomp = dblarr(dx,dy)+2
        rows = rebin(dindgen(dx)+1,dx,dy)
        cols = rebin(reform(dindgen(dy)+1,1,dy),dx,dy)
        dnuc = sqrt((rows-cx)^2d + (cols-cy)^2d)
        inuc = where(dnuc le dthresh/pscale)
        ncomp[inuc] = 0
        ncomp[9,17] = -1
        ncomp[10,17] = -1
        ncomp[10,18] = -1
        ncomp[10,19] = -1
        ncomp[14,13] = 3
        ncomp[14,14] = 3
        ncomp[13,14] = 3
        fluxsigthresh = 7d
        nad1thresh = 4d
        nad2thresh = 10d
        nadxref = 7
        nadyref = 13
        nadfitran = [6100d,6150d]
     endif else if bin eq 5 then begin
;    Binning = 5x5
        dx = 13
        dy = 15
        cx = 7d
        cy = 8d
        dthresh=0.5d            ; in arcseconds
        pscale = bin/10d           ; in arcseconds / pixel
        outstr = 'rb'+string(bin,format='(I0)')
        ncomp = dblarr(dx,dy)+2
        rows = rebin(dindgen(dx)+1,dx,dy)
        cols = rebin(reform(dindgen(dy)+1,1,dy),dx,dy)
        dnuc = sqrt((rows-cx)^2d + (cols-cy)^2d)
        inuc = where(dnuc le dthresh/pscale)
        ncomp[inuc] = 0
        fluxsigthresh = 7d
        nad1thresh = 5d
        nad2thresh = 10d
        nadxref = 4
        nadyref = 8
        nadfitran = [6100d,6150d]
     endif else if bin eq 9 then begin
;    Binning = 9x9
        dx = 7
        dy = 9
        cx = 4d
        cy = 5d
        dthresh=0.3d            ; in arcseconds
        pscale = bin/10d           ; in arcseconds / pixel
        outstr = 'rb'+string(bin,format='(I0)')
        ncomp = dblarr(dx,dy)+2
        rows = rebin(dindgen(dx)+1,dx,dy)
        cols = rebin(reform(dindgen(dy)+1,1,dy),dx,dy)
        dnuc = sqrt((rows-cx)^2d + (cols-cy)^2d)
        inuc = where(dnuc le dthresh/pscale)
        ncomp[inuc] = 0
        fluxsigthresh = 7d
        nad1thresh = 5d
        nad2thresh = 10d
        nadxref = 3
        nadyref = 5
        nadfitran = [6100d,6150d]
     endif else begin
        print,'GMOS_INITFIT_SPECTRA: This binning not set up.'
        stop
     endelse

     vdisp=75d
     outdir = '/Users/drupke/winds/gmos/specfits/'+gal+'/'+outstr+'/'
     infile=rootindir+'cubes/'+gal+'/'+gal+outstr+'.fits'
     ctoutfile=rootindir+'cubes/'+gal+'/'+gal+outstr+'_noqso.fits'
     qsoutfile=rootindir+'cubes/'+gal+'/'+gal+outstr+'_qsonly.fits'
     qsoutfile_ha=rootindir+'cubes/'+gal+'/'+gal+outstr+'_qsonly_ha.fits'
     qsotmpfile = '/Users/drupke/winds/gmos/red/'+gal+'/'+gal+'nuc_t.fits'
     zinit = dblarr(dx,dy,3)+0.0422d
     zinit[*,*,1:2] = 0.041d
     fitran_rest=[5750,6650]
;    9.69 A ~ 1000 km/s at 6840 A
     siglim = [0.8d,9.69d]
     sigguess = 0
     fcncontfit='gmos_fit_qso_continuum'
     fcnstrlines='gmos_plotstronglines'
     argstrlines=0
     fitord=3
     qsoord=3
     expterms=2
     argscontfit = {$
                   fitord:fitord,$
                   qsoord:qsoord,$
                   expterms:expterms,$
                   qsotmp:qsotmpfile,$
                   z:{star:0d,gas:0d} $
                   }
     plotstelfit='gmos_plotqsofit'
     qsocntargs = {$
                  fitord:fitord,$
                  qsoord:qsoord,$
                  expterms:expterms,$
                  qsotmp:qsotmpfile $
                  }
     argsinitpar = {siglim:siglim,polypars:[5,6770d,6920d]}
     nadfitord = 3
     nadrestcomp = 0
     argslinelist = 0
     outlines=0

  endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Mrk 273
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if gal eq 'mrk273' then begin
     if bin eq 3 then begin
;    Binning = 3x3
        dx = 20
        dy = 27
        cx = 9d
        cy = 13d
        pscale = bin/10d           ; in arcseconds / pixel
        outstr = 'rb'+string(bin,format='(I0)')
        ncomp = dblarr(dx,dy)+2
; 1-2 comp
        ncomp[0,3] = -2
        ncomp[0,5:6] = -2
        ncomp[0:3,22:26] = 1
        ncomp[2,0] = -2
        ncomp[6,13:14] = -2
        ncomp[10,1] = -2
        ncomp[13,16] = -2
        ncomp[15,0] = -2
        ncomp[19,0:9] = -2
; 3 comp

; discrete redshifted feature:
;  12,1-4
;  13,1-5
;  14,1-5
;  15,1-5
;  16,1
        ncomp[11,3] = 3
        ncomp[12,3:4] = 3
        ncomp[13,3:5] = 3
        ncomp[14,3:4] = 3

; discrete redshifted feature
;   7,8-9
;   8,7-10
;   9,7-10
;  10,9
        ncomp[6,7:8] = 3
        ncomp[7,6:9] = 3
        ncomp[8,6:9] = 3
        ncomp[9,8] = 3

;        ncomp[15,3:4] = 3
;        ncomp[16,3:4] = 3

; discrete blueshifted feature
        ncomp[13,22:24] = 3
; old
        ncomp[14,21:24] = 3
; newer try
        ; ncomp[14,21] = 3
        ; ncomp[14,22] = 4

; three dynamical components?
        ncomp[9,9:10] = 3
        ncomp[10,6:11] = 3
        ncomp[11,7:11] = 3
;        ncomp[12,7:10] = 3
;        ncomp[13,8:11] = 3

; combine two red components
        ncomp[4,10:11] = 3
        ncomp[4,14:16] = 3

; three dynamical components?
        ncomp[8,18:21] = 3
        ncomp[9,17:22] = 3
        ncomp[10,18:22] = 3
        ncomp[11,18:22] = 3
        ncomp[12,19:22] = 3
        ncomp[13,18:21] = 3
; old
        ncomp[14,18:20] = 3
; newer try
;        ncomp[14,15:20] = 3

; three dynamical components?
        ncomp[15,14:24] = 3
        ncomp[16,11:26] = 3
        ncomp[17,11:21] = 3
        ncomp[17,16] = -3
        ncomp[18,11:19] = 3
        ncomp[19,12:19] = 3

;?
        ncomp[13,25] = 3

; velocities
        zinit=dblarr(dx,dy,4)+0.037d
        zinit[0,3,1]=0.039d
        zinit[0,9,1]=0.039d
        zinit[1,10:12,1]=0.039d
        zinit[4,23:26,1]=0.040d
        zinit[5,20:24,1]=0.040d
        zinit[6,6:8,2]=0.039d
        zinit[7,9,2]=0.038d
        zinit[7,24:26,*]=0.0375d
        zinit[8,6:9,2]=0.038d
        zinit[8,15,1]=0.040d
        zinit[8,18:19,1]=0.033d        
        zinit[8,24:26,*]=0.0375d
        zinit[9,6:9,2]=0.038d
        zinit[9,17,0:2]=[0.036d,0.037d,0.040d]
        zinit[10,1,*]=0.038d
        zinit[10,4,1]=0.033d
        zinit[10,9,2]=0.033d
        zinit[10,18:22,2]=0.035d
        zinit[10,21,0:2]=[0.036d,0.037d,0.038d]
        zinit[11,0:2,0:1]=0.039d
        zinit[11,16:23,1]=0.040d
        zinit[11,17,*]=0.038d
        zinit[11,19,0:2]=[0.034d,0.036d,0.039d]
        zinit[11,26,0:1]=0.038d
        zinit[12,0:1,0:1]=0.039d
        zinit[12,4,0:2]=0.039d
        zinit[12,7:11,1]=0.035d
        zinit[12,16:23,1]=0.038d
        zinit[13,5,0:2]=[0.037d,0.037d,0.038d]
        zinit[13,8:11,0]=0.035d
        zinit[13,8:11,1]=0.038d
        zinit[13,16,0:1]=0.038d
        zinit[13,18:22,2] = 0.039d
        zinit[13,25,2] = 0.039d
; newer try
        ; zinit[14,15,0:2] = [0.038d,0.038d,0.039d]
        ; zinit[14,16,0:2] = [0.037d,0.038d,0.039d]
        zinit[14,18:19,1]=0.039d
; newer try
        ; zinit[14,21,0:2]=[0.037d,0.037d,0.038d]
        ; zinit[14,22,0:3]=[0.037d,0.038d,0.038d,0.039d]
        ; zinit[14,23,0:2]=[0.0385d,0.0385d,0.0385d]
        zinit[15,1,1]=0.038d
        zinit[15,1:6,1]=0.033d
        zinit[15,11,0]=0.035d
        zinit[15,16,0:2]=[0.036d,0.037d,0.038d]
        zinit[15,17,0:2]=[0.036d,0.037d,0.038d]
        zinit[15,18,0:2]=[0.036d,0.037d,0.038d]
        zinit[15,23:24,1:2]=0.039d
        zinit[16,0:4,1]=0.033d
        zinit[16,5,0]=0.033d
        zinit[16,5,1]=0.038d
        zinit[16,12,1]=0.035d
        zinit[16,15:16,1:2]=0.039d
        zinit[16,18,0:1]=0.036d
        zinit[16,19:22,2]=0.040d
        zinit[17,0:4,1]=0.033d
        zinit[17,11:21,0]=0.035d
        zinit[17,11:21,2]=0.039d
        zinit[18,0:9,1]=0.035d
        zinit[18,11:19,0]=0.035d
        zinit[18,11:19,2]=0.039d
        zinit[18,13,0:2]=0.037d
        zinit[19,0:9,1]=0.035d
        zinit[19,12:14,1]=0.035d
        zinit[19,20:26,1]=0.039d
        fluxsigthresh = 3d
        nad1thresh = 5d
        nad2thresh = 9999d
        nadxref = 9
        nadyref = 13
        nadfitran = [6050d,6150d]
     endif

     vdisp=75d
     outdir='/Users/drupke/winds/gmos/specfits/'+gal+'/'+outstr+'/'
     infile=rootindir+'cubes/'+gal+'/'+gal+outstr+'.fits'
     fcncontfit='gmos_fit_flat_continuum'
     plotstelfit='plotstelfit'
     fcnstrlines='gmos_plotstronglines'
     argstrlines=0
     argscontfit=0
     plotstelfit_args = {zbuf:1}
     fitran_rest=[5750,6700]
     siglim = [1.5d,15d]
     sigguess = 0
     ctoutfile=0
     qsoutfile=0
     qsoutfile_ha=0
     qsocntargs=0
     argsinitpar = {siglim:siglim}
     nadfitord = 3
     nadrestcomp = 0
     argslinelist = 0
     outlines=0

  endif


; COPY, 9 Feb 2012
;   if gal eq 'mrk273' then begin
;      if bin eq 3 then begin
; ;    Binning = 3x3
;         dx = 20
;         dy = 27
;         cx = 9d
;         cy = 13d
;         pscale = bin/10d           ; in arcseconds / pixel
;         outstr = 'rb'+string(bin,format='(I0)')
;         ncomp = dblarr(dx,dy)+2
; ; 1-2 comp
;         ncomp[0,3] = -2
;         ncomp[0,5:6] = -2
;         ncomp[0:3,22:26] = 1
;         ncomp[2,0] = -2
;         ncomp[6,13:14] = -2
;         ncomp[10,1] = -2
;         ncomp[13,16] = -2
;         ncomp[15,0] = -2
;         ncomp[19,0:9] = -2
; ; 3 comp

; ; discrete redshifted feature:
; ;  12,1-4
; ;  13,1-5
; ;  14,1-5
; ;  15,1-5
; ;  16,1
;         ncomp[11,3] = 3
;         ncomp[12,3:4] = 3
;         ncomp[13,3:5] = 3
;         ncomp[14,3:4] = 3

; ; discrete redshifted feature
; ;   7,8-9
; ;   8,7-10
; ;   9,7-10
; ;  10,9
;         ncomp[6,7:8] = 3
;         ncomp[7,6:9] = 3
;         ncomp[8,6:9] = 3
;         ncomp[9,8] = 3

; ;        ncomp[15,3:4] = 3
; ;        ncomp[16,3:4] = 3

; ; discrete blueshifted feature
;         ncomp[13,22:24] = 3
; ; old
;         ncomp[14,21:24] = 3
; ; newer try
;         ; ncomp[14,21] = 3
;         ; ncomp[14,22] = 4

; ; three dynamical components?
;         ncomp[9,9:10] = 3
;         ncomp[10,6:11] = 3
;         ncomp[11,7:11] = 3
; ;        ncomp[12,7:10] = 3
; ;        ncomp[13,8:11] = 3

; ; combine two red components
;         ncomp[4,10:11] = 3
;         ncomp[4,14:16] = 3

; ; three dynamical components?
;         ncomp[8,18:21] = 3
;         ncomp[9,17:22] = 3
;         ncomp[10,18:22] = 3
;         ncomp[11,18:22] = 3
;         ncomp[12,19:22] = 3
;         ncomp[13,18:21] = 3
; ; old
;         ncomp[14,18:20] = 3
; ; newer try
; ;        ncomp[14,15:20] = 3

; ; three dynamical components?
;         ncomp[15,14:24] = 3
;         ncomp[16,11:26] = 3
;         ncomp[17,11:21] = 3
;         ncomp[17,16] = -3
;         ncomp[18,11:19] = 3
;         ncomp[19,12:19] = 3

; ;?
;         ncomp[13,25] = 3

; ; velocities
;         zinit=dblarr(dx,dy,4)+0.037d
;         zinit[0,3,1]=0.039d
;         zinit[0,9,1]=0.039d
;         zinit[1,10:12,1]=0.039d
;         zinit[4,23:26,1]=0.040d
;         zinit[5,20:24,1]=0.040d
;         zinit[6,6:8,2]=0.039d
;         zinit[7,9,2]=0.038d
;         zinit[7,24:26,*]=0.0375d
;         zinit[8,6:9,2]=0.038d
;         zinit[8,15,1]=0.040d
;         zinit[8,18:19,1]=0.033d        
;         zinit[8,24:26,*]=0.0375d
;         zinit[9,6:9,2]=0.038d
;         zinit[9,17,0:2]=[0.036d,0.037d,0.040d]
;         zinit[10,1,*]=0.038d
;         zinit[10,4,1]=0.033d
;         zinit[10,9,2]=0.033d
;         zinit[10,18:22,2]=0.035d
;         zinit[10,21,0:2]=[0.036d,0.037d,0.038d]
;         zinit[11,0:2,0:1]=0.039d
;         zinit[11,16:23,1]=0.040d
;         zinit[11,17,*]=0.038d
;         zinit[11,19,0:2]=[0.034d,0.036d,0.039d]
;         zinit[11,26,0:1]=0.038d
;         zinit[12,0:1,0:1]=0.039d
;         zinit[12,4,0:2]=0.039d
;         zinit[12,7:11,1]=0.035d
;         zinit[12,16:23,1]=0.038d
;         zinit[13,5,0:2]=[0.037d,0.037d,0.038d]
;         zinit[13,8:11,0]=0.035d
;         zinit[13,8:11,1]=0.038d
;         zinit[13,16,0:1]=0.038d
;         zinit[13,18:22,2] = 0.039d
;         zinit[13,25,2] = 0.039d
; ; newer try
;         ; zinit[14,15,0:2] = [0.038d,0.038d,0.039d]
;         ; zinit[14,16,0:2] = [0.037d,0.038d,0.039d]
;         zinit[14,18:19,1]=0.039d
; ; newer try
;         ; zinit[14,21,0:2]=[0.037d,0.037d,0.038d]
;         ; zinit[14,22,0:3]=[0.037d,0.038d,0.038d,0.039d]
;         ; zinit[14,23,0:2]=[0.0385d,0.0385d,0.0385d]
;         zinit[15,1,1]=0.038d
;         zinit[15,1:6,1]=0.033d
;         zinit[15,11,0]=0.035d
;         zinit[15,16,0:2]=[0.036d,0.037d,0.038d]
;         zinit[15,17,0:2]=[0.036d,0.037d,0.038d]
;         zinit[15,18,0:2]=[0.036d,0.037d,0.038d]
;         zinit[15,23:24,1:2]=0.039d
;         zinit[16,0:4,1]=0.033d
;         zinit[16,5,0]=0.033d
;         zinit[16,5,1]=0.038d
;         zinit[16,12,1]=0.035d
;         zinit[16,15:16,1:2]=0.039d
;         zinit[16,18,0:1]=0.036d
;         zinit[16,19:22,2]=0.040d
;         zinit[17,0:4,1]=0.033d
;         zinit[17,11:21,0]=0.035d
;         zinit[17,11:21,2]=0.039d
;         zinit[18,0:9,1]=0.035d
;         zinit[18,11:19,0]=0.035d
;         zinit[18,11:19,2]=0.039d
;         zinit[18,13,0:2]=0.037d
;         zinit[19,0:9,1]=0.035d
;         zinit[19,12:14,1]=0.035d
;         zinit[19,20:26,1]=0.039d
;         fluxsigthresh = 3d
;         nad1thresh = 5d
;         nad2thresh = 9999d
;         nadxref = 9
;         nadyref = 13
;         nadfitran = [6050d,6150d]
;      endif

;      vdisp=75d
;      outdir='/Users/drupke/winds/gmos/specfits/'+gal+'/'+outstr+'/'
;      infile=rootindir+'cubes/'+gal+'/'+gal+outstr+'.fits'
;      fcncontfit='gmos_fit_flat_continuum'
;      plotstelfit='plotstelfit'
;      fcnstrlines='gmos_plotstronglines'
;      argstrlines=0
;      argscontfit=0
;      plotstelfit_args = {zbuf:1}
;      fitran_rest=[5750,6700]
;      siglim = [1.5d,15d]
;      sigguess = 0
;      ctoutfile=0
;      qsoutfile=0
;      qsoutfile_ha=0
;      qsocntargs=0
;      argsinitpar = {siglim:siglim}
;      nadfitord = 3
;      nadrestcomp = 0
;      argslinelist = 0
;      outlines=0

;   endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; F08572:NW
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if gal eq 'f08572nw' then begin
     if bin eq 3 then begin
        dx = 22
        dy = 16
        cx = 13d
        cy = 8d
     endif

     pscale = bin/10d           ; in arcseconds / pixel
     ncomp = dblarr(dx,dy)+1
     ncomp[6,15]=0
     ncomp[8,8]=-2
     ncomp[9,5:10]=-2
     ncomp[10,4:11]=-2
     ncomp[11,4:11]=-2
     ncomp[12,5:10]=-2
;     ncomp[13,5:8]=-2
     zinit=dblarr(dx,dy,2)
     zinit[*,*,0] = 0.058d
     zinit[*,*,1] = 0.052d
     vdisp=75d
     outstr = 'rb'+string(bin,format='(I0)')
     outdir='/Users/drupke/winds/gmos/specfits/'+gal+'/'+outstr+'/'
     infile=rootindir+'cubes/'+gal+'/'+gal+outstr+'.fits'
     fcncontfit='gmos_fit_flat_continuum'
     plotstelfit='plotstelfit'
     fcnstrlines='gmos_plotstronglines_cfoc'
     argstrlines={cfoc:2,woff:[-120d,20d]}
     argscontfit={fitord:3}
     plotstelfit_args = {zbuf:1}
     fitran_rest=[5400,6590]
     fluxsigthresh = 3d
     nad1thresh = 4d
     nad2thresh = 9999d
     nadxref = 10
     nadyref = 10
     nadfitran = [6180d,6280d]
     nadfitord = 1
     nadrestcomp = 1
     siglim = [1.5d,30d]
     ctoutfile=0
     qsoutfile=0
     qsoutfile_ha=0
     qsocntargs=0
     argsinitpar = {siglim:siglim,n2hafix:[0,0.9],sigfix:[0,1]}
     sigguess = dblarr(7,2)
     sigguess[*,0] = 2
     sigguess[*,1] = 15
     argslinelist = 0
     outlines=0

  endif



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; F08572:SE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if gal eq 'f08572se' OR gal eq 'f08572se' then begin
     if bin eq 3 then begin
        dx = 24
        dy = 16
        cx = 13d
        cy = 8d
     endif

     pscale = bin/10d           ; in arcseconds / pixel
     ncomp = dblarr(dx,dy)+1
     zinit=dblarr(dx,dy,2)+0.059d
     vdisp=75d
     outstr = 'rb'+string(bin,format='(I0)')
     outdir='/Users/drupke/winds/gmos/specfits/'+gal+'/'+outstr+'/'
     infile=rootindir+'cubes/'+gal+'/'+gal+outstr+'.fits'
     fcncontfit='gmos_fit_flat_continuum'
     plotstelfit='plotstelfit'
     fcnstrlines='gmos_plotstronglines'
     argscontfit={fitord:3}
     plotstelfit_args = {zbuf:1}
     argstrlines=0
     fitran_rest=[5400,6590]
     fluxsigthresh = 3d
     nad1thresh = 4d
     nad2thresh = 9999d
     nadxref = 10
     nadyref = 10
     nadfitran = [6180d,6280d]
     nadfitord = 1
     nadrestcomp = 1
     siglim=0
     ctoutfile=0
     qsoutfile=0
     qsoutfile_ha=0
     qsocntargs=0
     argsinitpar = {polypars:[2,6840d,6980d]}
     sigguess = 0
     argslinelist = 0
     outlines=0

  endif




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; VV705:NW
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if gal eq 'vv705nw' then begin
     if bin eq 3 then begin
        dx = 22
        dy = 16
        cx = 11d
        cy = 8d
     endif

     pscale = bin/10d           ; in arcseconds / pixel
     ncomp = dblarr(dx,dy)-2
     ncomp[0,5:15]=1
     ncomp[1,5:15]=1
     ncomp[2,10:15]=1
     ncomp[3:7,12:15]=1
     ncomp[6,10]=1
     ncomp[21,7]=1
     zinit=dblarr(dx,dy,2)
     zinit[*,*,0] = 0.040d
     zinit[*,*,1] = 0.040d
     zinit[3,5,*] = 0.041d
     zinit[5,0,*] = 0.041d
     zinit[18,6,*] = 0.041d
     vdisp=75d
     outstr = 'rb'+string(bin,format='(I0)')
     outdir='/Users/drupke/winds/gmos/specfits/'+gal+'/'+outstr+'/'
     infile=rootindir+'cubes/'+gal+'/'+gal+outstr+'.fits'
     fcncontfit='gmos_fit_flat_continuum'
     plotstelfit='plotstelfit'
     fcnstrlines='gmos_plotstronglines_cfoc'
     argstrlines = 0
     argscontfit={fitord:4}
     plotstelfit_args = {zbuf:1}
     fitran_rest=[5475,6625]
     fluxsigthresh = 3d
     nad1thresh = 4d
     nad2thresh = 9999d
     nadxref = 11
     nadyref = 8
     nadfitran = [6080d,6180d]
     nadfitord = 2
     nadrestcomp = 1
     siglim = [0.7d,10d]
     ctoutfile = 0
     qsoutfile = 0
     qsoutfile_ha = 0
     qsocntargs = 0
     argsinitpar = {siglim:siglim}
     sigguess = 0
     argslinelist = 0
     outlines=0

  endif




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; VV705:SE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if gal eq 'vv705se' then begin
     if bin eq 3 then begin
        dx = 22
        dy = 16
        cx = 8.5d
        cy = 5.5d
     endif

     pscale = bin/10d           ; in arcseconds / pixel
     ncomp = dblarr(dx,dy)+1
     ncomp[2,0:3]=-2
     ncomp[3:13,0:9]=-2
     ncomp[10:14,13:15]=-2
     ncomp[15:19,12:15]=-2
     ncomp[3:5,4]=1
     zinit=dblarr(dx,dy,2)
     zinit[*,*,0] = 0.0410d
     zinit[*,*,1] = 0.0405d
     zinit[7,4,*] = 0.0405d
     zinit[8,0:1,*] = 0.0405d
     zinit[10:14,13:15,*] = 0.0403d
     zinit[15:19,12:15,*] = 0.0403d
     vdisp=75d
     outstr = 'rb'+string(bin,format='(I0)')
     outdir='/Users/drupke/winds/gmos/specfits/'+gal+'/'+outstr+'/'
     infile=rootindir+'cubes/'+gal+'/'+gal+outstr+'.fits'
     fcncontfit='gmos_fit_flat_continuum'
     plotstelfit='plotstelfit'
     fcnstrlines='gmos_plotstronglines_cfoc'
     argstrlines = 0
     argscontfit={fitord:4}
     plotstelfit_args = {zbuf:1}
     fitran_rest=[5475,6625]
     fluxsigthresh = 3d
     nad1thresh = 4d
     nad2thresh = 9999d
     nadxref = 10
     nadyref = 5
     nadfitran = [6080d,6180d]
     nadfitord = 2
     nadrestcomp = 1
     siglim = [1.5d,15d]
     ctoutfile = 0
     qsoutfile = 0
     qsoutfile_ha = 0
     qsocntargs = 0
     argsinitpar = 0
     sigguess = 0
     argslinelist = 0
     outlines=0

  endif




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; F17207:E
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if gal eq 'f17207e' then begin
     if bin eq 3 then begin
        dx = 22
        dy = 16
        cx = 15d
        cy = 3d
        nadxref = 15
        nadyref = 3
     endif else if bin eq 6 then begin
        dx = 11
        dy = 8
        cx = 7.5d
        cy = 1.5d
        nadxref = 8
        nadyref = 2
     endif

     pscale = bin/10d           ; in arcseconds / pixel
     ncomp = dblarr(dx,dy)+1
     zinit=dblarr(dx,dy,2)+0.043
     vdisp=75d
     outstr = 'rb'+string(bin,format='(I0)')
     outdir='/Users/drupke/winds/gmos/specfits/'+gal+'/'+outstr+'/'
     infile=rootindir+'cubes/'+gal+'/'+gal+outstr+'.fits'
     fcncontfit='gmos_fit_flat_continuum'
     plotstelfit='plotstelfit'
     fcnstrlines='gmos_plotstronglines'
     argstrlines = 0
     argscontfit={fitord:4}
     plotstelfit_args = {zbuf:1}
     fitran_rest=[5400,6640]
     fluxsigthresh = 3d
     nad1thresh = 4d
     nad2thresh = 9999d
     nadfitran = [6090d,6190d]
     nadfitord = 2
     nadrestcomp = 1
     siglim = [1.5d,15d]
     ctoutfile = 0
     qsoutfile = 0
     qsoutfile_ha = 0
     qsocntargs = 0
     argsinitpar = 0
     sigguess = 0
     argslinelist = 0
     outlines = 0

  endif




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; F17207:W
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if gal eq 'f17207w' then begin
     if bin eq 3 then begin
        dx = 22
        dy = 16
        cx = 15d
        cy = 16d 
        nadxref = 15
        nadyref = 16
    endif else if bin eq 6 then begin
        dx = 11
        dy = 8
        cx = 7.5d
        cy = 8d
        nadxref = 8
        nadyref = 8
     endif

     pscale = bin/10d           ; in arcseconds / pixel
     ncomp = dblarr(dx,dy)+1
     zinit=dblarr(dx,dy,2)+0.043
     vdisp=75d
     outstr = 'rb'+string(bin,format='(I0)')
     outdir='/Users/drupke/winds/gmos/specfits/'+gal+'/'+outstr+'/'
     infile=rootindir+'cubes/'+gal+'/'+gal+outstr+'.fits'
     fcncontfit='gmos_fit_flat_continuum'
     plotstelfit='plotstelfit'
     fcnstrlines='gmos_plotstronglines'
     argstrlines = 0
     argscontfit={fitord:4}
     plotstelfit_args = {zbuf:1}
     fitran_rest=[5400,6640]
     fluxsigthresh = 3d
     nad1thresh = 4d
     nad2thresh = 9999d
     nadfitran = [6090d,6190d]
     nadfitord = 2
     nadrestcomp = 1
     siglim = [1.5d,15d]
     ctoutfile = 0
     qsoutfile = 0
     qsoutfile_ha = 0
     qsocntargs = 0
     argsinitpar = 0
     sigguess = 0
     argslinelist = 0
     outlines = 0

  endif





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; F10565
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if gal eq 'f10565' then begin
     if bin eq 3 then begin
        dx = 22
        dy = 26
        cx = 11d
        cy = 14d 
        nadxref = 12
        nadyref = 14
     endif

     pscale = bin/10d           ; in arcseconds / pixel
     ncomp = dblarr(dx,dy)+1
     ncomp[ 1,4:12]=-2
     ncomp[ 2,4:15]=-2
     ncomp[ 3,4:15]=-2
     ncomp[ 4,4:15]=-2
     ncomp[ 5,4:16]=-2
     ncomp[ 6,4:17]=-2
     ncomp[ 7,4:19]=-2
     ncomp[ 8,4:19]=-2
     ncomp[ 9,4:19]=-2
     ncomp[10,4:20]=-2
     ncomp[11,2:20]=-2
     ncomp[12,2:20]=-2
     ncomp[13,2:20]=-2
     ncomp[14,3:18]=-2
     ncomp[15,3:18]=-2
     ncomp[16,3:18]=-2
     ncomp[17,3:18]=-2
     ncomp[18,3:18]=-2
     ncomp[19,3:16]=-2
     ncomp[20,4:16]=-2
     ncomp[21,4:16]=-2
     zinit=dblarr(dx,dy,2)+0.043
     zinit[1,9,0]=0.042
     zinit[2,6,0]=0.041
     vdisp=75d
     outstr = 'rb'+string(bin,format='(I0)')
     outdir='/Users/drupke/winds/gmos/specfits/'+gal+'/'+outstr+'/'
     infile=rootindir+'cubes/'+gal+'/'+gal+outstr+'.fits'
     fcnstrlines='gmos_plotstronglines'
     argstrlines = 0
     ; fcncontfit='gmos_fit_continuum'
     ; argscontfit={fitord:4}
     ; fcncontfit='gmos_fit_flat_continuum_gaps'
     ; fitord=5
     ; expterms=2
     ; gaps=[[6015d,6055d],[6515d,6555d]]
     ; argscontfit = {fitord:fitord,$
     ;                expterms:expterms,$
     ;                gaps:gaps,$
     ;                z:{star:0d,gas:0d} $
     ;               }
     fcncontfit='gmos_fit_flat_continuum'
     argscontfit={fitord:4,refito1:1,z:{star:0d,gas:0d}}
     plotstelfit='plotstelfit'
     plotstelfit_args = {zbuf:1} 
     fitran_rest=[5400,6640]
     fluxsigthresh = 3d
     nad1thresh = 3d
     nad2thresh = 9999d
     nadfitran = [6090d,6190d]
     nadfitord = 2
     nadrestcomp = 1
;     ncompnad = dblarr(dx,dy)+1
     siglim = [1.5d,30d]
     ctoutfile = 0
     qsoutfile = 0
     qsoutfile_ha = 0
     qsocntargs = 0
     argsinitpar = 0
;     argsinitpar = {polypars:[3,6750d,6930d]}
     sigguess = 0
     argslinelist = 0
     outlines = 0

  endif

; Make sure data exists  
  if ~ file_test(infile) then begin
     print,"GMOS_INITFIT_SPECTRA: Data cube not found."
     stop
  endif

  init = {zinit            : zinit, $
          ncomp            : ncomp, $
          infile           : infile, $
          ctoutfile        : ctoutfile, $
          qsoutfile        : qsoutfile, $
          qsoutfile_ha     : qsoutfile_ha, $
          outdir           : outdir, $
          fcninitpar       : fcninitpar, $
          fcnlinefit       : fcnlinefit, $
          fcncontfit       : fcncontfit, $
          fcnstrlines      : fcnstrlines, $
          argsinitpar      : argsinitpar, $
          argscontfit      : argscontfit, $
          argstrlines      : argstrlines, $
          plotstelfit      : plotstelfit, $
          qsocntargs       : qsocntargs, $
          argslinelist     : argslinelist, $
          fitran_rest      : fitran_rest, $
          vdisp            : vdisp, $
          startempfile     : startempfile,$
          dx               : dx,$
          dy               : dy, $
          cx               : cx,$
          cy               : cy, $
          nad1thresh       : nad1thresh,$
          nad2thresh       : nad2thresh,$
          nadxref          : nadxref,$
          nadyref          : nadyref,$
          nadfitran        : nadfitran,$
          nadfitord        : nadfitord,$
          nadrestcomp      : nadrestcomp,$
          bin              : bin, $
          fluxsigthresh    : fluxsigthresh, $
          outlines         : outlines, $
          siglim           : siglim, $
          sigguess         : sigguess $
         }

  return,init

end
