; docformat = 'rst'
;
;+
;
; This procedure is the core routine to fit the continuum and emission
; lines of a spectrum. The 
;
;      fcninitpar: in, required, type=string
;        Name of function for initializing continuum.
;
;      fitran_rest: in, required, type=dblarr(2)
;        Range of fitting, in rest frame
;      infile: in, required, type=string
;        Filename of input data cube.
;      linetie: in, required, type=strarr(nlines)
;        Name of emission line to which each emission line is tied
;        (in redshift and linewidth).
;      ncomp: in, required, type=dblarr(ncols,nrows,nlines)
;        For each spaxel and emission line, # of components to fit.
;      outdir: in, required, type=string
;        Directory for output save (.xdr) file
;      zinit_stars: in, required, type=double
;        Redshift used to shift any stellar templates to observed
;        frame.
;      zinit_gas: in, required, type=dblarr(ncols,nrows,nlines,ncomp)
;        Initial redshift guesses for each spaxel, emission line, and
;        component.
;      startempfile: in, optional, type=structure
;        File containing IDL save file (usually ending in .xdr) of
;        stellar templates. Tags:
;          lambda: in, required, type=dblarr(nwave)
;          flux: in, required, type=dblarr(nwave,ntemplates)
;      argsaddpoly2temp: in, optional, type=structure
;        Arguments for UHSF_ADDPOLY2TEMP call.
;      argscontfit: in, optional, type=structure
;        Arguments for continuum fit routine.
;      argsinitpar: in, optional, type=structure
;        Arguments for parameter initialization routine.
;      argslinefit: in, optional, type=structure
;        Arguments for line fitting routine
;      argslinelist: in, optional, type=structure
;        Arguments for line selection routine
;      argsoptstelz: in, optional, type=structure
;        Arguments for stellar redshift optimization.
;      fcncontfit: in, optional, type=string
;        Name of continuum fitting function. If not specified,
;        continuum is not fit.
;      fcnlinefit: in, optional, type=string
;        Name of line fitting function. Default: UHSF_MANYGAUSS
;      fcnoptstelz: in, optional, type=string
;        Name of routine to optimize stellar redshift.
;      fitran: in, optional, type=dblarr(2)
;        Range of fitting, in observed frame. If not set, default is
;        entire range of data / template intersection.
;      keepnad: in, optional, type=string
;        Set to not remove NaD region from fit.
;      loglam: in, optional, type=byte
;        Set if data has constant log(lambda) dispersion.
;      masklines: in, optional, type=strarr
;        Lines to mask. If not set, all lines in the LINELIST
;        parameter are selected.
;      maskwidths: in, optional, type=dblarr
;        Width, in km/s, of regions to mask from continuum fit. If not
;        set, routine defaults to +/- 500 km/s. If parameter has one
;        value, then this half-width is applied to all emission
;        lines. If it has multiple values, it should have exactly the
;        same number of elements as lines that are being fit.
;      nomaskran: in, optional, type=dblarr(2)
;        Wavelength region *not* to mask.
;      peakinit: in, optional, type=dblarr(nlines,ncomp)
;        Initial peak flux guesses.
;      siginit_gas: in, optional, type=dblarr(nlines,ncomp)
;        Initial line width guesses, in sigma and km/s.
;      siginit_stars: in, optional, type=double
;        Initial sigma value, in km/s, for a Gaussian kernel for
;        convolving with stellar template. Convolution only performed
;        if this param is set.
;      sigfitvals: in, optional, type=dblarr
;        If this param is set, routine cross-correlates data with
;        continua convolved with each sigma value in this array, and
;        chooses the sigma with the highest correlation coeff.
;      dividecont: in, optional, type=byte
;        Set this param to divide the data by the continuum
;        fit. Default is to subtract.
;      vacuum: in, optional, type=byte
;        Set this param to shift stellar templates from air to
;        vacuum wavelengths.
;
; :Categories:
;    UHSPECFIT
;
; :Returns:
;    A structure that contains the fit and much else ...
;
; :Params:
; 
; :Author:
;    David Rupke
;
; :History:
;    Change History::
;      2013sep, DSNR, complete re-write

function uhsf_gm_init_f05189,gal,bin,sky=sky

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

;
  fcninitpar='uhsf_gm_initpar'
  fcncontfit='uhsf_fitcont' ; fit polynomial
  plotstelfit='uhsf_pltcont'
  fcnstrlines='uhsf_pltskylin'
  plotstelfit_args = {zbuf:1}
  argscontfit = {fitord:2}
  argstrlines = 0
  argslinelist = 0
  argstelcc = 0
;
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
  siglim_gas = [0.699d,0.701d]
  sigguess = dblarr(7)+0.7d
  argsinitpar = {siglim_gas:siglim_gas}
  outlines = ['[OI]6300']
  dodisptemp = 0
  nomaskran = 0

  galname = strmid(gal,0,6)
  exptag = strmid(gal,6,1)
  exptagmatch = {a: '1a', b: '1b', c: '1c', d: '1d', e: '3a', f: '3b', g: '4a'}
  tagnames = tag_names(exptagmatch)
  if galname eq 'f05189' AND exptag ne '' then begin
     zinit = dblarr(743)
     ncomp = intarr(743)+1
     itag = where(strcmp(exptag,tagnames,/fold_case),cttag)
     infile=rootindir+'red/'+galname+'/ctexrdat_'+exptagmatch.(itag)+'.fits'
     outdir = '/Users/drupke/winds/gmos/specfits/'+galname+'/o1_5577/exp'+exptag+'/'
     outlines = ['[OI]5577']
     fitran_rest=[5550,5600]
     siglim_gas = [0.65d,0.75]
     sigguess = dblarr(10)+0.7d
     argslinelist = {addo1_5577:1}
     argsinitpar = {siglim_gas:siglim_gas}
     argstrlines = {outlines:outlines,argslinelist:argslinelist}
  endif



     outstr = 'rb'+string(bin,format='(I0)')
     outdir='/Users/drupke/winds/gmos/specfits/'+gal+'/'+outstr+'/'
     infile=rootindir+'cubes/'+gal+'/'+gal+outstr+'.fits'

     if bin eq 2 then begin
        dx = 28
        dy = 27
        cx = 14d
        cy = 14d 
        nadxref = 14
        nadyref = 14
     endif
     
     pscale = bin/10d           ; in arcseconds / pixel
 
     nlines = 19
     maxncomp = 3
     ncomp = dblarr(dx,dy,nlines) + maxncomp
     linetie = strarr(nlines)+'Halpha'
     zinit_gas=dblarr(dx,dy,nlines,maxncomp) + 0.0425d
     siginit_gas = dblarr(nlines,3)+150d
     siglim_gas = [299792d/3000d/2.35d,2000d]
;; ;    col = 14, row = 8
;;      zinitstar=dblarr(dx,dy) + 0.043d
;; ;    Balmer lines, low-ion. colliosional lines
;;      zinit_gas[*,*,7:nlines-1,1] = 0.040d
;;      zinit_gas[*,*,7:nlines-1,2] = 0.038d
;;      sigguess[7:nlines-1,2] = 20d
;; ;    [NI] lines
;;      ncomp[*,*,10:11] = 1
;;      linetie[10:11] = '[NI]5198'
;; ;    HeI lines
;;      ncomp[*,*,5:6] = 0
;;      linetie[5:6] = 'HeI6678'
;;      zinit_gas[*,*,5:6,0] = 0.040d
;;      sigguess[5:6,0] = 20d
;; ;    HeII line
;;      ncomp[*,*,4] = 0
;;      linetie[4] = 'HeII4686'
;; ;    iron lines
;;      ncomp[*,*,0:3] = 0
;;      linetie[0:3] = '[FeVII]6087'

;; ;    col = 14, row = 20
;;      vdisp=50d
     ;; zinitstar=dblarr(dx,dy) + 0.0427d
;; ;    Balmer lines, low-ion. colliosional lines
;;      zinit_gas[*,*,7:nlines-1,1] = 0.040d
;;      zinit_gas[*,*,7:nlines-1,2] = 0.038d
;;      sigguess[7:nlines-1,2] = 20d
;; ;    [NI] lines
;;      ncomp[*,*,10:11] = 1
;;      linetie[10:11] = '[NI]5198'
;; ;    HeI lines
;;      ncomp[*,*,5:6] = 0
;;      linetie[5:6] = 'HeI6678'
;;      zinit_gas[*,*,5:6,0] = 0.040d
;;      sigguess[5:6,0] = 20d
;; ;    HeII line
;;      ncomp[*,*,4] = 0
;;      linetie[4] = 'HeII4686'
;; ;    iron lines
;;      ncomp[*,*,0:3] = 0
;;      linetie[0:3] = '[FeVII]6087'

;    col = 14, row = 18
     zinit_stars=dblarr(dx,dy) + 0.043d
;    Balmer lines, low-ion. colliosional lines
     zinit_gas[*,*,7:nlines-1,1] = 0.040d
     zinit_gas[*,*,7:nlines-1,2] = 0.038d
     siginit_gas[7:nlines-1,2] = 1000d
;    [NI] lines
     ncomp[*,*,10:11] = 1
     linetie[10:11] = '[NI]5198'
;; ;    [OIII] lines
;;      ncomp[*,*,7:8] = 2
;;      linetie[7:8] = '[OIII]5007'
;;      zinit_gas[*,*,7:8,0] = 0.040d
;;      zinit_gas[*,*,7:8,1] = 0.038d
;;      sigguess[7:8,1] = 20d
;    HeI lines
     ncomp[*,*,5:6] = 0
     linetie[5:6] = 'HeI6678'
     zinit_gas[*,*,5:6,0] = 0.040d
     siginit_gas[5:6,0] = 1000d
;    HeII line
     ncomp[*,*,4] = 2
     linetie[4] = 'HeII4686'
     zinit_gas[*,*,4,0] = 0.040d
     zinit_gas[*,*,4,1] = 0.038d
     siginit_gas[4,1] = 1000d
;    iron lines
     ncomp[*,*,0:3] = 1
     linetie[0:3] = '[FeVII]6087'
     zinit_gas[*,*,0:3,0] = 0.040d
     siginit_gas[0:3,0] = 1000d

;    col = 14, row = 14
     ;; zinitstar=dblarr(dx,dy) + 0.0427d
     ;; vdisp=100d
;; ;    Balmer lines, low-ion. colliosional lines
;;      zinit_gas[*,*,7:nlines-1,1] = 0.040d
;;      zinit_gas[*,*,7:nlines-1,2] = 0.038d
;;      sigguess[7:nlines-1,2] = 20d
;; ;    HeI lines
;;      ncomp[*,*,5:6] = 1
;;      linetie[5:6] = 'HeI6678'
;;      zinit_gas[*,*,5:6,0] = 0.040d
;;      sigguess[5:6,0] = 20d
;; ;    HeII line
;;      ncomp[*,*,4] = 2
;;      linetie[4] = 'HeII4686'
;;      zinit_gas[*,*,4,0] = 0.040d
;;      zinit_gas[*,*,4,1] = 0.038d
;;      sigguess[4,1] = 20d
;; ;    iron lines
;;      ncomp[*,*,0:3] = 2
;;      linetie[0:3] = '[FeVII]6087'
;;      zinit_gas[*,*,0:3,0] = 0.040d
;;      zinit_gas[*,*,0:3,1] = 0.038d
;;      sigguess[0:3,0] = 10d
;;      sigguess[0:3,1] = 10d

     fcncontfit='uhsf_fitcont'
     plotstelfit='plotstelfit'
     pltlin='uhsf_pltlin'
     fcnstelcc='gmos_stelcc'
     linoth = strarr(2,6)
     linoth[0,2] = '[OIII]4959'
     linoth[*,3] = ['[OI]6364','[FeX]6375']
     linoth[*,4] = ['[NII]6548','[NII]6583']
     linoth[*,5] = ['HeI6678','[SII]6716']
     pltlin_par1 = {nx: 3, ny: 2,$
                    label: ['HeII4686','Hbeta','[OIII]5007',$
                            '[OI]6300','Halpha','[SII]6731'],$
                    wave: [4686,4861,5007,6300,6563,6731],$
                    off: [[-120,90],[-80,50],[-130,50],$
                          [-80,120],[-95,70],[-95,50]],$
                    linoth: linoth}
     linoth = strarr(3,6)
     linoth[*,0] = ['[NI]5198','[NI]5200','[FeVII]5159']
     pltlin_par2 = {nx: 3, ny: 2,$
                    label: ['[FeVII]5159','[FeVII]5721','[FeVII]6087',$
                            'HeI7065','',''],$
                    wave: [5159,5721,6087,7065,0,0],$
                    off: [[-120,90],[-120,90],[-120,90],$
                          [-90,80],[-90,80],[-90,80]],$
                    linoth: linoth}
     plotstelfit_args = {zbuf:1} 
     ;; fitran_rest=[4400,7120]
     fitran_rest=[4400,6800]
     fluxsigthresh = 3d
     nad1thresh = 3d
     nad2thresh = 9999d
     nadfitran = [6000d,6250d]
     nadfitord = 2
     nadrestcomp = 1
     ;; startempfile=''
     dodisptemp = 1
     nomaskran = [5075,5100]
     ;; tmp_refit = {ran: [[4950,5100],[5250,5450],[5850,6000],$
     ;;                    [6200,6400],[6500,6700],[6725,6925],$
     ;;                    [6925,7100],[7250,7400]],$
     ;;              ord: [2,2,2,2,2,1,1,2]}
     tmp_refit = {ran: [[6500,6700]],$
                  ord: [2]}
     ;; argscontfit={fitord:2,refit:tmp_refit}
     argscontfit={no_dust:1,refit:tmp_refit}
     argsinitpar = {siglim: siglim_gas}
     argslinelist = {twoslit: 1, felines: 1}
     argstelcc = {lrange: [5200,5550]}
     outlines = 0
;    
     siginit_stars = 100
     sigstep = 25
     maxsig = 500     
     sigfitvals = dindgen(fix(maxsig/sigstep)+1)*sigstep
     maskwidths = 500

  endif

; Make sure data exists  
  if ~ file_test(infile) then begin
     print,"GMOS_INITFIT_SPECTRA: Data cube not found."
     stop
  endif

  init = {zinit_stars      : zinit_stars, $
          zinit_gas        : zinit_gas, $
          siglim_gas       : siglim_gas, $
          siginit_gas      : siginit_gas, $
          siginit_stars    : siginit_stars, $
          sigfitvals       : sigfitvals, $
          maskwidths       : maskwidths, $
;
          ncomp            : ncomp, $
          infile           : infile, $
          ctoutfile        : ctoutfile, $
          qsoutfile        : qsoutfile, $
          qsoutfile_ha     : qsoutfile_ha, $
          outdir           : outdir, $
          fcninitpar       : fcninitpar, $
          fcncontfit       : fcncontfit, $
          fcnstelcc        : fcnstelcc, $
          plotstelfit      : plotstelfit, $
          argsinitpar      : argsinitpar, $
          argscontfit      : argscontfit, $
          argstrlines      : argstrlines, $
          qsocntargs       : qsocntargs, $
          argslinelist     : argslinelist, $
          argstelcc        : argstelcc, $
          pltlin           : pltlin, $
          pltlin_par1      : pltlin_par1, $
          pltlin_par2      : pltlin_par2, $

          fitran_rest      : fitran_rest, $
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
          linetie          : linetie, $
          dodisptemp       : dodisptemp, $
          nomaskran        : nomaskran $
         }

  return,init

end
