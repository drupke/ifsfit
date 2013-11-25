; docformat = 'rst'
;
;+
;
; This function initializes the fitting parameters.
;
; :Categories:
;    UHSPECFIT/GMOS
;
; :Returns:
; A structure with the following required tags::
;
;   fcninitpar: in, required, type=string
;     Name of function for initializing continuum.
;   infile: in, required, type=string
;     Filename of input data cube.
;   linetie: in, required, type=strarr(dx,dy,nlines)
;     Name of emission line to which each emission line is tied
;     (in redshift and linewidth).
;   ncomp: in, required, type=dblarr(ncols,nrows,nlines)
;     For each spaxel and emission line, # of components to fit.
;   outdir: in, required, type=string
;     Directory for output save (.xdr) file
;   zinit_stars: in, required, type=double
;     Redshift used to shift any stellar templates to observed
;     frame.
;   zinit_gas: in, required, type=dblarr(ncols,nrows,nlines,ncomp)
;     Initial redshift guesses for each spaxel, emission line, and
;     component.
; 
; Also possibly included are the following optional tags::
;
;   argsaddpoly2temp: in, optional, type=structure
;     Arguments for UHSF_ADDPOLY2TEMP call.
;   argscontfit: in, optional, type=structure
;     Arguments for continuum fit routine.
;   argsinitpar: in, optional, type=structure
;     Arguments for parameter initialization routine.
;   argslinefit: in, optional, type=structure
;     Arguments for line fitting routine
;   argslinelist: in, optional, type=structure
;     Arguments for line selection routine
;   argsoptstelz: in, optional, type=structure
;     Arguments for stellar redshift optimization.
;   argspltlin1: in, optional, type=structure
;     Arguments for first line plot
;   argspltlin2: in, optional, type=structure
;     Arguments for first line plot
;   dividecont: in, optional, type=byte
;     Set this param to divide the data by the continuum
;     fit. Default is to subtract.
;   fcncontfit: in, optional, type=string
;     Name of continuum fitting function. If not specified,
;     continuum is not fit.
;   fcnlinefit: in, optional, type=string
;     Name of line fitting function. Default: UHSF_MANYGAUSS
;   fcnoptstelsig: in, optional, type=string
;     Name of routine to optimize stellar dispersion.
;   fcnoptstelz: in, optional, type=string
;     Name of routine to optimize stellar redshift. If not specified,
;     redshift is not optimized.
;   fcnpltcont: in, optional, type=string
;     Name of continuum plotting function. Default: UHSF_PLTCONT
;   fcnpltlin: in, optional, type=string
;     Name of line plotting function. Default: UHSF_PLTLIN
;   fitran: in, optional, type=dblarr(2)
;     Range of fitting, in observed frame. If not set, default is
;     entire range of data / template intersection.
;   keepnad: in, optional, type=string
;     Set to not remove NaD region from fit.
;   loglam: in, optional, type=byte
;     Set if data has constant log(lambda) dispersion.
;   maskwidths: in, optional, type=dblarr
;     Width, in km/s, of regions to mask from continuum fit. If not
;     set, routine defaults to +/- 500 km/s. If parameter has one
;     value, then this half-width is applied to all emission
;     lines. If it has multiple values, it should have exactly the
;     same number of elements as lines that are being fit.
;   nomaskran: in, optional, type=dblarr(2)
;     Wavelength region *not* to mask.
;   outlines: in, optional, type=strarr
;     Labels of emission lines for which to print line fluxes to output.
;   peakinit: in, optional, type=dblarr(nlines,ncomp)
;     Initial peak flux guesses.
;   siginit_gas: in, optional, type=dblarr(nlines,ncomp)
;     Initial line width guesses, in sigma and km/s.
;   siginit_stars: in, optional, type=double
;     Initial sigma value, in km/s, for a Gaussian kernel for
;     convolving with stellar template. Convolution only performed
;     if this param is set.
;   sigfitvals: in, optional, type=dblarr
;     If this param is set, routine cross-correlates data with
;     continua convolved with each sigma value in this array, and
;     chooses the sigma with the highest correlation coeff.
;   startempfile: in, optional, type=structure
;     File containing IDL save file (usually ending in .xdr) of
;     stellar templates. Tags:
;       lambda: in, required, type=dblarr(nwave)
;       flux: in, required, type=dblarr(nwave,ntemplates)
;   vacuum: in, optional, type=byte
;     Set this param to shift stellar templates from air to
;     vacuum wavelengths.
;
; :Params:
;    bin: in, required, type=integer
; 
; :Author:
;    David Rupke
;
; :History:
;    Change History::
;      2013sep, DSNR, complete re-write
;-
function uhsf_gm_init_f05189,bin

  gal = 'f05189'
  bin = double(bin)

  if bin eq 2 then begin
     dx = 28
     dy = 27
     cx = 14d
     cy = 14d 
  endif else begin 
     print,'UHSF_GM_INIT: Bin = ',string(bin,format='(I0)'),' not allowed.'
     return,0
  endelse

  outstr = 'rb'+string(bin,format='(I0)')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Required pars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Input file
  infile='/Users/drupke/winds/gmos/cubes/'+gal+'/'+gal+outstr+'.fits'
  if ~ file_test(infile) then begin
     print,"UHSF_GM_INIT: Data cube not found."
     return,0
  endif

; Initialize line ties, n_comps, and z_inits.
  nlines = 19
  maxncomp = 3
  linetie = strarr(dx,dy,nlines) + 'Halpha'
  ncomp = dblarr(dx,dy,nlines) + maxncomp
  zinit_gas=dblarr(dx,dy,nlines,maxncomp) + 0.0425d
  zinit_stars=dblarr(dx,dy) + 0.043d
; iron lines
  ncomp[*,*,0:3] = 1
  ncomp[13,13,0:3] = 2
  linetie[*,*,0:3] = '[FeVII]6087'
  zinit_gas[13,13,0:3,0] = 0.038d
  zinit_gas[13,13,0:3,1] = 0.040d
; HeII line
  ncomp[*,*,4] = 2
  linetie[*,*,4] = 'HeII4686'
  zinit_gas[*,*,4,0] = 0.040d
  zinit_gas[*,*,4,1] = 0.038d
; HeI lines
  ncomp[*,*,5:6] = 0
  ncomp[13,13,5:6] = 1
  linetie[*,*,5:6] = 'HeI6678'
; Balmer lines, low-ion. colliosional lines
  zinit_gas[*,*,7:nlines-1,1] = 0.040d
  zinit_gas[*,*,7:nlines-1,2] = 0.038d
; [NI] lines
  ncomp[*,*,10:11] = 1
  linetie[*,*,10:11] = '[NI]5198'


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Optional pars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Parameters for continuum fit
  ;; refit = {ran: [[4950,5100],[5250,5450],[5850,6000],$
  ;;                    [6200,6400],[6500,6700],[6725,6925],$
  ;;                    [6925,7100],[7250,7400]],$
  ;;              ord: [2,2,2,2,2,1,1,2]}
  refit = {ran: [[6500,6700]],$
           ord: [2]}

; Parameters for emission line plotting
  linoth = strarr(2,6)
  linoth[0,2] = '[OIII]4959'
  linoth[*,3] = ['[OI]6364','[FeX]6375']
  linoth[*,4] = ['[NII]6548','[NII]6583']
  linoth[*,5] = ['HeI6678','[SII]6716']
  argspltlin1 = {nx: 3, ny: 2,$
                 label: ['HeII4686','Hbeta','[OIII]5007',$
                         '[OI]6300','Halpha','[SII]6731'],$
                 wave: [4686,4861,5007,6300,6563,6731],$
                 off: [[-120,90],[-80,50],[-130,50],$
                       [-80,120],[-95,70],[-95,50]],$
                 linoth: linoth}
  linoth = strarr(3,6)
  linoth[*,0] = ['[NI]5198','[NI]5200','[FeVII]5159']
  argspltlin2 = {nx: 3, ny: 2,$
                 label: ['[FeVII]5159','[FeVII]5721','[FeVII]6087',$
                         'HeI7065','',''],$
                 wave: [5159,5721,6087,7065,0,0],$
                 off: [[-120,90],[-120,90],[-120,90],$
                       [-90,80],[-90,80],[-90,80]],$
                 linoth: linoth}

; Velocity dispersion guesses
  siginit_gas = dblarr(nlines,3)+150d
; iron lines
  siginit_gas[0:3,0] = 1000d
; HeII line
  siginit_gas[4,1] = 1000d
; HeI lines
  siginit_gas[5:6,0] = 1000d
; Balmer lines, low-ion. colliosional lines
  siginit_gas[7:nlines-1,2] = 1000d

; Velocity dispersion limits
  siglim_gas = [299792d/3000d/2.35d,2000d]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output structure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  init = {$
; Required pars
         fcninitpar: 'uhsf_gm_initpar',$
         fitran: [4600,7100],$
         infile: infile,$
         linetie: linetie,$
         ncomp: ncomp,$
         outdir: '/Users/drupke/winds/gmos/specfits/'+gal+'/'+outstr+'/',$
         zinit_stars: zinit_stars,$
         zinit_gas: zinit_gas,$
; Optional pars
         argscontfit: {refit: refit},$
         argsinitpar: {siglim: siglim_gas},$
         argslinelist: {felines: 1},$
         argsoptstelz: {lrange: [5200,5550]},$
         argspltlin1: argspltlin1,$
         argspltlin2: argspltlin2,$
         fcncontfit: 'uhsf_fitcont',$
         fcnoptstelsig: 'uhsf_optstelsig',$
         fcnoptstelz: 'uhsf_optstelz',$
         nomaskran: [5075,5100],$
         siglim_gas: siglim_gas,$
         siginit_gas: siginit_gas,$
         siginit_stars: 100d,$
;        first # is max sig, second is step size
         sigfitvals: dindgen(fix(500d/25d)+1)*25d,$
         startempfile: '/Users/drupke/src/idl/uhspecfit/stellar_models/'+$
         'gonzalezdelgado/SSPGeneva_z020.sav' $
         }

  return,init

end
