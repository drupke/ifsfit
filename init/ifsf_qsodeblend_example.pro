; docformat = 'rst'
;
;+
;
; This function initializes the fitting parameters for an example galaxy. 
; It uses the full nuclear spectrum (including all emission lines) as 
; the PSF template. 
; 
; :Categories:
;    IFSF
;
; :Returns:
;    A structure with tags specified in INITTAGS.txt.
;
; :Params:
; 
; :Keywords:
;    initmaps: out, optional, type=structure
;      Parameters for map making.
;    initnad: out, optional, type=structure
;      Parameters for NaD fitting.
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
;      2019nov22, DSNR, copied from ifsf_f13342.pro
;    
; :Copyright:
;    Copyright (C) 2019 David S. N. Rupke
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
function ifsf_qsodeblend_example,initmaps=initmaps,initnad=initnad

  bad=1d99

  gal = 'galaxyshorthand'
  bin = 3d
  ncols = 14
  nrows = 17
  centcol = 7.014
  centrow = 8.982
  platescale = 0.3d
  outstr = 'rb'+string(bin,format='(I0)')
  fitrange = [5400,8200]

; distance from central pixel
  x_pix = rebin(indgen(ncols)+1,ncols,nrows)
  y_pix = rebin(transpose(indgen(nrows)+1),ncols,nrows)
  rad_pix = sqrt((double(x_pix-centcol))^2d + (double(y_pix-centrow))^2d)

; Regions for setting components
  ireg0  = where(rad_pix lt 3d,ctreg0)
  ireg2 = where(rad_pix ge 7d,ctreg2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Required pars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Input file
  infile='/path/datacube.fits'
  if ~ file_test(infile) then message,'Data cube not found.'

; Lines to fit.
  lines = ['Halpha','Hbeta','HeII4686',$
           '[OI]6300','[OI]6364','[OIII]4959','[OIII]5007',$
           '[NI]5198','[NI]5200','[NII]6548','[NII]6583',$
           '[SII]6716','[SII]6731']
  nlines = n_elements(lines)

; Max no. of components.
  maxncomp = 2

; Initialize line ties, n_comps, z_inits, and sig_inits.
  linetie = hash(lines,'Halpha')
  ncomp = hash(lines)
  zinit_gas = hash(lines)
  siginit_gas = hash(lines)
; note that siginit_gas is technically optional, put here for convenience
  foreach i,lines do begin
     ncomp[i] = dblarr(ncols,nrows)+2
     zinit_gas[i] = dblarr(ncols,nrows,maxncomp) + 0.1797d
     siginit_gas[i] = dblarr(maxncomp) + 100d
     zinit_gas[i,*,*,1] = 0.1797d
     if ctreg2 gt 0 then for j=0,ctreg2-1 do $
       ncomp[i,x_pix[ireg2[j]]-1,y_pix[ireg2[j]]-1] = 1
  endforeach
  zinit_stars=dblarr(ncols,nrows) + 0.1797d

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Optional pars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Tweaked regions are around Hbeta/[OIII], HeI/NaD, [OI], Ha/[NII], and [SII]
;
;; Parameters for continuum fit
;  tweakcntfit = dblarr(ncols,nrows,3,5)
;; Default fitting order
;  tweakcntfit[*,*,2,*] = 2
;; Number of wavelength regions to re-fit
;  nregions = 5
;; Lower wavelength for re-fit
;  tweakcntfit[*,*,0,0:nregions-1] = $
;     rebin(reform([5650,6850,7400,7600,7850],1,1,1,nregions),$
;           ncols,nrows,1,nregions)
;; Upper wavelength for re-fit
;  tweakcntfit[*,*,1,0:nregions-1] = $
;     rebin(reform([6000,7050,7550,7850,8000],1,1,1,nregions),$
;           ncols,nrows,1,nregions)
;; Order for re-fit
;  tweakcntfit[*,*,2,0:nregions-1] = $
;     rebin(reform([1,1,1,1,1],1,1,1,nregions),$
;           ncols,nrows,1,nregions)


; Parameters for emission line plotting
  linoth = strarr(2,6)
  linoth[0,2] = '[OIII]4959'
  linoth[0,3] = '[OI]6364'
  linoth[0:1,4] = ['[NII]6548','[NII]6583']
  linoth[0,5] = '[SII]6716'
  argspltlin1 = {nx: 3, ny: 2,$
                 label: ['HeII4686','Hbeta','[OIII]5007',$
                         '[OI]6300','Halpha','[SII]6731'],$
                 wave: [4686,4861,5007,6300,6563,6731],$
                 off: [[-120,90],[-80,50],[-130,50],$
                       [-80,120],[-95,70],[-95,50]],$
                 linoth: linoth}
  linoth = strarr(1,6)
  linoth[0,0] = '[NI]5198'
  linoth[0,1] = '[NII]5755'
  argspltlin2 = {nx: 3, ny: 2,$
                 label: ['[NI]5200','','','','',''],$
                 wave: [5159,0,0,0,0,0],$
                 off: [[-120,90],[-120,90],[-120,90],$
                       [-90,80],[-90,80],[-90,80]],$
                 linoth: linoth}

; Velocity dispersion limits and fixed values
  siglim_gas = dblarr(ncols,nrows,2)
  siglim_gas[*,*,0] = 5d
  siglim_gas[*,*,1] = 1000d
;  lratfix=hash()
;  lratfix['[NI]5200/5198'] = [1.5d]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output structure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  init = {$
;        Required pars
         fcninitpar: 'ifsf_gmos',$
         fitran: fitrange,$
         fluxunits: 1d-15, $ ; erg/s/cm^2/arcsec^2
         infile: infile,$
         label: gal,$
         lines: lines,$
         linetie: linetie,$
         maxncomp: maxncomp,$
         name: 'galaxyname',$
         ncomp: ncomp,$
         mapdir: '/path/mapdir/'+gal+'/'+outstr+'/',$
         outdir: '/path/fitdir/'+gal+'/'+outstr+'/',$
         platescale: platescale,$
         positionangle: 0d,$
         zinit_stars: zinit_stars,$
         zinit_gas: zinit_gas,$
         zsys_gas: 0.1797d,$
;        Optional pars
         argscheckcomp: {sigcut: 3,$
                         ignore: ['[OI]6300','[OI]6364',$
                                  '[SII]6716','[SII]6731']},$ ;struct
         argscontfit: {qsoxdr: '/path/nucleartemplate.xdr',$
                       siginit_stars: 50d,$
                       uselog: 1b $
                       },$
; Comment out for step 2; uncomment for step 7
;                       refit: 1b},$ ;comment out for step (2)
; Comment out for step 2
;         startempfile: '/path/stellar_models/'+$
;                       'gonzalezdelgado/SSPGeneva_z020.sav',$
; Comment out for step 2; uncomment for step 7
;         startempfile: '/path/starlighttemplate.xdr', $
         argspltlin1: argspltlin1,$
         argspltlin2: argspltlin2,$
         donad: 1,$
         decompose_qso_fit: 1b,$
         maskwidths_def: 500d,$
         fcncheckcomp: 'ifsf_checkcomp',$
         fcncontfit: 'ifsf_fitqsohost',$
 ;        tweakcntfit: tweakcntfit,$
         emlsigcut: 2d,$
         logfile: '/path/fitdir/'+gal+'/'+outstr+'/'+$
                  gal+'_fitlog.txt',$
         batchfile: '/path/to/ifsfit/common/ifsf_fitloop.pro',$
         batchdir: '/path/to/batch/directory/',$
;         cutrange: [[6920,6965]], $
         nocvdf: 1b,$
         siglim_gas: siglim_gas,$
         siginit_gas: siginit_gas,$
         siginit_stars: 100d, $
; Comment out for step 7
         host: {dat_fits: '/path/cubes/'+gal+'/'+gal+'starlight.fits'} $
; Uncomment for step 7
;         host: {dat_fits: '/path/cubes/'+gal+'/'+gal+'starlight_iter1.fits'} $
        }

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Arguments for maps
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
   if keyword_set(initmaps) then begin

      sorttype = hash()
      sorttype['Halpha'] = 'sigma'
      
      contourlevels = hash()
      contourlevels['Halpha_vpk'] = $
         [-200,-150,-100,-50,0,50,100,150,200]
      contourlevels['[NII]6583_v%50c1'] = contourlevels['Halpha_vpk']
      contourlevels['[NII]6583_v%50c2'] = [0]
      contourlevels['[OIII]5007_v%50c1'] = contourlevels['Halpha_vpk']
      contourlevels['[OIII]5007_v%50c2'] = [0]

      initmaps = {$
                  ebv: {calc: ['ftot','fc1','fc2'],$
                        titles: ['F$\downtot$','F$\downc1$','F$\downc2$'],$
                        apply: 1},$
                  lr: {calc: ['ftot','fc1','fc2']},$
                  cvdf: {flux_maps: [750,-750],$
                         sigcut: 1},$
                  emlplot: {ftags: ['ftot','fc1','fc2','fv-750'],$
                            vtags: ['vsigc1','vsigc2','v%50c1','v%50c2','v%98c2'],$
                            radtags: ['ftot'],$
                            ftitles: ['F$\uptot$',$
                                      'F$\upc1$',$
                                      'F$\upc2$',$
                                      'F$\down\\lambda$$\up-750km/s$'],$
                            vtitles: ['$\sigma$$\upc1$',$
                                      '$\sigma$$\upc2$',$
                                      'v$\down50$$\upc1$',$
                                      'v$\down50$$\upc2$',$
                                      'v$\down98$$\upc2$'],$
                            radtitles: ['F$\uptot$']},$
                  aspectrat: double(nrows)/double(ncols),$
                  center_axes: [centcol,centrow],$
                  center_nuclei: [centcol,centrow],$
                  contourlevels: contourlevels,$
                  ctradprof_psffwhm: 0.6d,$
                  nadabsweq_snrthresh: 3d,$
                  nademweq_snrthresh: 3d,$
                  ct: {sumrange: fitrange,$
                       sumrange_hstcomp: [7000,8230],$
                       scllim: [0,1],$
                       scllim_rad: [-2.5,0],$
                       stretch: 1,$
                       fitifspeak: 1b,$
                       fitifspeakwin_kpc: 5d},$
                  hst: {refcoords: [985.305,1443.95],$
                        subim_sm: 7d,$
                        subim_big: 20d,$
                        smoothfwhm: 6,$
                        fithstpeak: 1b,$
                        fithstpeakwin_kpc: 2d},$
                  hstrd: {file: '/path/hst.fits',$
                          label: 'WFPC2/F814W',$
                          scllim: [0.01,50],$
                          sclargs_sm: {beta: 0.05,stretch: 5},$
                          sclargs_big: {beta: 0.05,stretch: 5},$
                          sclargs_fov: {beta: 0.05,stretch: 5},$
                          photflam: 3.5224d-19,$
                          photplam: 7995.943,$
                          platescale: 0.10d,$
                          nucoffset: [0d,0d]},$
                  hstrdsm: {scllim: [0,20],$
                            sclargs: {beta: 0.5,stretch: 5}},$
                  fcnsortcomp: 'ifsf_sortcomp',$
                  rangefile: '/path/ranges.txt',$
                  fluxunits: 1d-15, $
                  outflow_eml: {R: 10.7d,$
                                use_spaxel_elecden: 1}, $ ;, usered: 1b},$
                  outflow_abs: {R: 10.7d} $
                 }
   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Parameters for NaD + HeI 5876 fit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if keyword_set(initnad) then begin

      normnadlo = [6810,6910]
      normnadhi = [6960,7060]
      pltnormnad = [6810,7060]
      maxncomp = 2

;     Initialize n_comps, z_inits, and sig_inits.
;     Use 1 HeI component
      heitie = strarr(ncols,nrows) + 'Halpha'
      hei_zinit = dblarr(ncols,nrows,maxncomp)
      hei_siginit = dblarr(ncols,nrows,maxncomp)

      nnadabs = intarr(ncols,nrows) + 1
      nadabs_zinit = dblarr(ncols,nrows,maxncomp)+0.1797d
      nadabs_siginit = dblarr(ncols,nrows,maxncomp)+50d
      nadabs_siglim = [5d,1000d]
      nadabs_cfinit = dblarr(ncols,nrows,maxncomp)+0.5d
      nadabs_tauinit = dblarr(ncols,nrows,maxncomp)+0.5d

      initnad = {$
                 argsinitpar: {siglimhei: [5d,1000d]},$
                 argsnadweq: {autowavelim: [6930,6960,6960,7000],$
                              autoindices:1,$
                              snrabsthresh: 1d},$
                 argsnormnad: {fitranlo: normnadlo,$
                               fitranhi: normnadhi,$
                               snavg_thresh: 1d},$
                 argspltnormnad: {fitranlo: normnadlo,$
                                  fitranhi: normnadhi,$
                                  pltran: pltnormnad,$
                                  fitord: 3},$
                 noautoreject: 1b,$
                 nadfitran: [6910,6960],$
                 argspltfitnad: {yran: [0,2]},$
                 fcnfitnad: 'ifsf_nadfcn',$
                 fcninitpar: 'ifsf_initnad',$
                 maxncomp: maxncomp,$
                 mcniter: 200, $
                 zref: 0.1797d, $
;                NaD absorption
                 nnadabs: nnadabs,$
                 nadabs_cfinit: nadabs_cfinit,$
                 nadabs_tauinit: nadabs_tauinit,$
                 nadabs_zinit: nadabs_zinit,$
                 nadabs_siginit: nadabs_siginit,$
                 nadabs_siglim: nadabs_siglim,$
;                NaD emission
                 nnadem: intarr(ncols,nrows),$
;                HeI
                 hei_zinit: hei_zinit,$
                 hei_siginit: hei_siginit,$
                 heitie: heitie $
                }
   endif
                  
   return,init

end
