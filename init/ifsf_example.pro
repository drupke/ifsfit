; docformat = 'rst'
;
;+
;
; This function initializes the fitting parameters for an example galaxy
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
;      2017apr04, DSNR, created
;    
; :Copyright:
;    Copyright (C) 2017 David S. N. Rupke
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
function ifsf_example,initmaps=initmaps,initnad=initnad

  bad=1d99

  gal = 'example'
  bin = 3d
  ncols = 15
  nrows = 23
  centcol = 8d
  centrow = 12d
  platescale = 0.3d
  outstr = 'rb'+string(bin,format='(I0)')
  fitrange = [3750,6400]

; distance from central pixel
  x_pix = rebin(indgen(ncols)+1,ncols,nrows)
  y_pix = rebin(transpose(indgen(nrows)+1),ncols,nrows)
  rad_pix = sqrt((double(x_pix-centcol))^2d + (double(y_pix-centrow))^2d)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Required pars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Input file
  infile='/Location/of/data/cube/example.fits'
  if ~ file_test(infile) then message,'Data cube not found.'

; Lines to fit.
  lines = ['[OII]3726','[OII]3729','[NI]5198','[NI]5200',$
           'Hdelta','Hbeta','[OIII]4959','[OIII]5007']
  nlines = n_elements(lines)

; Max no. of components.
  maxncomp = 2

; Initialize line ties, n_comps, z_inits, and sig_inits.
  linetie = hash(lines,'[OIII]5007')
  ncomp = hash(lines)
  zinit_gas = hash(lines)
  siginit_gas = hash(lines)
; note that siginit_gas is technically optional, put here for convenience
  foreach i,lines do begin
     ncomp[i] = dblarr(ncols,nrows)+1
     zinit_gas[i] = dblarr(ncols,nrows,maxncomp) + 0.1234d
     siginit_gas[i] = dblarr(ncols,nrows,maxncomp) + 50d
     ncomp[i,5:10,10:15] = 2
     zinit_gas[i,5:10,10:15,1] = 0.1233d
     siginit_gas[i,5:10,10:15,1] = 200d
     ncomp[i,7:10,19:22] = 2
     zinit_gas[i,7:10,19:22,1] = 0.1232d
     siginit_gas[i,7:10,19:22,1] = 50d
     ncomp[i,7,13]=1
  endforeach
  zinit_stars=dblarr(ncols,nrows) + 0.1234d

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Optional pars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Tweaked regions are around ...

; Number of wavelength regions to re-fit
  nregions = 3
  ; Parameters for continuum fit
  tweakcntfit = dblarr(ncols,nrows,3,nregions)
  ; Default fitting order
  tweakcntfit[*,*,2,*] = 2
; Lower wavelength for re-fit
  tweakcntfit[*,*,0,0:nregions-1] = $
     rebin(reform([3788,4975,6065],1,1,1,nregions),$
           ncols,nrows,1,nregions)
; Upper wavelength for re-fit
  tweakcntfit[*,*,1,0:nregions-1] = $
     rebin(reform([3950,5275,6250],1,1,1,nregions),$
           ncols,nrows,1,nregions)
; Order for re-fit
  tweakcntfit[*,*,2,0:nregions-1] = $
     rebin(reform([2,2,2],1,1,1,nregions),$
           ncols,nrows,1,nregions)
;
; Parameters for emission line plotting
  linoth = strarr(1,6)
  linoth[0,0] = '[OII]3726'
  linoth[0,4] = '[OIII]4959'
  linoth[0,5] = '[NI]5198'
  argspltlin1 = {nx: 3, ny: 2,$
                 label: ['[OII]3729','[NeIII]3869','Hdelta',$
                         'Hbeta','[OIII]5007','[NI]5200'],$
                 wave: [3727,3869,4101,4861,4985,5199],$
                 off: [[-100,100],[-100,100],[-100,100],$
                       [-100,100],[-100,100],[-100,100]],$
                 linoth: linoth}

; Velocity dispersion limits and fixed values
  siglim_gas = [5d,1500d]
  lratfix = hash()
; 1 corresponds to n ~ 400 cm^-3; Pradhan et al. 2006, MNRAS, 366, L6
  lratfix['[OII]3729/3726']=[1.2d,1.2d]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output structure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  init = {$
; Required pars
         fcninitpar: 'ifsf_gmos',$
         fitran: fitrange,$
         fluxunits: 1d-15, $ ; erg/s/cm^2/arcsec^2
         infile: infile,$
         label: gal,$
         lines: lines,$
         linetie: linetie,$
         maxncomp: maxncomp,$
         name: 'Example',$
         ncomp: ncomp,$
         mapdir: '/Location/of/maps/',$
         outdir: '/Location/of/fit/outputs/',$
         platescale: platescale,$
         positionangle: 123d,$
         zinit_stars: zinit_stars,$
         zinit_gas: zinit_gas,$
         zsys_gas: 0.1234d,$
; Optional pars
         argscheckcomp: {sigcut: 3d},$
         argsinitpar: {lratfix: lratfix},$
         argspltlin1: argspltlin1,$
         fcncheckcomp: 'ifsf_checkcomp',$
         maskwidths_def: 500d,$
         tweakcntfit: tweakcntfit,$
         emlsigcut: 2d,$
         logfile: '/Location/of/logfile',$
         batchfile: '/Location/of/ifsfit/common/ifsf_fitloop.pro',$
         batchdir: '/Location/for/temp/batch/files/',$
         cvdf_vlimits: [-2d3,2d3],$
         cvdf_vstep: 1d,$
         siglim_gas: siglim_gas,$
         siginit_gas: siginit_gas,$
         siginit_stars: 50d $
       }

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Arguments for maps
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
   if keyword_set(initmaps) then begin

      initmaps = {$
                  cvdf: {flux_maps: [500,-500],$
                         sigcut: 1},$
                  emlplot: {ftags: ['ftot','fc1','fc2','fv-500'],$
                            vtags: ['vsig','v%50c1','v%50c2','v%98c2'],$
                            radtags: ['ftot'],$
                            ftitles: ['F$\uptot$',$
                                      'F$\upc1$',$
                                      'F$\upc2$',$
                                      'F$\down\\lambda$$\up-500km/s$'],$
                            vtitles: ['$\sigma$',$
                                      'v$\down50$$\upc1$',$
                                      'v$\down50$$\upc2$',$
                                      'v$\down98$$\upc2$'],$
                            radtitles: ['F$\uptot$']},$
                  aspectrat: (double(nrows)/double(ncols)),$
                  center_axes: [centcol,centrow],$
                  center_nuclei: [centcol,centrow],$
                  ctradprof_psffwhm: 0.6d,$
                  contourlevels: contourlevels,$f
                  ct: {sumrange: fitrange,$
                       sumrange_hstcomp: [3751,4800],$
                       scllim: [0,1],$
                       scllim_rad: [-3,0],$
                       stretch: 1,$
                       fitifspeak: 1b,$
                       fitifspeakwin_kpc: 3d},$
                  hst: {refcoords: [2825,2711],$
                        subim_sm: 7d,$
                        subim_big: 20d,$
                        smoothfwhm: 17},$
                  hstrd: {file: '/path/HSTimage.fits',$
                          label: 'ACS/F435W',$
                          scllim: [0.01,10],$
                          sclargs_sm: {beta: 0.05,stretch: 5},$
                          sclargs_big: {beta: 0.05,stretch: 5},$
                          sclargs_fov: {beta: 0.05,stretch: 5},$
                          platescale: 0.05d,$
                          nucoffset: [0d,0d]},$
                  fcnsortcomp: 'ifsf_sortcomp',$
                  rangefile: '/path/ranges.txt',$
                  fluxunits: 1d-15 $
                 }
   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Parameters for NaD + HeI 5876 fit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   initnad = {}

   return,init

end
