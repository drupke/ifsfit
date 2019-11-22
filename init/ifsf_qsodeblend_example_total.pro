; docformat = 'rst'
;
;+
;
; This function initializes the fitting parameters for an example total spectrum.
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
;      2019nov22, DSNR, copied from ifsf_f13342host.pro
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
function ifsf_qsodeblend_example_total,initmaps=initmaps,initnad=initnad

  bad=1d99

  gal = 'galaxyshorthand'
  nrows = 1
  ncols = 1
  fitran = [5400,8200]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Required pars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; Input file
  infile='/path/totalspectrum.fits'
  if ~ file_test(infile) then message,'Data cube not found.'

; Lines to fit.
  lines = ['Halpha','Hbeta',$
           '[OI]6300','[OI]6364','[OIII]4959','[OIII]5007',$
           '[NII]6548','[NII]6583','[SII]6716','[SII]6731',$
           '[NI]5198','[NI]5200']
  nlines = n_elements(lines)

; Max no. of components.
  maxncomp = 3

; Initialize line ties, n_comps, z_inits, and sig_inits.
  linetie = hash(lines,'Halpha')
  ncomp = hash(lines)
  zinit_gas = hash(lines)
  siginit_gas = hash(lines)
; note that siginit_gas is technically optional, put here for convenience
  foreach i,lines do begin
     ncomp[i] = dblarr(ncols)+maxncomp
     zinit_gas[i] = dblarr(ncols,maxncomp) + 0.1797d
     siginit_gas[i] = dblarr(ncols,maxncomp) + 75d
  endforeach
  zinit_stars=dblarr(ncols) + 0.1797d



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Optional pars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Tweaked regions are around Hbeta/[OIII], HeI/NaD, [OI], Ha/[NII], and [SII]

; Parameters for continuum fit
  tweakcntfit = dblarr(ncols,nrows,3,5)
; Default fitting order
  tweakcntfit[*,*,2,*] = 2
; Number of wavelength regions to re-fit
  nregions = 5
; Lower wavelength for re-fit
  tweakcntfit[*,*,0,0:nregions-1] = $
     rebin(reform([5650,6850,7400,7600,7850],1,1,1,nregions),$
           ncols,nrows,1,nregions)
; Upper wavelength for re-fit
  tweakcntfit[*,*,1,0:nregions-1] = $
     rebin(reform([6000,7050,7550,7850,8000],1,1,1,nregions),$
           ncols,nrows,1,nregions)
; Order for re-fit
  tweakcntfit[*,*,2,0:nregions-1] = $
     rebin(reform([1,1,1,1,1],1,1,1,nregions),$
           ncols,nrows,1,nregions)


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
  siglim_gas = [5d,1000d]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output structure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  init = {$
;        Required pars
         fcninitpar: 'ifsf_gmos',$
         fitran: fitrange,$
         infile: infile,$
         label: gal,$
         lines: lines,$
         linetie: linetie,$ ;hash
         maxncomp: maxncomp,$
         name: 'galaxyname',$
         ncomp: ncomp,$ ;hash
         outdir: '/path/fitdir/'+gal+'/total/',$
         zinit_stars: zinit_stars,$ ;hash
         zinit_gas: zinit_gas,$ ;hash
         zsys_gas: 0.1797d,$
;        Optional pars
         argscheckcomp: {sigcut: 3d,$
                         ignore: ['[OI]6300','[OI]6364',$
                                  '[SII]6716','[SII]6731']},$ ;struct
         argspltlin1: argspltlin1,$
         decompose_qso_fit: 1b,$
         donad: 1b,$
         emlsigcut: 2d,$
         fcncheckcomp: 'ifsf_checkcomp',$
         fcncontfit: 'ifsf_fitqsohost',$
         argscontfit: {qsoxdr: '/path/nucleartemplate.xdr',$
                       siginit_stars: 50d,$
                       uselog: 1b,$
                       refit: 1b,$
                       add_poly_degree: 50d,$
                       print_output: 1b},$ ;struct
         tweakcntfit: tweakcntfit, $
         startempfile: '/path/stellar_models/'+$
                       'gonzalezdelgado/SSPGeneva_z020.sav',$
         maskwidths_def: 500d,$
         siglim_gas: siglim_gas,$
         siginit_gas: siginit_gas, $ ; hash
         cutrange: [[6920,6965]], $
         host: {dat_fits: '/path/cubes/'+gal+'/'+gal+'starlight_fromtotal.fits', $
                singlespec: 1b} $
        }

   initmaps = {}
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Parameters for NaD + HeI 5876 fit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if keyword_set(initnad) then begin

      normnadlo = [6810,6910]
      normnadhi = [6960,7060]
      pltnormnad = [6810,7060]
      maxncomp = 2

;     Initialize n_comps, z_inits, and sig_inits.
      heitie = strarr(ncols,nrows) + 'HeI5876'
      hei_zinit = dblarr(ncols,nrows,maxncomp)+0.1797d
      hei_siginit = dblarr(ncols,nrows,maxncomp)+10d
      hei_zinit[*,*,0] = 0.1793d
      hei_zinit[*,*,1] = 0.1801d

      nnadabs = intarr(ncols,nrows) + 1
      nadabs_zinit = dblarr(ncols,nrows,maxncomp)+0.1797d
      nadabs_siginit = dblarr(ncols,nrows,maxncomp)+50d
      nadabs_siglim = [5d,1000d]
      nadabs_cfinit = dblarr(ncols,nrows,maxncomp)+0.5d
      nadabs_tauinit = dblarr(ncols,nrows,maxncomp)+0.5d

      initnad = {$
                 argsinitpar: {siglimhei: [5d,1000d]},$
                 argsnadweq: {autowavelim: [6930,6960,6960,7000],$
                              autoindices:1},$
                 argsnormnad: {fitranlo: normnadlo,$
                               fitranhi: normnadhi,$
                               snavg_thresh: 1d},$
                 argspltnormnad: {fitranlo: normnadlo,$
                                  fitranhi: normnadhi,$
                                  pltran: pltnormnad,$
                                  fitord: 3},$
                 argspltfitnad: {yran: [0,2]},$
                 fcnfitnad: 'ifsf_nadfcn',$
                 fcninitpar: 'ifsf_initnad',$
                 maxncomp: maxncomp,$
                 mcniter: 1000, $
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
                 nhei: dblarr(ncols,nrows)+2d,$
                 hei_zinit: hei_zinit,$
                 hei_siginit: hei_siginit,$
                 heitie: heitie $
                }
   endif
                                    
   return,init

end
