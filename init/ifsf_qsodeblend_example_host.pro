; docformat = 'rst'
;
;+
;
; This function initializes the fitting parameters for an example host spectrum.
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
function ifsf_qsodeblend_example_host,initmaps=initmaps,initnad=initnad

  bad=1d99

  gal = 'galaxyshorthand'
  ncols=1
  nrows=1
  fitran = [5400,8200]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Required pars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Input file
; Comment out for step 10.
  infile='/path/starlight_spatiallyintegrated.fits'
; Uncomment for step 10
;  infile='/path/starlight_spatiallyintegrated_iter1.fits'
  if ~ file_test(infile) then message,'Data cube not found.'

; Lines to fit.
  lines = ['Halpha']
  nlines = n_elements(lines)

; Max no. of components.
  maxncomp = 0

; Initialize line ties, n_comps, z_inits, and sig_inits.
  linetie = hash(lines,'Halpha')
  ncomp = hash(lines)
  zinit_gas = hash(lines)
  siginit_gas = hash(lines)
; note that siginit_gas is technically optional, put here for convenience
  foreach i,lines do ncomp[i] = dblarr(ncols,nrows)+maxncomp
  zinit_stars=dblarr(ncols,nrows) + 0.1797d

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output structure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  init = {$
; Required pars
         fcninitpar: 'ifsf_gmos',$
         fitran: fitran,$
         infile: infile,$
         label: gal,$
         name: 'galaxyname',$
         lines: lines,$
         linetie: linetie,$
         maxncomp: maxncomp,$
         ncomp: ncomp,$
         outdir: '/path/fitdir/'+gal+'/host/',$
         zinit_stars: zinit_stars,$
         zinit_gas: zinit_gas,$
         zsys_gas: 0.1797,$
; Optional pars
         donad: 1b,$
         noemlinfit: 1,$
; Regions to ignore: Hb/[OIII], HeI/NaD
         cutrange: [[5800,5950],[6920,6965]], $
         siginit_stars: 50d,$
         startempfile: '/Users/drupke/Documents/stellar_models/'+$
                       'gonzalezdelgado/SSPGeneva_z020.sav', $
         fcncontfit: 'ppxf', $
         decompose_ppxf_fit: 1b,$
         argscontfit: {print_output: 1b, add_poly_degree: 50} $
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
                 zref: 0.1797d,$
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
