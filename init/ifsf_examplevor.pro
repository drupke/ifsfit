; docformat = 'rst'
;
;+
;
; This function initializes the fitting parameters for PG1700
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
;      2019jan22, DSNR, created
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
function ifsf_examplevor,initmaps=initmaps,initnad=initnad

  bad=1d99

  gal = 'example'
  bin = 1d
  ncols = 15
  nrows = 23
  centcol = 8d
  centrow = 12d
  platescale = 0.3d
  outstr = 'vor'
  fitrange = [3750,6400]

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
  ncomp = hash(lines)
  linetie = hash(lines)
  zinit_gas = hash(lines)
  siginit_gas = hash(lines)
  zinit_stars=dblarr(ncols,nrows) + 0.1234d

  foreach i,lines do begin
     linetie[i] = 'Hbeta'
     ncomp[i] = dblarr(ncols,nrows)+maxncomp
     siginit_gas[i] = dblarr(ncols,nrows,maxncomp) + 100d
     zinit_gas[i] = dblarr(ncols,nrows,maxncomp) + 0.1234d
  endforeach

  vormap = mrdfits(infile,2,/silent)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Optional pars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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
         datext: -1,$
         varext: 1,$
         dqext: -1,$
         vormap: vormap,$
         argscheckcomp: {sigcut: 3d},$
         argspltlin1: argspltlin1,$
         fcncheckcomp: 'ifsf_checkcomp',$
         maskwidths_def: 500d,$
         tweakcntfit: tweakcntfit,$
         emlsigcut: 2d,$
         logfile: '/Users/drupke/specfits/kcwi/'+gal+'/'+outstr+'/'+$
                  gal+'_fitlog.txt',$
         batchfile: '/Users/drupke/Dropbox/git/ifsfit/common/ifsf_fitloop.pro',$
         batchdir: '/Users/drupke/src/idl/batch/',$
         siglim_gas: siglim_gas,$
         siginit_gas: siginit_gas,$
         siginit_stars: 50d, $
         cvdf_vlimits: [-2d3,2d3],$
         cvdf_vstep: 1d $
        }

  initmaps = {}
  initnad = {}
     
  return,init

end
