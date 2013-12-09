; docformat = 'rst'
;
;+
;
; This function initializes the fitting parameters for F05189-2524,
; 2011 GMOS data.
;
; :Categories:
;    IFSF
;
; :Returns:
;    A structure with tags specified in INITTAGS.txt.
;
; :Params:
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
;      2013sep, DSNR, complete re-write
;      2013nov25, DSNR, renamed, added copyright and license; moved
;                       description of tags to INITTAGS.txt file.
;      2013nov26, DSNR, changed line arrays to hashes to prevent
;                       bookkeeping errors
;    
; :Copyright:
;    Copyright (C) 2013 David S. N. Rupke
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
function ifsf_f05189

  gal = 'f05189'
  bin = 2d

  dx = 28
  dy = 27
  cx = 14d
  cy = 14d 

  outstr = 'rb'+string(bin,format='(I0)')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Required pars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Input file
  infile='/Users/drupke/ifs/gmos/cubes/'+gal+'/'+gal+outstr+'.fits'
  if ~ file_test(infile) then begin
     print,"ERROR: Data cube not found."
     return,0
  endif

; Lines to fit.
  lines = ['Halpha','Hbeta','HeI6678','HeI7065','HeII4686',$
           '[OI]6300','[OI]6364','[OIII]4959','[OIII]5006',$
           '[NI]5198','[NI]5200','[NII]6548','[NII]6583',$
           '[SII]6716','[SII]6731','[FeVII]5159','[FeVII]5721',$
           '[FeVII]6087','[FeX]6375']
  nlines = n_elements(lines)

; Max no. of components.
  maxncomp = 3

; Initialize line ties, n_comps, and z_inits.
  linetie = hash(lines,strarr(nlines)+'Halpha')
  ncomp = hash()
  zinit_gas = hash()
  siginit_gas = hash() ; note that this is technically optional, put here for convenience
  foreach i,lines do begin
     ncomp[i] = dblarr(dx,dy)+maxncomp
     zinit_gas[i] = dblarr(dx,dy,maxncomp) + 0.0425d
     siginit_gas[i] = dblarr(maxncomp) + 150d
  endforeach
  zinit_stars=dblarr(dx,dy) + 0.043d
; iron lines
  tmplines = ['[FeVII]5159','[FeVII]5721','[FeVII]6087','[FeX]6375']
  foreach i,tmplines do begin
     linetie[i] = '[FeVII]6087'
     ncomp[i,*,*] = 1
     ncomp[i,13,13] = 2
     zinit_gas[i,13,13,0] = 0.038d
     zinit_gas[i,13,13,1] = 0.040d
     siginit_gas[i,0] = 1000d
  endforeach
; HeII line
  tmplines = ['HeII4686']
  foreach i,tmplines do begin
     linetie[i] = 'HeII4686'
     ncomp[i,*,*] = 2
     zinit_gas[i,*,*,0] = 0.040d
     zinit_gas[i,*,*,1] = 0.038d
     siginit_gas[i,1] = 1000d
  endforeach
; HeI lines
  tmplines = ['HeI6678','HeI7065']
  foreach i,tmplines do begin
     linetie[i] = 'HeI6678'
     ncomp[i,*,*] = 0
     ncomp[i,13,13] = 1
     siginit_gas[i,0] = 1000d
  endforeach
; Balmer lines, low-ion. colliosional lines
  tmplines = ['Halpha','Hbeta','[OI]6300','[OI]6364',$
              '[OIII]4959','[OIII]5006','[NII]6548','[NII]6583',$
              '[SII]6716','[SII]6731']
  foreach i,tmplines do begin
     zinit_gas[i,*,*,1] = 0.040d
     zinit_gas[i,*,*,2] = 0.038d
     siginit_gas[i,2] = 1000d
  endforeach
; [NI] lines
  tmplines = ['[NI]5198']
  foreach i,tmplines do begin
     linetie[i] = '[NI]5198'
     ncomp[i,*,*] = 1
     siginit_gas[i,2] = 1000d
  endforeach

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

; Velocity dispersion limits
  siglim_gas = [299792d/3000d/2.35d,2000d]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output structure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  init = {$
; Required pars
         fcninitpar: 'ifsf_gmos',$
         fitran: [4600,7100],$
         infile: infile,$
         lines: lines,$
         linetie: linetie,$
         maxncomp: maxncomp,$
         ncomp: ncomp,$
         outdir: '/Users/drupke/specfits/gmos/'+gal+'/'+outstr+'/',$
         zinit_stars: zinit_stars,$
         zinit_gas: zinit_gas,$
; Optional pars
         argscontfit: {refit: refit},$
         argsinitpar: {siglim: siglim_gas},$
         argslinelist: {felines: 1},$
         argsoptstelz: {lrange: [5200,5550]},$
         argspltlin1: argspltlin1,$
         argspltlin2: argspltlin2,$
         fcncontfit: 'ifsf_fitcont',$
         fcnoptstelsig: 'ifsf_optstelsig',$
         fcnoptstelz: 'ifsf_optstelz',$
         nomaskran: [5075,5100],$
         siglim_gas: siglim_gas,$
         siginit_gas: siginit_gas,$
         siginit_stars: 100d,$
;        first # is max sig, second is step size
         sigfitvals: dindgen(fix(500d/25d)+1)*25d,$
         startempfile: '/Users/drupke/Documents/stellar_models/'+$
         'gonzalezdelgado/SSPGeneva_z020.sav' $
         }

  return,init

end