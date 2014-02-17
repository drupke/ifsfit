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
  bin = 1d

  outstr = 'rb'+string(bin,format='(I0)')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Required pars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Input file
  infile='/Users/drupke/winds/gmos/cubes/'+gal+'/'+gal+outstr+'.fits'
  if ~ file_test(infile) then begin
     print,"ERROR: Data cube not found."
     return,0
  endif

; Initialize line ties, n_comps, and z_inits.
  nlines = 2
  linetie = strarr(dx,dy,nlines) + 'Halpha'
  ncomp = dblarr(dx,dy,nlines) + 1
  zinit_gas=dblarr(dx,dy,nlines,maxncomp)
  zinit_stars=dblarr(dx,dy)

; Parameters for emission line plotting
  linoth = strarr(2,2)
  argspltlin1 = {nx: 2, ny: 1,$
                 label: ['[OI]6300','[OI]6363']
                 wave: [6300,6363]
                 off: [[-120,90],[-80,50]]
                 linoth: linoth}

; Velocity dispersion guesses
  siginit_gas = dblarr(nlines,3)+299792d/3000d/0.7d

; Velocity dispersion limits
  siglim_gas = [299792d/3000d/0.699d,299792d/3000d/0.701d]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output structure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  init = {$
; Required pars
         fcninitpar: 'uhsf_gm_initpar',$
         fitran: [6285d,6380d],$
         infile: infile,$
         linetie: linetie,$
         ncomp: ncomp,$
         outdir: '/Users/drupke/winds/gmos/specfits/'+$
         galname+'/o1_5577/exp'+exptag+'/',$
         zinit_stars: zinit_stars,$
         zinit_gas: zinit_gas,$
; Optional pars
         argsinitpar: {siglim: siglim_gas},$
         argspltlin1: argspltlin1,$
         fcncontfit: 'ifsf_fitcont',$
         siglim_gas: siglim_gas,$
         siginit_gas: siginit_gas,$
         }

  return,init

end

  galname = strmid(gal,0,6)
  exptag = strmid(gal,6,1)
  exptagmatch = {a: '1a', b: '1b', c: '1c', d: '1d', e: '3a', f: '3b', g: '4a'}
  tagnames = tag_names(exptagmatch)
  if galname eq 'f05189' AND exptag ne '' then begin
     zinit = dblarr(743)
     ncomp = intarr(743)+1
     itag = where(strcmp(exptag,tagnames,/fold_case),cttag)
     infile=rootindir+'red/'+galname+'/ctexrdat_'+exptagmatch.(itag)+'.fits'
     outlines = ['[OI]5577']
     fitran_rest=[5550,5600]
     siglim_gas = [0.65d,0.75]
     argsinitpar = {siglim_gas:siglim_gas}
     argstrlines = {outlines:outlines,argslinelist:argslinelist}
  endif

end
