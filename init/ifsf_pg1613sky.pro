; docformat = 'rst'
;
;+
;
; This function initializes the fitting parameters for the sky lines in 
; F13342, 2012 GMOS data.
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
;      2014jun09, DSNR, created
;
; :Copyright:
;    Copyright (C) 2014 David S. N. Rupke
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
function ifsf_pg1613sky,skyexp=skyexp

  ncols=749
  gal='pg1613'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Required pars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  infile='/Users/drupke/ifs/gmos/red/'+gal+'/ctexrdat_'+skyexp+'.fits'
  lines = ['[OI]5577']
  maxncomp = 1
  linetie = hash(lines,'[OI]5577')
  ncomp = hash(lines)
  zinit_gas = hash(lines)
  siginit_gas = hash(lines)
  foreach i,lines do begin
    ncomp[i] = dblarr(ncols)+maxncomp
    zinit_gas[i] = dblarr(ncols,maxncomp)
    siginit_gas[i] = dblarr(maxncomp) + 40d
  endforeach
  zinit_stars=dblarr(ncols)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Optional pars
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; Parameters for emission line plotting
  argspltlin1 = {nx: 1, ny: 1,$
                 label: ['[OI]5577'],$
                 wave: [5577],$
                 off: [[-20,20]]}

  ; Velocity dispersion limits
  siglim_gas = [35,50]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Output structure
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  init = {$
    ; Required pars
    fcninitpar: 'ifsf_gmos',$
    fitran: [5550,5600],$
    infile: infile,$
    label: gal,$
    lines: lines,$
    linetie: linetie,$
    maxncomp: maxncomp,$
    ncomp: ncomp,$
    outdir: '/Users/drupke/specfits/gmos/'+gal+'/o1_5577/exp'+skyexp+'/',$
    zinit_stars: zinit_stars,$
    zinit_gas: zinit_gas,$
    ; Optional pars
    argscontfit: {nobvls: 1, fitord: 1},$
    argsinitpar: {siglim: siglim_gas},$
    argspltlin1: argspltlin1,$
    datext: 2,$
    dqext: 4,$
    fcncontfit: 'ifsf_fitpoly',$
    siglim_gas: siglim_gas,$
    siginit_gas: siginit_gas, $
    varext: 3 $
  }

  return,init

end
