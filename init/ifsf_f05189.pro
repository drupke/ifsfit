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
;      2013sep, DSNR, complete re-write
;      2013nov25, DSNR, renamed, added copyright and license; moved
;                       description of tags to INITTAGS.txt file.
;      2013nov26, DSNR, changed line arrays to hashes to prevent
;                       bookkeeping errors
;      2013dec10, DSNR, testing and bug fixes
;      2013dec17, DSNR, renamed variables dx, dy, cx, cy;
;                       moved from unordered to ordered hashes; 
;                       turn hashes into structures before passing to IFSF
;      2013jan13, DSNR, updated to pass hashes for many parameters into IFSF, 
;                       instead of structures
;      2014jan16, DSNR, fixed one wrong wavelength label
;      2014feb26, DSNR, replaced ordered hashes with hashes
;      2014feb27, DSNR, added zsys_gas, platescale, specres
;      2014apr21, DSNR, added arguments for line ratio maps / VO plots
;      2014may23, DSNR, added arguments for plotting continuum images;
;                       changed way that map making and NaD parameters are 
;                       treated
;    
; :Copyright:
;    Copyright (C) 2013-2014 David S. N. Rupke
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
function ifsf_f05189,initmaps=initmaps,initnad=initnad

  bad=1d99

  gal = 'f05189'
  bin = 2d
  ncols = 28
  nrows = 27
  centcol = 14
  centrow = 14
  outstr = 'rb'+string(bin,format='(I0)')

; distance from central pixel
  x_pix = rebin(indgen(ncols)+1,ncols,nrows)
  y_pix = rebin(transpose(indgen(nrows)+1),ncols,nrows)
  rad_pix = sqrt((double(x_pix-centcol))^2d + (double(y_pix-centrow))^2d)

; Regions for setting components
  inuc0  = where(rad_pix le 3d,ctnuc0)
  inuc1  = where(rad_pix gt 3d AND rad_pix le 6d,ctnuc1)
  idisk0 = where(rad_pix gt 8d,ctdisk0)
  iedge0 = where(rad_pix ge 12d,ctedge0)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Required pars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Input file
;  infile='/Users/drupke/ifs/gmos/cubes/'+gal+'/'+gal+outstr+'.fits'
  infile='/Users/drupke/ifs/gmos/cubes/'+gal+'/rb2_7exp/'+gal+outstr+'.fits'
  if ~ file_test(infile) then begin
     print,"ERROR: Data cube not found."
     return,0
  endif

; Lines to fit.
  lines = ['Halpha','Hbeta','HeI6678','HeI7065','HeII4686',$
           '[OI]6300','[OI]6364','[OIII]4959','[OIII]5007',$
           '[NI]5198','[NI]5200','[NII]6548','[NII]6583',$
           '[SII]6716','[SII]6731','[FeVII]5159','[FeVII]5721',$
           '[FeVII]6087','[FeX]6375']
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
     ncomp[i] = dblarr(ncols,nrows)+maxncomp
     zinit_gas[i] = dblarr(ncols,nrows,maxncomp) + 0.0425d
     siginit_gas[i] = dblarr(maxncomp) + 150d
  endforeach
  zinit_stars=dblarr(ncols,nrows) + 0.043d
; iron lines
  tmplines = ['[FeVII]5159','[FeVII]5721','[FeVII]6087','[FeX]6375']
  foreach i,tmplines do begin
     linetie[i] = '[FeVII]6087'
     ncomp[i,*,*] = 1
     zinit_gas[i,*,*,0] = 0.041d
     siginit_gas[i,0] = 1000d
;     if ctnuc0 gt 0 then for j=0,ctnuc0-1 do begin
;        ncomp[i,x_pix[inuc0[j]]-1,y_pix[inuc0[j]]-1] = 2
;        zinit_gas[i,x_pix[inuc0[j]]-1,y_pix[inuc0[j]]-1,1] = 0.038d
;     endfor
     if ctedge0 gt 0 then for j=0,ctedge0-1 do $
        ncomp[i,x_pix[iedge0[j]]-1,y_pix[iedge0[j]]-1] = 0
  endforeach
; HeII line
  tmplines = ['HeII4686']
  foreach i,tmplines do begin
     linetie[i] = 'HeII4686'
     ncomp[i,*,*] = 1
     zinit_gas[i,*,*,0] = 0.041d
     if ctnuc0 gt 0 then for j=0,ctnuc0-1 do begin
        ncomp[i,x_pix[inuc0[j]]-1,y_pix[inuc0[j]]-1] = 2
        zinit_gas[i,x_pix[inuc0[j]]-1,y_pix[inuc0[j]]-1,1] = 0.039d
        siginit_gas[i,1] = 1000d
     endfor
;     if ctnuc1 gt 0 then for j=0,ctnuc1-1 do $
;       ncomp[i,x_pix[inuc1[j]]-1,y_pix[inuc1[j]]-1] = 1
;     if ctedge0 gt 0 then for j=0,ctedge0-1 do $
;        ncomp[i,x_pix[iedge0[j]]-1,y_pix[iedge0[j]]-1] = 0
  endforeach
; HeI lines
  tmplines = ['HeI6678','HeI7065']
  foreach i,tmplines do begin
     linetie[i] = 'HeI6678'
     ncomp[i,*,*] = 0
     if ctnuc0 gt 0 then for j=0,ctnuc0-1 do begin
        ncomp[i,x_pix[inuc0[j]]-1,y_pix[inuc0[j]]-1] = 1
        siginit_gas[i,0] = 500d
     endfor
  endforeach
; [OIII] lines
  tmplines = ['[OIII]4959','[OIII]5007']
  foreach i,tmplines do begin
    ncomp[i,*,*] = 2
    linetie[i] = '[OIII]5007'
    zinit_gas[i,*,*,0] = 0.041d
    zinit_gas[i,*,*,1] = 0.039d
    siginit_gas[i,1] = 1000d
  endforeach
; Balmer lines, low-ion. colliosional lines
  tmplines = ['Halpha','Hbeta',$
              '[NII]6548','[NII]6583',$
              '[SII]6716','[SII]6731']
  foreach i,tmplines do begin
     zinit_gas[i,*,*,1] = 0.041d
     zinit_gas[i,*,*,2] = 0.039d
     siginit_gas[i,2] = 1000d
;     if ctedge0 gt 0 then for j=0,ctedge0-1 do $
;        ncomp[i,x_pix[iedge0[j]]-1,y_pix[iedge0[j]]-1] = 1
     if ctdisk0 gt 0 then for j=0,ctdisk0-1 do $
        ncomp[i,x_pix[idisk0[j]]-1,y_pix[idisk0[j]]-1] = 2
  endforeach
;; Note that if we want to allow Hbeta to vary independently of Halpha, we also
;; have to turn off the line ratio constraint in IFSF_GMOS.
;; Hbeta
;  tmplines = ['Hbeta']
;  foreach i,tmplines do begin
;     linetie[i] = 'Hbeta'
;     zinit_gas[i,*,*,1] = 0.041d
;     zinit_gas[i,*,*,2] = 0.039d
;     siginit_gas[i,2] = 1000d
;     if ctedge0 gt 0 then for j=0,ctedge0-1 do $
;        ncomp[i,x_pix[iedge0[j]]-1,y_pix[iedge0[j]]-1] = 1
;     if ctdisk0 gt 0 then for j=0,ctdisk0-1 do $
;        ncomp[i,x_pix[idisk0[j]]-1,y_pix[idisk0[j]]-1] = 2
;  endforeach
; [OI] lines
  tmplines = ['[OI]6300','[OI]6364']
  foreach i,tmplines do begin
     linetie[i] = '[OI]6300'
     ncomp[i,*,*] = 2
     zinit_gas[i,*,*,1] = 0.041d
     if ctdisk0 gt 0 then for j=0,ctdisk0-1 do $
        ncomp[i,x_pix[idisk0[j]]-1,y_pix[idisk0[j]]-1] = 1
  endforeach
; [NI] lines
  tmplines = ['[NI]5198','[NI]5200']
  foreach i,tmplines do begin
     linetie[i] = '[NI]5200'
     ncomp[i,*,*] = 1
;     if ctedge0 gt 0 then for j=0,ctedge0-1 do $
;        ncomp[i,x_pix[iedge0[j]]-1,y_pix[iedge0[j]]-1] = 0
  endforeach

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Optional pars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Parameters for continuum fit
  tweakcntfit = dblarr(ncols,nrows,3,10)
; Default fitting order
  tweakcntfit[*,*,2,*] = 2
; Number of wavelength regions to re-fit
  nregions = 7
; Lower wavelength for re-fit
  tweakcntfit[*,*,0,0:nregions-1] = $
     rebin(reform([4950,5250,5850,6200,6500,6725,6925],1,1,1,nregions),$
           ncols,nrows,1,nregions)
; Upper wavelength for re-fit
  tweakcntfit[*,*,1,0:nregions-1] = $
     rebin(reform([5100,5450,6000,6400,6700,6925,7100],1,1,1,nregions),$
           ncols,nrows,1,nregions)
; Order for re-fit
  tweakcntfit[*,*,2,0:nregions-1] = $
     rebin(reform([2,2,2,2,2,1,1],1,1,1,nregions),$
           ncols,nrows,1,nregions)

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

; Velocity dispersion limits and fixed values
  siglim_gas = [299792d/3000d/2.35d,2000d]
  sigfix=hash()
  sigfix['[FeVII]6087'] = 725d
  lratfix=hash()
  lratfix['[NI]5200/5198'] = [1.5d]
  lratfix['[NII]6583/Ha'] = [bad,1.80,2.14]
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output structure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  init = {$
; Required pars
         fcninitpar: 'ifsf_gmos',$
         fitran: [4600,7100],$
         infile: infile,$
         label: gal,$
         lines: lines,$
         linetie: linetie,$
         mapdir: '/Users/drupke/ifs/gmos/maps/'+gal+'/'+outstr+'/',$
         maxncomp: maxncomp,$
         ncomp: ncomp,$
         outdir: '/Users/drupke/specfits/gmos/'+gal+'/'+outstr+'/',$
         platescale: 0.2d,$
         specres: 1.6d,$
         positionangle: 0d,$
         zinit_stars: zinit_stars,$
         zinit_gas: zinit_gas,$
         zsys_gas: 0.04275d,$
; Optional pars
;         argscheckcomp: {sigcut: 2},$
         argsinitpar: {siglim: siglim_gas,$
                       sigfix: sigfix,$
                       lratfix: lratfix},$
         argspltlin1: argspltlin1,$
         argspltlin2: argspltlin2,$
         donad: 1,$
;         dored: 1,$
         fcncheckcomp: 'ifsf_checkcomp',$
         fcncontfit: 'ppxf',$
         mapcent: [centcol,centrow],$
         nomaskran: [5075,5100],$
         siglim_gas: siglim_gas,$
         siginit_gas: siginit_gas,$
         siginit_stars: 100d,$
;        first # is max sig, second is step size
         startempfile: '/Users/drupke/Documents/stellar_models/'+$
         'gonzalezdelgado/SSPGeneva_z020.sav', $
         tweakcntfit: tweakcntfit $
        }

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Arguments for maps
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
   if keyword_set(initmaps) then begin
      argslinratmaps_comp = hash()
      argslinratmaps_comp['lrat1'] = [['1_n2ha','2_n2ha'],$
                                      ['1_o3hb','2_o3hb'],$
                                      ['1_n2ha_vs_o3hb','2_n2ha_vs_o3hb']]
      argslinratmaps_comp['ebv'] = ['1_ebv','2_ebv']
      argslinratmaps_cvdf = hash()
      argslinratmaps_cvdf['lrat1'] = $
         [['ftot_n2ha','fpk_n2ha','fv50_n2ha','fv84_n2ha','fv98_n2ha'],$
          ['ftot_o3hb','fpk_o3hb','fv50_o3hb','fv84_o3hb','fv98_o3hb'],$
          ['ftot_n2ha_vs_o3hb','fpk_n2ha_vs_o3hb','fv50_n2ha_vs_o3hb',$
           'fv84_n2ha_vs_o3hb','fv98_n2ha_vs_o3hb']]
      argslinratmaps_cvdf['ebv'] = $
         ['ftot_ebv','fpk_ebv','fv50_ebv','fv84_ebv','fv98_ebv']

      badnademp = bytarr(ncols,nrows)
      badnademp[0,*]=1b
      badnademp[*,nrows-1]=1b
      badnademp[5:8,1]=1b
      badnademp[1:2,2]=1b

      initmaps = {$
                  aspectrat: 1.05d,$
                  center_axes: [centcol,centrow],$
                  center_nuclei: [centcol,centrow],$
                  rangefile: '/Users/drupke/ifs/gmos/maps/'+$
                             'f05189/rb2/ranges.txt',$
                  argslinratmaps_comp: argslinratmaps_comp,$
                  argslinratmaps_cvdf: argslinratmaps_cvdf,$
                  badnademp: badnademp,$
                  doemlinradprof: 1,$
                  emlinradprof_psffwhm: 0.6d,$
;                 Base units are 10^-15 erg/s/cm^2/spaxel
;                 Multiplying fluxes by 10 gives fluxes in units of 10^-16 erg/s/cm^2/spaxel
;                 Dividing fluxes by 0.04 gives fluxes in units of 10^-16 erg/s/cm^2/arcsecond
;                 15jan26 -- DSNR -- oops! Was multiplying by 20 instead of 25.
                  fluxfactor: 10d*25d,$
;                  applyebv: [1,0,0],$
                  nadabsweq_snrthresh: 3d,$
                  nademweq_snrthresh: 3d,$
                  nademflux_cbint: 0.5d,$
                  fcn_oplots: 'ifsf_makemaps_f05189',$
                  tags_oplots: ['nadcube',$
                                'nadfit',$
                                'initnad',$
                                'nadabsncomp',$
                                'map_rkpc_hst',$
                                'map_rkpc_bhst',$
                                'map_rkpc_rhst',$
                                'bhst_fov_ns',$
                                'rhst_fov_ns',$
                                'bhst_big',$
                                'rhst_big',$
                                'hst_big_ifsfov',$
                                'cshst_fov_s',$
                                'chst_fov_ns',$
                                'cshst_fov_ns',$
                                'cshst_fov_rb',$
                                'ctcube',$
                                'contcube',$
                                'nadabsnh','errnadabsnh',$
                                'nadabscnh','errnadabscnh',$
                                'nadabssig','nademsig',$
                                'nadabsvel','nademvel',$
                                'nadabsv98','nademv98',$
                                'errnadabsvel','errnademvel',$
                                'nadabscf','errnadabscf',$
                                'nadabstau','errnadabstau'],$
                  col: {sumrange: [4900,5000,6650,6750],$
                        scllim: [-0.1,0.2],$
                        stretch: 1},$
                  ct: {sumrange: [5600,6400],$
                       scllim: [0,1],$
                       stretch: 1},$
; This coordinate is in zero-offset pixels; i.e., the central pixel as measured in
; DS9 minus 1 (for PA=0; could be different for other PAs). It is chosen to 
; align the two continuum maps in *cont.eps, and to center 
; the HST map for plotting. The nuclear offsets below give the nuclear 
; coordinates of the red and blue maps for making galactocentric radius arrays, 
; also in single-offset pixels.
                  hst: {refcoords: [3261,2708],$
                        subim_sm: 7d,$
                        subim_big: 25d,$
                        smoothfwhm: 12},$
                  hstbl: {file: '/Users/drupke/ifs/gmos/ancillary/hst/'+$
                                'f05189/f05189_acs_435w.fits',$
                          scllim: [0.01,100],$
                          sclargs_sm: {beta: 0.05,stretch: 5},$
                          sclargs_big: {beta: 0.05,stretch: 5},$
                          sclargs_fov: {beta: 1,stretch: 5},$
                          photflam: 3.1840413d-19,$
                          photplam: 4318.8672,$
                          platescale: 0.05d,$
                          nucoffset: [-0.5d,0d]},$
                  hstblsm: {scllim: [0,10],$
                            sclargs: {beta: 0.5, stretch: 5}},$
                  hstrd: {file: '/Users/drupke/ifs/gmos/ancillary/hst/'+$
                                'f05189/f05189_acs_814w.fits',$
                          scllim: [0.01,100],$
                          sclargs_sm: {beta: 0.05,stretch: 5},$
                          sclargs_big: {beta: 0.05,stretch: 5},$
                          sclargs_fov: {beta: 1,stretch: 5},$
                          photflam: 7.0331898e-20,$
                          photplam: 8056.944,$
                          platescale: 0.05d,$
                          nucoffset: [-0.75d,-0.75d]},$
                  hstrdsm: {scllim: [0,20],$
                            sclargs: {beta: 0.5,stretch: 5}}, $
                  hstcol: {scllim: [0.5,2.5],$
                           sclargs: {stretch: 1},$
                           ncbdiv: 4},$
                  hstcolsm: {scllim: [0.8,1.8],$
                             sclargs: {stretch: 1},$
                             ncbdiv: 5}$
                 }
   endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Parameters for NaD + HeI 5876 fit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if keyword_set(initnad) then begin
  
      normnadlo = [6040,6090]
      normnadhi = [6170,6220]
      pltnormnad = [6040,6220]
      nad_maxncomp = 2

;     Initialize n_comps, z_inits, and sig_inits.
;     Use 1 HeI component w/in a circular region
      heitie = strarr(ncols,nrows)
      heitie[11,13:16]='HeI6678'
      heitie[12,12:17]='HeI6678'
      heitie[13,11:18]='HeI6678'
      heitie[14,11:18]='HeI6678'
      heitie[15,12:18]='HeI6678'
      heitie[16,12:16]='HeI6678'
      heitie[17,14:15]='HeI6678'
      heitiecol = intarr(ncols,nrows)+14
      heitierow = intarr(ncols,nrows)+14
      hei_zinit = dblarr(ncols,nrows,nad_maxncomp)
      hei_siginit = dblarr(ncols,nrows,nad_maxncomp)

      nnadabs = dblarr(ncols,nrows)
      nadabs_zinit = dblarr(ncols,nrows,nad_maxncomp)+0.042
      nadabs_zinit[*,*,1] = 0.041
      nadabs_siginit = dblarr(ncols,nrows,nad_maxncomp)+100d
      nadabs_siginit[*,*,1] = 300d
      nadabs_siglim = [299792d/3000d/2.35d,1000d]
      nadabs_fix = bytarr(ncols,nrows,nad_maxncomp,4)
      nadabs_cfinit = dblarr(ncols,nrows,nad_maxncomp)+0.5d
      nadabs_tauinit = dblarr(ncols,nrows,nad_maxncomp)+0.5d

;     These parameters are from the average of spaxels [17,13], [17,14],
;     [16,14], and [16,15]. Note that one also has to set nadem_fitinit=0 for 
;     this to work properly.
;      nadabs_fix[16,17,*,*]=1b
;      nadabs_cfinit[16,17,0:1]=[0.15d,0.25d]
;      nadabs_tauinit[16,17,0:1]=[1.48,0.1d]
;      nadabs_zinit[16,17,0:1]=[0.0424d,0.0414d]
;      nadabs_siginit[16,17,0:1]=[69d,288d]

      nadabs_zinit[15,8,0] = 0.0425d
      nadabs_siginit[15,8,0] = 75d      
      nadabs_zinit[17,10:12,0] = 0.0425d
      nadabs_siginit[17,10:12,0] = 75d
      nadabs_zinit[18,9:13,0] = 0.0425d
      nadabs_siginit[18,9:13,0] = 75d

      nnadabs[3,8:17] = 1
      nnadabs[4,6:11] = 1
      nnadabs[4,12:15] = 2
      nnadabs[4,16:17] = 1
      nnadabs[5,5:10] = 1
      nnadabs[5,11:17] = 2
      nnadabs[5,18:19] = 1
      nnadabs[6,5:7] = 1
      nnadabs[6,8:17] = 2
      nnadabs[6,18:19] = 1
      nnadabs[7,5:7] = 1
      nnadabs[7,8:18] = 2
      nnadabs[7,19:20] = 1
      nnadabs[8,4:5] = 1
      nnadabs[8,6:19] = 2
      nnadabs[8,20] = 1
      nnadabs[9,3:5] = 1
      nnadabs[9,6:19] = 2
      nnadabs[9,20:21] = 1
      nnadabs[10,3:5] = 1
      nnadabs[10,6:19] = 2
      nnadabs[10,20:21] = 1
      nnadabs[11,3:5] = 1
      nnadabs[11,6:19] = 2
      nnadabs[11,20:21] = 1
      nnadabs[12,4:19] = 2
      nnadabs[12,20:21] = 1      
      nnadabs[13,4] = 1
      nnadabs[13,5:18] = 2
      nnadabs[13,19:20] = 1
      nnadabs[14,3:6] = 1
      nnadabs[14,7:18] = 2
      nnadabs[14,19] = 1
      nnadabs[15,3:5] = 1
      nnadabs[15,6:18] = 2
      nnadabs[15,19] = 1
      nnadabs[16,2:5] = 1
      nnadabs[16,6:17] = 2
      nnadabs[16,18:19] = 1
      nnadabs[17,1:6] = 1
      nnadabs[17,7:16] = 2
      nnadabs[17,17:18] = 1
      nnadabs[18,2:8] = 1
      nnadabs[18,9:13] = 2
      nnadabs[18,14:18] = 1
      nnadabs[19,6:16] = 1
      nnadabs[20,6:15] = 1
      nnadabs[21,8:13] = 1

      nnadem = dblarr(ncols,nrows)
      nadem_zinit = dblarr(ncols,nrows,nad_maxncomp)+0.044d
;      nadem_siginit = dblarr(ncols,nrows,nad_maxncomp)+150d
      nadem_siginit = dblarr(ncols,nrows,nad_maxncomp)+75d
      nadem_finit = dblarr(ncols,nrows,nad_maxncomp)+0.1d
      nadem_rinit = dblarr(ncols,nrows,nad_maxncomp)+1.5d
      nadem_siglim = [299792d/3000d/2.35d,750d]
      nadem_fix = bytarr(ncols,nrows,nad_maxncomp,4)

;      nadem_rinit[*,*,*] = 2.005d
      nadem_rinit[*,*,*] = 1d
      nadem_fix[*,*,*,3] = 1b

;     Regions where we don't fix the emission line ratio
      nadem_fix[5:6,2,*,3] = 0b
      nadem_fix[6:8,1,*,3] = 0b

      nadem_fix[3,22:25,*,3] = 0b
      nadem_fix[4,20:24,*,3] = 0b
      nadem_fix[5,21:24,*,3] = 0b
      nadem_fix[6,21:24,*,3] = 0b
      nadem_fix[7,22:25,*,3] = 0b
      nadem_fix[8,22:23,*,3] = 0b
      nadem_fix[8:10,25,*,3] = 0b

;      nadem_fix[13,23,*,3] = 0b
;      nadem_fix[15,24,*,3] = 0b
;      nadem_fix[17,22,*,3] = 0b
      nadem_fix[19,23,*,3] = 0b
      nadem_fix[20,20:23,*,3] = 0b
      nadem_fix[21,21:22,*,3] = 0b
      nadem_fix[22,15:17,*,3] = 0b
      nadem_fix[22,19:20,*,3] = 0b
      nadem_fix[22,22,*,3] = 0b
      nadem_fix[23,17:18,*,3] = 0b
      nadem_fix[24,16:21,*,3] = 0b
      nadem_fix[20:24,24,*,3] = 0b

      nadem_zinit[0:8,0:2,0]=0.043d
      nadem_zinit[0:4,0:5,0]=0.043d
      nadem_zinit[0:8,19:nrows-1,0]=0.043d
      nadem_zinit[9:11,24:nrows-1,0]=0.043d
      nadem_zinit[12:14,23:nrows-1,0]=0.043d
      nadem_zinit[15:17,20:nrows-1,0]=0.043d
      nadem_zinit[18:22,18:nrows-1,0]=0.043d
      nadem_zinit[23:26,*,0]=0.043d

;      nnadem[0,9:15]=1
      nnadem[1,12:17]=1
      nnadem[2,9:19]=1
      nnadem[3:9,0:25]=1
      nnadem[4:8,1]=0
      nnadem[10:15,0:6]=1
      nnadem[14,6]=0
      nnadem[10:15,22:25]=1
      nnadem[13,20:21]=1
;      nnadem[14,7]=1
      nnadem[14,18:21]=1
      nnadem[15,15:21]=1
      nnadem[16,0:2]=1
      nnadem[16,4:5]=1
      nnadem[16,13]=1
      nnadem[16,15:25]=1
      nnadem[17,0:5]=1
      nnadem[17,13:25]=1
      nnadem[18,0:5]=1
      nnadem[18,12:25]=1
      nnadem[19,0:4]=1
      nnadem[19,12:25]=1
      nnadem[20,0:5]=1
      nnadem[20,11:25]=1
      nnadem[21,2:7]=1
      nnadem[21,11:25]=1
      nnadem[22,2:4]=1
      nnadem[22,11:25]=1
      nnadem[23,3:25]=1
      nnadem[24,5:8]=1
      nnadem[24,12:25]=1
      nnadem[25,16:22]=1
;      nnadem[26,15:24]=1

;     Trying to understand upper limits for errors due to 
;     absorption/emission mixing
;      nadem_zinit[16,*,0]=0.043d
;      nadem_fix[16,*,0,0]=1b
;      nnadabs[16,*]=1

      initnad = {$
                 argsnadweq: {autowavelim: [6110,6160,6140,6180],$
                              autoindices:1},$
                 argsnormnad: {fitranlo: normnadlo,$
                               fitranhi: normnadhi},$
                 argspltnormnad: {fitranlo: normnadlo,$
                                  fitranhi: normnadhi,$
                                  pltran: pltnormnad},$
                 argspltfitnad: {yran: [0,2]},$
                 fcnfitnad: 'ifsf_nadfcn',$
                 fcninitpar: 'ifsf_initnad',$
                 maxncomp: nad_maxncomp,$
                 mcniter: 1000,$
;                 outxdr: '/Users/drupke/specfits/gmos/f05189/rb2/'+$
;                          'f05189.nadfit_emoptthin.xdr',$
;                 outxdr: '/Users/drupke/specfits/gmos/f05189/rb2/'+$
;                          'f05189.nadfit_emoptthick.xdr',$
;                NaD absorption
                 nnadabs: nnadabs,$
                 nadabs_cfinit: nadabs_cfinit,$
                 nadabs_tauinit: nadabs_tauinit,$
                 nadabs_zinit: nadabs_zinit,$
                 nadabs_siginit: nadabs_siginit,$
                 nadabs_siglim: nadabs_siglim,$
                 nadabs_fix: nadabs_fix,$
;                NaD emission
                 nnadem: nnadem,$
                 nadem_fitinit: 1,$
                 nadem_zinit: nadem_zinit,$
                 nadem_siginit: nadem_siginit,$
                 nadem_siglim: nadem_siglim,$
                 nadem_finit: nadem_finit,$
                 nadem_rinit: nadem_rinit,$
                 nadem_fix: nadem_fix,$
;                HeI
                 hei_zinit: hei_zinit,$
                 hei_siginit: hei_siginit,$
                 heitiecol: heitiecol,$
                 heitierow: heitierow,$
                 heitie: heitie $
                }
   endif
                  
   return,init

end
