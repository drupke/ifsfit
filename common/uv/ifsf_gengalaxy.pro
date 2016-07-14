FUNCTION IFSF_GENGALAXYNbgV, directoryname, gal, zgal, profileshifts, profilesig, _initmaps=initmaps, initnad=initnad


;readcol, initfile, gal, ncols, nrows, centcol, centrow, format = '(A,D,D,D,D,D,D,D)'
;bad=1d99
  bad = 1d99
  gal = gal
  bin = 2d
  ncols = 1
  nrows = 1
  centcol = 1
  centrow = 1
  outstr = 'rb'+string(bin,format='(I0)')
;  FOR K =0,(FILE_LINES('galaxyinitproc')-1) DO BEGIN
;  readcol,gal[K]+'n.txt', ignore, relativeflux, SKIPLINE=34
;  readcol,gal[K]+'.txt', wavelength, ignore, error
;  readcol,gal[K]+'f.txt', ignore, continuum, SKIPLINE=34
  readcol,directoryname+gal+'n.txt', ignore, relativeflux, SKIPLINE=34
  readcol,directoryname+gal+'.txt', wavelength, ignore, error
  readcol,directoryname+gal+'f.txt', ignore, continuum, SKIPLINE=34

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
;  infile='/Users/to/ToIRAF/spectra/galaxies/'+gal[K]+'.txt'
  infile=directoryname+gal+'.txt'
  if ~ file_test(infile) then begin
    print,"ERROR: Data cube not found."
    return,0
  endif

  ; Lines to fit.
  lines = ['Halpha','Hbeta','HeI6678','HeI7065','HeII4686',$
    '[OI]6300','[OI]6364','[OIII]4959','[OIII]5007',$
    '[NI]5198','[NI]5200','[NII]6548','[NII]6583',$
    '[SII]6716','[SII]6731','[FeVII]5159','[FeVII]5721',$
    '[FeVII]6087','[FeX]6375','[OVI1]1032','[OVI2]1038',$
    '[LyB]1026','[LyA]1216','[NV1]1239','[NV2]1243', $
    '[CII]1347', '[CII]1335.7', '[CII]1335.6', '[CII]1335', $
    '[CI]1329', '[NiII]1317']
    

;  lines = hash()
;  waves = hash()
;  lineshorts = hash()
;
;  FOR I = 0, (N_ELEMENTS(linename)-1) DO BEGIN
;    waves[linename[I]] = linewave[I]
;    lines[linename[I]] = linename[I]
;    lineshorts[linename[I]] = lineshort[I]
;  ENDFOR
  
  winit = hash(lines)
  signit = hash(lines)
  cfinit = hash(lines)
  tauinit = hash(lines)
  nlines = n_elements(lines)

  ; Max no. of components.
  maxncomp = 4

  ; Initialize line ties, n_comps, z_inits, and sig_inits.
  linetie = hash(lines,'[CII]1347')
  ncomp = hash(lines)
  zinit_gas = hash(lines)
  siginit_gas = hash(lines)
  ; note that siginit_gas is technically optional, put here for convenience
;  foreach j,lines do begin
;    ncomp[j] = dblarr(ncols,nrows)+maxncomp
;    zinit_gas[j] = dblarr(ncols,nrows,maxncomp) + zgal[K]
;    siginit_gas[j] = dblarr(maxncomp)
;  endforeach
  foreach i,lines do begin
  ncomp[i] = dblarr(ncols,nrows)+maxncomp
  zinit_gas[i] = dblarr(ncols,nrows,maxncomp) + zgal
  siginit_gas[i] = dblarr(maxncomp)
  endforeach
;  zinit_stars=dblarr(ncols,nrows) + zgal[K]
  zinit_stars=dblarr(ncols,nrows) + zgal
;  linetie = hash(lines,'[NV1]1239')
;  ncomp = hash(lines)
;  zinit_gas = hash(lines)
;  siginit_gas = hash(lines)
;  ; note that siginit_gas is technically optional, put here for convenience
;  foreach i,lines do begin
;    ncomp[i] = dblarr(ncols,nrows)+maxncomp
;    zinit_gas[i] = dblarr(ncols,nrows,maxncomp) + 0.0882d
;    siginit_gas[i] = dblarr(maxncomp) + 150d
;  endforeach
;  zinit_stars=dblarr(ncols,nrows) + 0.0882d
;  linetie = hash(lines,'[NV2]1243')
;  ncomp = hash(lines)
;  zinit_gas = hash(lines)
;  siginit_gas = hash(lines)
;  ; note that siginit_gas is technically optional, put here for convenience
;  foreach i,lines do begin
;    ncomp[i] = dblarr(ncols,nrows)+maxncomp
;    zinit_gas[i] = dblarr(ncols,nrows,maxncomp) + 0.0882d
;    siginit_gas[i] = dblarr(maxncomp) + 150d
;  endforeach
  zinit_stars=dblarr(ncols,nrows) + zgal
  ; Carbon
;  tmplines = ['[NV1]1239','[NV2]1243','[CII]1335','[CII]1335.6','[CII]1347', '[CII]1335.7']
;  foreach i,tmplines do begin
;    linetie[i] = '[LyA]1216'
;    ncomp[i,*,*] = 2
;    zinit_gas[i,*,*,0] = 0.041d
;    siginit_gas[i,0] = 1000d
;    ;     if ctnuc0 gt 0 then for j=0,ctnuc0-1 do begin
;    ;        ncomp[i,x_pix[inuc0[j]]-1,y_pix[inuc0[j]]-1] = 2
;    ;        zinit_gas[i,x_pix[inuc0[j]]-1,y_pix[inuc0[j]]-1,1] = 0.038d
;    ;     endfor
;    if ctedge0 gt 0 then for j=0,ctedge0-1 do $
;      ncomp[i,x_pix[iedge0[j]]-1,y_pix[iedge0[j]]-1] = 0
;  endforeach


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Optional pars
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; Parameters for continuum fit
  tweakcntfit = dblarr(ncols,nrows,3,10)
  ; Default fitting order
  tweakcntfit[*,*,2,*] = 2
  ; Number of wavelength regions to re-fit
  nregions = 1
  ; Lower wavelength for re-fit
  tweakcntfit[*,*,0,0:nregions-1] = $
    rebin(reform([1100],1,1,1,nregions),$
    ncols,nrows,1,nregions)
  ; Upper wavelength for re-fit
  tweakcntfit[*,*,1,0:nregions-1] = $
    rebin(reform([1350],1,1,1,nregions),$
    ncols,nrows,1,nregions)
  ; Order for re-fit
  tweakcntfit[*,*,2,0:nregions-1] = $
    rebin(reform([3],1,1,1,nregions),$
    ncols,nrows,1,nregions)

  ; Parameters for emission line plotting
  linoth = strarr(3,6)
  linoth[0,2] = '[CII]1347'
;  linoth[*,3] = '[NiII]1317'
  argspltlin1 = {nx: 3, ny: 2,$
    label: ['[LyA]1216','[NV]1239','[NV]1243',$
    '','',''],$
    wave: [1216, 1239, 1243, 0, 0, 0],$
    off: [[-120,90],[-80,50],[-130,50],$
    [-80,120],[-95,70],[-95,50]],$
    linoth: linoth}
  linoth = strarr(1,6)
  linoth[*,0] = '[NiII]1317'
  argspltlin2 = {nx: 3, ny: 2,$
    label: ['[LyA]1216','[NV]1239','[NV]1243',$
    '','',''],$
    wave: [1216,1239,1243,0,0,0],$
    off: [[-120,90],[-120,90],[-120,90],$
    [-90,80],[-90,80],[-90,80]],$
    linoth: linoth}

  ; Velocity dispersion limits and fixed values
  siglim_gas = [6d,2000d]
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
    fcninitpar: 'ifsf_gmos',$$
    infile: infile,$
    label: gal,$
    lines: lines,$
    linetie: linetie,$
    mapdir: 'files/',$
    maxncomp: maxncomp,$
    ncomp: ncomp,$
    outdir: 'files/',$
    platescale: 0.11d,$
    specres: 1.6d,$
    positionangle: 0d,$
    zinit_stars: zinit_stars,$
    zinit_gas: zinit_gas,$
    zsys_gas: 0.0882d,$
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
    nomaskran: [1100,1300],$
    siglim_gas: siglim_gas,$
    siginit_gas: siginit_gas,$
    siginit_stars: 100d,$
    ;        first # is max sig, second is step size
    startempfile: 'files/',$
    tweakcntfit: tweakcntfit $
  }

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Arguments for maps
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  if keyword_set(initmaps) then begin
;    argslinratmaps_comp = hash()
;    argslinratmaps_comp['lrat1'] = [['1_n2ha','2_n2ha'],$
;      ['1_o3hb','2_o3hb'],$
;      ['1_n2ha_vs_o3hb','2_n2ha_vs_o3hb']]
;    argslinratmaps_comp['ebv'] = ['1_ebv','2_ebv']
;    argslinratmaps_cvdf = hash()
;    argslinratmaps_cvdf['lrat1'] = $
;      [['ftot_n2ha','fpk_n2ha','fv50_n2ha','fv84_n2ha','fv98_n2ha'],$
;      ['ftot_o3hb','fpk_o3hb','fv50_o3hb','fv84_o3hb','fv98_o3hb'],$
;      ['ftot_n2ha_vs_o3hb','fpk_n2ha_vs_o3hb','fv50_n2ha_vs_o3hb',$
;      'fv84_n2ha_vs_o3hb','fv98_n2ha_vs_o3hb']]
;    argslinratmaps_cvdf['ebv'] = $
;      ['ftot_ebv','fpk_ebv','fv50_ebv','fv84_ebv','fv98_ebv']
;
;    badnademp = bytarr(ncols,nrows)
;    badnademp[0,*]=1b
;    badnademp[*,nrows-1]=1b
;    badnademp[5:8,1]=1b
;    badnademp[1:2,2]=1b
;
;    initmaps = {$
;      aspectrat: 1.05d,$
;      center_axes: [centcol,centrow],$
;      center_nuclei: [centcol,centrow],$
;      rangefile: '/Users/drupke/ifs/gmos/maps/'+$
;      'f05189/rb2/ranges.txt',$
;      argslinratmaps_comp: argslinratmaps_comp,$
;      argslinratmaps_cvdf: argslinratmaps_cvdf,$
;      badnademp: badnademp,$
;      doemlinradprof: 1,$
;      emlinradprof_psffwhm: 0.6d,$
;      ;                 Base units are 10^-15 erg/s/cm^2/spaxel
;      ;                 Multiplying fluxes by 10 gives fluxes in units of 10^-16 erg/s/cm^2/spaxel
;      ;                 Dividing fluxes by 0.04 gives fluxes in units of 10^-16 erg/s/cm^2/arcsecond
;      ;                 15jan26 -- DSNR -- oops! Was multiplying by 20 instead of 25.
;      fluxfactor: 10d*25d,$
;      ;                  applyebv: [1,0,0],$
;      nadabsweq_snrthresh: 3d,$
;      nademweq_snrthresh: 3d,$
;      nademflux_cbint: 0.5d,$
;      fcn_oplots: 'ifsf_makemaps_f05189',$
;      tags_oplots: ['nadcube',$
;      'nadfit',$
;      'initnad',$
;      'nadabsncomp',$
;      'map_rkpc_hst',$
;      'map_rkpc_bhst',$
;      'map_rkpc_rhst',$
;      'bhst_fov_ns',$
;      'rhst_fov_ns',$
;      'bhst_big',$
;      'rhst_big',$
;      'hst_big_ifsfov',$
;      'cshst_fov_s',$
;      'chst_fov_ns',$
;      'cshst_fov_ns',$
;      'cshst_fov_rb',$
;      'ctcube',$
;      'contcube',$
;      'nadabsnh','errnadabsnh',$
;      'nadabscnh','errnadabscnh',$
;      'nadabssig','nademsig',$
;      'nadabsvel','nademvel',$
;      'nadabsv98','nademv98',$
;      'errnadabsvel','errnademvel',$
;      'nadabscf','errnadabscf',$
;      'nadabstau','errnadabstau'],$
;      col: {sumrange: [1300,1330,1310,1340],$
;      scllim: [-0.1,0.2],$
;      stretch: 1},$
;      ct: {sumrange: [1000,1400],$
;      scllim: [0,1],$
;      stretch: 1},$
;      ; This coordinate is in zero-offset pixels; i.e., the central pixel as measured in
;      ; DS9 minus 1 (for PA=0; could be different for other PAs). It is chosen to
;      ; align the two continuum maps in *cont.eps, and to center
;      ; the HST map for plotting. The nuclear offsets below give the nuclear
;      ; coordinates of the red and blue maps for making galactocentric radius arrays,
;      ; also in single-offset pixels.
;      hst: {refcoords: [3261,2708],$
;      subim_sm: 7d,$
;      subim_big: 25d,$
;      smoothfwhm: 12},$
;      hstbl: {file: '/Users/drupke/ifs/gmos/ancillary/hst/'+$
;      'f05189/f05189_acs_435w.fits',$
;      scllim: [0.01,100],$
;      sclargs_sm: {beta: 0.05,stretch: 5},$
;      sclargs_big: {beta: 0.05,stretch: 5},$
;      sclargs_fov: {beta: 1,stretch: 5},$
;      photflam: 3.1840413d-19,$
;      photplam: 4318.8672,$
;      platescale: 0.05d,$
;      nucoffset: [-0.5d,0d]},$
;      hstblsm: {scllim: [0,10],$
;      sclargs: {beta: 0.5, stretch: 5}},$
;      hstrd: {file: '/Users/drupke/ifs/gmos/ancillary/hst/'+$
;      'f05189/f05189_acs_814w.fits',$
;      scllim: [0.01,100],$
;      sclargs_sm: {beta: 0.05,stretch: 5},$
;      sclargs_big: {beta: 0.05,stretch: 5},$
;      sclargs_fov: {beta: 1,stretch: 5},$
;      photflam: 7.0331898e-20,$
;      photplam: 8056.944,$
;      platescale: 0.05d,$
;      nucoffset: [-0.75d,-0.75d]},$
;      hstrdsm: {scllim: [0,20],$
;      sclargs: {beta: 0.5,stretch: 5}}, $
;      hstcol: {scllim: [0.5,2.5],$
;      sclargs: {stretch: 1},$
;      ncbdiv: 4},$
;      hstcolsm: {scllim: [0.8,1.8],$
;      sclargs: {stretch: 1},$
;      ncbdiv: 5}$
;  }
;endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Parameters for NaD + HeI 5876 fit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if keyword_set(initnad) then begin

  normnadlo = [1300,1330]
  normnadhi = [1310,1340]
  pltnormnad = [1300,1340]
  nad_maxncomp = 20

  ;     Initialize n_comps, z_inits, and sig_inits.
  ;     Use 1 HeI component w/in a circular region
  heitie = strarr(ncols,nrows)
  heitie[0,0]='[LyA]1216'
  heitiecol = intarr(ncols,nrows)+14
  heitierow = intarr(ncols,nrows)+14
  hei_zinit = dblarr(ncols,nrows,nad_maxncomp)
  hei_siginit = dblarr(ncols,nrows,nad_maxncomp)

  nnadabs = dblarr(ncols,nrows)
;  nadabs_zinit = dblarr(ncols,nrows,nad_maxncomp)+zgal[K]
;  nadabs_zinit[*,*,0] += (6/1238.8210)
;  nadabs_zinit[*,*,1] += (5/1238.8210)
;  nadabs_zinit[*,*,2] += (3.2/1238.8210)
;  nadabs_zinit[*,*,3] += (1.6/1238.8210)
  nadabs_zinit = dblarr(ncols,nrows,nad_maxncomp)+zgal
  FOR I = 0, N_ELEMENTS(profilesig)-1 DO BEGIN
    nadabs_zinit[*,*,I] += profileshifts[I]/1238.8210
  ENDFOR
;  nadabs_zinit[*,*,0] += .00484261501d
;  nadabs_zinit[*,*,3] += .00395480225d
;  nadabs_zinit[*,*,4] += .002582728d
;  nadabs_zinit[*,*,7] += .001291364d
  
  
  nadabs_siginit = dblarr(ncols,nrows,nad_maxncomp)+100d
  ;PROPOSED
;  nadabs_siginit[*,*,0] = 10d
;  nadabs_siginit[*,*,1] = 100d
;  nadabs_siginit[*,*,2] = 10d
;  nadabs_siginit[*,*,3] = 100d
  ;PROPOSED
  
  FOR I = 0, N_ELEMENTS(profilesig)-1 DO BEGIN
    nadabs_siginit[*,*,I] += profilesig[I]
  ENDFOR
;  nadabs_siginit[*,*,0] = 10d  
;  nadabs_siginit[*,*,1] = 100d
;  nadabs_siginit[*,*,2] = 10d
;  nadabs_siginit[*,*,3] = 100d
  
  nadabs_siglim = [6d,1000d]
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


  nadabs_siginit[0,0,0] = 75d

  nnadabs[*,0] = 8

  nnadem = dblarr(ncols,nrows)
  nadem_zinit = dblarr(ncols,nrows,nad_maxncomp)+zgal
  ;      nadem_siginit = dblarr(ncols,nrows,nad_maxncomp)+150d
  nadem_siginit = dblarr(ncols,nrows,nad_maxncomp)+75d
  nadem_finit = dblarr(ncols,nrows,nad_maxncomp)+0.1d
  nadem_rinit = dblarr(ncols,nrows,nad_maxncomp)+1.5d
;  nadem_siglim = [299792d/3000d/2.35d,750d]
  nadem_siglim = [2d, 750d]
  nadem_fix = bytarr(ncols,nrows,nad_maxncomp,4)

  ;      nadem_rinit[*,*,*] = 2.005d
  nadem_rinit[*,*,*] = 1d
  nadem_fix[*,*,*,3] = 1b

  ;     Regions where we don't fix the emission line ratio
  nadem_fix[0,0,*,3] = 0b

  nadem_fix[0,0,*,3] = 0b
  ;      nadem_fix[13,23,*,3] = 0b
  ;      nadem_fix[15,24,*,3] = 0b
  ;      nadem_fix[17,22,*,3] = 0b


  nadem_zinit[0,0,0]=zgal


  ;      nnadem[0,9:15]=1
  nnadem[0,0]=0

  ;      nnadem[26,15:24]=1

  ;     Trying to understand upper limits for errors due to
  ;     absorption/emission mixing
  ;      nadem_zinit[16,*,0]=0.043d
  ;      nadem_fix[16,*,0,0]=1b
  ;      nnadabs[16,*]=1

  initnad = {$
    argsnadweq: {autowavelim: [1315,1320,1321,1325],$
    autoindices:1},$
    argsnormnad: {fitranlo: normnadlo,$
    fitranhi: normnadhi},$
    argspltnormnad: {fitranlo: normnadlo,$
    fitranhi: normnadhi,$
    pltran: pltnormnad},$
    argspltfitnad: {yran: [0,2]},$
    fcnfitnad: 'ifsf_uvabsfcn',$
    fcninitpar: 'ifsf_inituvabs',$
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
    heitie: heitie, $
    galaxy: gal, $
    wavelength: wavelength, $
    relativeflux: relativeflux, $
    error: error $
  }
endif

return,init
;  ENDFOR

END