FUNCTION IFSF_pg0953ovi, directoryname, gal, zgal, profileshifts, profilesig, coveringfactor,$
  opticaldepth, _initmaps=initmaps, initnad=initnad


  ;readcol, initfile, gal, ncols, nrows, centcol, centrow, format = '(A,D,D,D,D,D,D,D)'
  ;bad=1d99
  bad = 1d99
  gal = gal[0]
  bin = 2d
  ncols = 1
  nrows = 1
  centcol = 1
  centrow = 1
  zgal=zgal[0]
  outstr = 'rb'+string(bin,format='(I0)')
  comps=N_ELEMENTS(profileshifts)
  ;  FOR K =0,(FILE_LINES('galaxyinitproc')-1) DO BEGIN
  ;  readcol,gal[K]+'n.txt', ignore, relativeflux, SKIPLINE=34
  ;  readcol,gal[K]+'.txt', wavelength, ignore, error
  ;  readcol,gal[K]+'f.txt', ignore, continuum, SKIPLINE=34
  readcol,directoryname+'/'+gal+'/'+gal+'.txt', wavelength, flux, error, FORMAT='D,D,D'

  ; Finding the index to fit over
  linefit=CALL_FUNCTION('ifsf_linelist',['[OVI1]1032','[OVI2]1038'])
  ;  checkindex1=lines['[NV1]1239']-50
  ;  checkindex2=lines['[NV2]1243']+50
  zgalint=DOUBLE(zgal)
  doubletregion=[1262,1283]
  checkindex=[1267,1281]
  plotrange=[VALUE_LOCATE(wavelength, checkindex[0]),VALUE_LOCATE(wavelength,checkindex[1])]
  fitrange=[VALUE_LOCATE(wavelength, doubletregion[0]),VALUE_LOCATE(wavelength,doubletregion[1])]
  Doubletabsorptionregion=[VALUE_LOCATE(wavelength,doubletregion[0]),VALUE_LOCATE(wavelength,doubletregion[1])]
  hydrogenabsorptionregion=[VALUE_LOCATE(wavelength,1215),VALUE_LOCATE(wavelength,1228)]
  spectramin=VALUE_LOCATE(wavelength,MIN(wavelength))
  spectramax=VALUE_LOCATE(wavelength,MAX(wavelength))
  indextoplot=[INDGEN(VALUE_LOCATE(wavelength,1259),START=spectramin), $
    INDGEN(VALUE_LOCATE(wavelength,1264)-VALUE_LOCATE(wavelength,1261),START=VALUE_LOCATE(wavelength,1261)),$
    INDGEN(VALUE_LOCATE(wavelength,1268.5)-VALUE_LOCATE(wavelength,1265.5),START=VALUE_LOCATE(wavelength,1265.5)),$
    INDGEN(VALUE_LOCATE(wavelength,1270.4)-VALUE_LOCATE(wavelength,1269.5),START=VALUE_LOCATE(wavelength,1269.5)),$
    INDGEN(VALUE_LOCATE(wavelength,1271.6)-VALUE_LOCATE(wavelength,1271),START=VALUE_LOCATE(wavelength,1271)),$
    INDGEN(VALUE_LOCATE(wavelength,1272.8)-VALUE_LOCATE(wavelength,1272),START=VALUE_LOCATE(wavelength,1272)),$
    INDGEN(VALUE_LOCATE(wavelength,1275.7)-VALUE_LOCATE(wavelength,1273.5),START=VALUE_LOCATE(wavelength,1273.5)),$
    INDGEN(VALUE_LOCATE(wavelength,1276.27)-VALUE_LOCATE(wavelength,1276.14),START=VALUE_LOCATE(wavelength,1276.14)),$
    INDGEN(VALUE_LOCATE(wavelength,1277)-VALUE_LOCATE(wavelength,1276.65),START=VALUE_LOCATE(wavelength,1276.65)),$
    INDGEN(N_ELEMENTS(wavelength)-VALUE_LOCATE(wavelength,1277.7),START=VALUE_LOCATE(wavelength,1279))]
  weight=1d/error^2
  fitreg=[[MIN(wavelength),1267],[1267,1280],[1280,MAX(wavelength)]]
  fitfcn=['ifsf_fitspline','ifsf_fitspline','ifsf_fitspline']
  fitargs=HASH()
  fitargs['reg1'] = {argsbkpts:{everyn:40}}
  fitargs['reg2'] = {argsbkpts:{everyn:40}}
  fitargs['reg3'] = {argsbkpts:{everyn:50}}

  cgplot, wavelength, flux, XRAN=checkindex, YRAN=[-.3*MAX(flux[plotrange[0]:plotrange[1]]),1.5*MAX(flux[plotrange[0]:plotrange[1]])],$
    XSTYLE=1,YSTYLE=1,backg='Black',axiscolor='White',color='White',$
    xtit='Wavelength ($\Angstrom$)',ytit='Flux (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)'
    continuum=CALL_FUNCTION('ifsf_fitmulticont', wavelength, flux, weight, ignored, indextoplot,0,fitreg=fitreg,$
    fitfcn=fitfcn, fitargs=fitargs)
  cgoplot, wavelength, continuum, color='Red',thick=4
  img = cgsnapshot(filename=directoryname+'/'+gal+'/'+gal+'_continuum',/jpeg,/nodialog,quality=100)

  relativeflux=flux/continuum

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
  infile=directoryname+'/'+gal+'/'+gal+'.txt'
  if ~ file_test(infile) then begin
    print,"ERROR: Data cube not found."
    return,0
  endif

  ; Lines to fit.
  lines = ['[OVI1]1032','[OVI2]1038',$
    '[LyB]1026','[LyA]1216','[NV1]1239','[NV2]1243']


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
  maxncomp = comps

  ; Initialize line ties, n_comps, z_inits, and sig_inits.
  linetie = hash(lines,'[LyB]1026')
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
  tweakcntfit[*,*,2,*] = 4
  ; Number of wavelength regions to re-fit
  nregions = 1
  ; Lower wavelength for re-fit
  tweakcntfit[*,*,0,0:nregions-1] = $
    rebin(reform(1031.613,1,1,1,nregions),$
    ncols,nrows,1,nregions)
  ; Upper wavelength for re-fit
  tweakcntfit[*,*,1,0:nregions-1] = $
    rebin(reform(1037.912,1,1,1,nregions),$
    ncols,nrows,1,nregions)
  ; Order for re-fit
  tweakcntfit[*,*,2,0:nregions-1] = $
    rebin(reform([3],1,1,1,nregions),$
    ncols,nrows,1,nregions)

  ; Parameters for emission line plotting
  ;  linoth = strarr(3,6)
  ;  linoth[0,2] = '[CII]1347'
  ;;  linoth[*,3] = '[NiII]1317'
  ;  argspltlin1 = {nx: 3, ny: 2,$
  ;    label: ['[LyA]1216','[NV]1239','[NV]1243',$
  ;    '','',''],$
  ;    wave: [1216, 1239, 1243, 0, 0, 0],$
  ;    off: [[-120,90],[-80,50],[-130,50],$
  ;    [-80,120],[-95,70],[-95,50]],$
  ;    linoth: linoth}
  ;  linoth = strarr(1,6)
  ;  linoth[*,0] = '[NiII]1317'
  ;  argspltlin2 = {nx: 3, ny: 2,$
  ;    label: ['[LyA]1216','[NV]1239','[NV]1243',$
  ;    '','',''],$
  ;    wave: [1216,1239,1243,0,0,0],$
  ;    off: [[-120,90],[-120,90],[-120,90],$
  ;    [-90,80],[-90,80],[-90,80]],$
  ;    linoth: linoth}


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
    mapdir: directoryname+'/'+gal+'/',$
    maxncomp: maxncomp,$
    ncomp: ncomp,$
    outdir: directoryname+'/'+gal+'/',$
    zsys_gas: zgal,$
    ; Optional pars
    ;        first # is max sig, second is step size
    startempfile: directoryname+'/'+gal+'/',$
    tweakcntfit: tweakcntfit $
  }


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Parameters for N V 1239 fit
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if keyword_set(initnad) then begin

    normnadlo = [1300,1330]
    normnadhi = [1310,1340]
    pltnormnad = [1300,1360]
    nad_maxncomp = comps

    ;     Initialize n_comps, z_inits, and sig_inits.
    ;     Use 1 HeI component w/in a circular region
    heitie = strarr(ncols,nrows)
    heitie[0,0]='[LyB]1026'
    heitiecol = intarr(ncols,nrows)+14
    heitierow = intarr(ncols,nrows)+14
    hei_zinit = dblarr(ncols,nrows,nad_maxncomp)
    hei_siginit = dblarr(ncols,nrows,nad_maxncomp)

    nnadabs = dblarr(ncols,nrows)
    nadabs_zinit = dblarr(ncols,nrows,nad_maxncomp) + zgal
    FOR I = 0, comps-1 DO BEGIN
      nadabs_zinit[*,*,I] += profileshifts[I]/1031.912d
    ENDFOR


    nadabs_siginit = dblarr(ncols,nrows,nad_maxncomp)
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
    nadabs_cfinit = coveringfactor
    nadabs_tauinit = opticaldepth

    ;     These parameters are from the average of spaxels [17,13], [17,14],
    ;     [16,14], and [16,15]. Note that one also has to set nadem_fitinit=0 for
    ;     this to work properly.
    ;      nadabs_fix[16,17,*,*]=1b
    ;      nadabs_cfinit[16,17,0:1]=[0.15d,0.25d]
    ;      nadabs_tauinit[16,17,0:1]=[1.48,0.1d]
    ;      nadabs_zinit[16,17,0:1]=[0.0424d,0.0414d]
    ;      nadabs_siginit[16,17,0:1]=[69d,288d]


    nadabs_siginit[0,0,0] = 75d

    nnadabs[*,0] =comps

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
      argsnadweq: {autowavelim: [1200,1250,1321,1325],$
      autoindices:1},$
      plotindex:plotrange,$
      fitindex:fitrange,$
      argsnormnad: {fitranlo: normnadlo,$
      fitranhi: normnadhi},$
      argspltnormnad: {fitranlo: normnadlo,$
      fitranhi: normnadhi,$
      pltran: pltnormnad},$
      argspltfitnad: {yran: [0,2]},$
      fcnfitnad: 'ifsf_uvabsfcnOVI',$
      fcninitpar: 'ifsf_inituvabsOVI',$
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
      error: error, $
      continuum: continuum, $
      flux: flux $
    }
  endif

  return,init
  ;  ENDFOR

END