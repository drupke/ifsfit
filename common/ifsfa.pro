; docformat = 'rst'
;
;+
;
; This procedure is the core routine to plot the continuum and emission
; lines fits to a spectrum.
;
; As input, it requires a structure of initialization parameters and the output 
; XDR file from IFSF. The tags for the initialization structure can be found in 
; INITTAGS.txt.
;
; Output plots:
;
; Output data files:
;
; Continuum data file: initdat.outdir+initdat.label+'.cont.xdr'
; Contains structure:
;   CONTCUBE
; 
; NaD data file: file=initdat.outdir+initdat.label+'.cont.xdr'
; Contains structure:
;   NADCUBE
; 
; Emission line data file: initdat.outdir+initdat.label+'.lin.xdr'
; Contains hashes:
;   EMLWAV[compkeys,lines,dx,dy]
;   EMLWAVERR[compkeys,lines,dx,dy]
;   EMLSIG[compkeys,lines,dx,dy]
;   EMLSIGERR[compkeys,lines,dx,dy]
;   EMLFLX[fluxkeys,lines,dx,dy]
;   EMLFLXERR[fluxkeys,lines,dx,dy]
;   EMLWEQ[fluxkeys,lines,dx,dy]
;   EMLCVDF[cvdfkeys]
; where
;   lines = emission line names [Halpha, Hbeta, ...]
;   compkeys = 'c1','c2',...,'cN' [N = max. no. of components]
;   fluxkeys = 'ftot' [total flux in line]
;            = 'fc1', 'fc2', ... [total flux in each component]
;            = 'fc1pk', 'fc2pk', ... [peak flux in each component]
;   cvdfkeys = 'vel' [M-element array of velocities]
;            = 'flux' (hash [lines,dx,dy,M] -- flux density profile)
;            = 'fluxerr' (hash [lines,dx,dy,M] -- flux density error profile)
;            = 'cumfluxnorm' (hash [lines,dx,dy,M] -- actual CVDF)
;            
; :Categories:
;    IFSFIT
;
; :Returns:
;    None.
;
; :Params:
;    initproc: in, required, type=string
;      Name of procedure to initialize the fit.
;
; :Keywords:
;    cols: in, optional, type=intarr, default=all
;      Columns to fit, in 1-offset format. Either a scalar or a
;      two-element vector listing the first and last columns to fit.
;    rows: in, optional, type=intarr, default=all
;      Rows to fit, in 1-offset format. Either a scalar or a
;      two-element vector listing the first and last rows to fit.
;    noplots: in, optional, type=byte
;      Disable plotting.
;    oned: in, optional, type=byte
;      Data is assumed to be in a 2d array; choose this switch to
;      input data as a 1d array.
;    verbose: in, optional, type=byte
;      Print error and progress messages. Propagates to most/all
;      subroutines.
; 
; :Author:
;    David Rupke
;
; :History:
;    Change History::
;      2009may13, DSNR, created
;      2013oct04, DSNR, started re-write for new data
;      2013oct09, DSNR, documented
;      2013nov25, DSNR, renamed, added copyright and license; changed
;                       required parameters from 'gal' and 'bin' to
;                       'initproc', and optional parameter 'fibers' to
;                       'oned', to make it more general
;      2014jan13, DSNR, propagated use of hashes; 
;                       updated for new linelist routine; 
;                       re-wrote IFSF_PRINTLINPAR and created IFSF_PRINTFITPAR
;      2014jan16, DSNR, bugfixes
;      2014jan29, DSNR, added _extra parameter to permit passing parameters
;                       to initialization routine; added some lines to deal
;                       properly with case of 1d data "cube"
;      2014feb26, DSNR, replaced ordered hashes with hashes
;      2014apr15, DSNR, bugfix in file sanity check
;      2014mayXY, DSNR, added NaD outputs
;      2014jun04, DSNR, added extra dimension to output LINMAP hash to get
;                       peak flux in emission lines
;      2014jul29, DSNR, updated treatment of emission line upper limits
;      2014decXY, DSNR, added output of details from SPS fit and calculation
;                       of relative contributions of SPS models
;                       and polynomial components to continuum fit
;      2015jan05, DSNR, added calculation and output of NaI D equivalent width
;                       from SPS fit
;      2015jun03, DSNR, updated for QSO/host decomposition in continuum and
;                       NaD fits
;      2015sep20, DSNR, output line fluxes/errors summed over components
;      2016jan06, DSNR, allow no emission line fit with initdat.noemlinfit
;      2016feb02, DSNR, handle cases with QSO+stellar continuum fits
;      2016sep02, DSNR, modified decomposition of QSO/host so that
;                       continuum tweaks are split between the two models;
;                       still need to implement for QSO+stellar case
;      2016sep13, DSNR, added internal logic to check if emission-line fit present
;      2016sep22, DSNR, adjusted quasar/host decomposition to account for new
;                       quasar fitting options
;      2016oct10, DSNR, added option to combine doublet fluxes
;      2018feb08, DSNR, updated call to IFSF_READCUBE
;      2018feb23, DSNR, write host spectrum with wavelength extension
;      2018jun26, DSNR, added option to re-weight NaD spectra by stellar chi-squared
;      2019jan24, DSNR, added option to pass arguments to continuum plotting
;      2019mar20, DSNR, compute Weq
;
; :Copyright:
;    Copyright (C) 2013--2019 David S. N. Rupke
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
pro ifsfa,initproc,cols=cols,rows=rows,noplots=noplots,oned=oned,$
          verbose=verbose,_extra=_extra

  bad = 1d99
  fwhmtosig = 2d*sqrt(2d*alog(2d))
  if keyword_set(verbose) then quiet=0 else quiet=1
  if keyword_set(oned) then oned=1 else oned=0

; Get fit initialization
  if keyword_set(_extra) then $
     initdat = call_function(initproc,_extra=_extra) $
  else $
     initdat = call_function(initproc)   
  if tag_exist(initdat,'donad') then begin
    initnad={dumy: 1}
    if keyword_set(_extra) then $
      initdat = call_function(initproc,initnad=initnad,_extra=_extra) $
    else $
      initdat = call_function(initproc,initnad=initnad)
  endif

  if ~ tag_exist(initdat,'noemlinfit') then begin
;    Get linelist
     if tag_exist(initdat,'argslinelist') then $
        linelist = ifsf_linelist(initdat.lines,_extra=initdat.argslinelist) $
     else linelist = ifsf_linelist(initdat.lines)
     nlines = linelist.count()
;    Linelist with doublets to combine
     emldoublets = [['[SII]6716','[SII]6731'],$
                    ['[OII]3726','[OII]3729'],$
                    ['[NI]5198','[NI]5200'],$
                    ['[NeIII]3869','[NeIII]3967'],$
                    ['[NeV]3345','[NeV]3426'],$
                    ['MgII2796','MgII2803']]
     sdoub = size(emldoublets)
     if sdoub[0] eq 1 then ndoublets = 1 else ndoublets = sdoub[2]
     lines_with_doublets = initdat.lines
     for i=0,ndoublets-1 do begin
        if linelist.haskey(emldoublets[0,i]) AND $
           linelist.haskey(emldoublets[1,i]) then begin
           dkey = emldoublets[0,i]+'+'+emldoublets[1,i]
           lines_with_doublets = [lines_with_doublets,dkey]
        endif
     endfor
     if tag_exist(initdat,'argslinelist') then $
        linelist_with_doublets = $
           ifsf_linelist(lines_with_doublets,_extra=initdat.argslinelist) $
     else  linelist_with_doublets = $
           ifsf_linelist(lines_with_doublets)
  endif

  if tag_exist(initdat,'fcnpltcont') then $
     fcnpltcont=initdat.fcnpltcont $
  else fcnpltcont='ifsf_pltcont'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if not tag_exist(initdat,'datext') then datext=1 else datext=initdat.datext
  if not tag_exist(initdat,'varext') then varext=2 else varext=initdat.varext
  if not tag_exist(initdat,'dqext') then dqext=3 else dqext=initdat.dqext

  header=1b
  if tag_exist(initdat,'argsreadcube') then $
     cube = ifsf_readcube(initdat.infile,quiet=quiet,oned=oned,header=header,$
                          datext=datext,varext=varext,dqext=dqext,$
                          _extra=initdat.argsreadcube) $
  else $
     cube = ifsf_readcube(initdat.infile,quiet=quiet,oned=oned,header=header,$
                          datext=datext,varext=varext,dqext=dqext)                     

  if tag_exist(initdat,'vormap') then begin
     vormap=initdat.vormap
     nvorcols = max(vormap)
     vorcoords = intarr(nvorcols,2)
     for i=1,nvorcols do begin
        ivor = where(vormap eq i,ctivor)
        xyvor = array_indices(vormap,ivor[0])
        vorcoords[i-1,*] = xyvor
     endfor
  endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Initialize output files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if ~ tag_exist(initdat,'noemlinfit') then $
     ifsf_printlinpar,lines_with_doublets,linlun,$
                      outfile=initdat.outdir+initdat.label+'.lin.dat'
  ifsf_printfitpar,fitlun,$
                   outfile=initdat.outdir+initdat.label+'.fit.dat'
                   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Initialize line hash
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if ~ tag_exist(initdat,'noemlinfit') then begin
      emlwav = hash()
      emlwaverr = hash()
      emlsig = hash()
      emlsigerr = hash()
      emlweq = hash()
      emlflx = hash()
      emlflxerr = hash()
      emlweq['ftot'] = hash()
      emlflx['ftot'] = hash()
      emlflxerr['ftot'] = hash()
      for k=0,initdat.maxncomp-1 do begin
         cstr='c'+string(k+1,format='(I0)')
         emlwav[cstr]=hash()
         emlwaverr[cstr]=hash()
         emlsig[cstr]=hash()
         emlsigerr[cstr]=hash()
         emlweq['f'+cstr]=hash()
         emlflx['f'+cstr]=hash()
         emlflxerr['f'+cstr]=hash()
         emlflx['f'+cstr+'pk']=hash()
         emlflxerr['f'+cstr+'pk']=hash()
      endfor
      foreach line,lines_with_doublets do begin
         emlweq['ftot',line]=dblarr(cube.ncols,cube.nrows)+bad
         emlflx['ftot',line]=dblarr(cube.ncols,cube.nrows)+bad
         emlflxerr['ftot',line]=dblarr(cube.ncols,cube.nrows)+bad
         for k=0,initdat.maxncomp-1 do begin
            cstr='c'+string(k+1,format='(I0)')
            emlwav[cstr,line]=dblarr(cube.ncols,cube.nrows)+bad
            emlwaverr[cstr,line]=dblarr(cube.ncols,cube.nrows)+bad
            emlsig[cstr,line]=dblarr(cube.ncols,cube.nrows)+bad
            emlsigerr[cstr,line]=dblarr(cube.ncols,cube.nrows)+bad
            emlweq['f'+cstr,line]=dblarr(cube.ncols,cube.nrows)+bad
            emlflx['f'+cstr,line]=dblarr(cube.ncols,cube.nrows)+bad
            emlflxerr['f'+cstr,line]=dblarr(cube.ncols,cube.nrows)+bad
            emlflx['f'+cstr+'pk',line]=dblarr(cube.ncols,cube.nrows)+bad
            emlflxerr['f'+cstr+'pk',line]=dblarr(cube.ncols,cube.nrows)+bad
         endfor
      endforeach
   endif

   if tag_exist(initdat,'flipsort') then begin
     flipsort = bytarr(cube.ncols,cube.nrows)
     sizefs = size(initdat.flipsort)
     for i=0,sizefs[2]-1 do $
       flipsort[initdat.flipsort[0,i]-1,initdat.flipsort[1,i]-1] = 1b
   endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Loop through spaxels
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Switch to track when first NaD normalization done
  firstnadnorm=1
; Switch to track when first continuum processed
  firstcontproc=1

  if not keyword_set(cols) then cols=[1,cube.ncols] $
  else if n_elements(cols) eq 1 then cols = [cols,cols]
  cols = fix(cols)
  for i=cols[0]-1,cols[1]-1 do begin

     if keyword_set(verbose) then $
        print,'Column ',i+1,' of ',cube.ncols,format='(A,I0,A,I0)'

     if not keyword_set(rows) then rows=[1,cube.nrows] $
     else if n_elements(rows) eq 1 then rows = [rows,rows]
     rows = fix(rows)
     for j=rows[0]-1,rows[1]-1 do begin

        novortile=0b
        if oned then begin
           flux = cube.dat[*,i]
           err = sqrt(abs(cube.var[*,i]))
           dq = cube.dq[*,i]
           labin = string(i+1,format='(I04)')
           labout = labin
        endif else begin
           if keyword_set(verbose) then $
              print,'  Row ',j+1,' of ',cube.nrows,format='(A,I0,A,I0)'
           if tag_exist(initdat,'vormap') then begin
              if finite(initdat.vormap[i,j]) AND $
                        initdat.vormap[i,j] ne bad then begin
                 iuse = vorcoords[initdat.vormap[i,j]-1,0]
                 juse = vorcoords[initdat.vormap[i,j]-1,1]
              endif else begin
                 novortile=1b
              endelse
           endif else begin
              iuse = i
              juse = j
           endelse
           if ~ novortile then begin
              flux = reform(cube.dat[iuse,juse,*],cube.nz)
              err = reform(sqrt(abs(cube.var[iuse,juse,*])),cube.nz)
              dq = reform(cube.dq[iuse,juse,*],cube.nz)
              labin = string(iuse+1,'_',juse+1,format='(I04,A,I04)')
              labout = string(i+1,'_',j+1,format='(I04,A,I04)')
           endif
        endelse

;;       Apply flux calibration
;        if tag_exist(initdat,'fluxunits') then begin
;           flux*=initdat.fluxunits
;           err*=initdat.fluxunits
;        endif
;        ;  Normalize
;        fluxmed = median(flux)
;        flux/=fluxmed
;        err/=fluxmed

;       Restore fit, after a couple of sanity checks
;       TODO: Mark zero-ed spectra as bad earlier!
        if ~ novortile then begin
           infile = initdat.outdir+initdat.label+'_'+labin+'.xdr'
           outfile = initdat.outdir+initdat.label+'_'+labout
           nodata = where(flux ne 0d,ct)
           filepresent = file_test(infile)
        endif else begin
           filepresent = 0b
           ct = 0d
        endelse
        nofit = 0b
        if ~ filepresent OR ct le 0 then begin
           nofit = 1b
           badmessage = 'No data for '+string(i+1,format='(I0)')+$
                        ', '+string(j+1,format='(I0)')+'.'
           message,badmessage,/cont
        endif else begin
         
        restore,file=infile

;       Restore original error.
;       TODO: Is this necessary?
        struct.spec_err = err[struct.fitran_indx]

        if ~ struct.noemlinfit then begin
;          Get line fit parameters
           tflux=1b
           linepars = ifsf_sepfitpars(linelist,struct.param,struct.perror,$
                                      struct.parinfo,tflux=tflux,$
                                      doublets=emldoublets)
           lineweqs = ifsf_cmpweq(struct,linelist,$
                                  doublets=emldoublets)
        endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plot emission-line data and print data to a file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        if not keyword_set(noplots) then begin

;          Plot emission lines
           if ~ struct.noemlinfit then begin
              if not linepars.nolines then begin
                 if tag_exist(initdat,'fcnpltlin') then $
                    fcnpltlin=initdat.fcnpltlin else fcnpltlin='ifsf_pltlin'
                 if tag_exist(initdat,'argspltlin1') then $
                    call_procedure,fcnpltlin,struct,initdat.argspltlin1,$
                                   outfile+'_lin1'
                 if tag_exist(initdat,'argspltlin2') then $
                    call_procedure,fcnpltlin,struct,initdat.argspltlin2,$
                                   outfile+'_lin2'
              endif
           endif

        endif 
              
;       Print fit parameters to a text file
        ifsf_printfitpar,fitlun,i+1,j+1,struct

        if ~ struct.noemlinfit then begin
;          First get correct number of components in this spaxel
           thisncomp = 0
           thisncompline = ''
           foreach line,lines_with_doublets do begin
             sigtmp = linepars.sigma[line,*]
             fluxtmp = linepars.flux[line,*]
             igd = where(sigtmp ne 0d AND sigtmp ne bad AND $
                         fluxtmp ne 0d AND fluxtmp ne bad,ctgd)
             if ctgd gt thisncomp then begin
                thisncomp = ctgd
                thisncompline = line
             endif
;            Assign total fluxes
             if ctgd gt 0 then begin
                emlweq['ftot',line,i,j]=lineweqs.tot[line]
                emlflx['ftot',line,i,j]=tflux.tflux[line]
                emlflxerr['ftot',line,i,j]=tflux.tfluxerr[line]
             endif
           endforeach
           if thisncomp eq 1 then begin
              isort = 0
              if tag_exist(initdat,'flipsort') then begin
                 if flipsort[i,j] then begin
                    message,'Flipsort set for spaxel ['+$
                            string(i+1,format='(I0)')+','+$
                            string(j+1,format='(I0)')+'] but '+$
                            'only 1 component. Setting to 2 components and '+$
                            'flipping anyway.',/cont
                    isort = [1,0]
                 endif
              endif
           endif else if thisncomp ge 2 then begin
;             Sort components
              igd = dindgen(thisncomp)
              indices = lindgen(initdat.maxncomp)
              sigtmp = linepars.sigma[thisncompline,*]
              fluxtmp = linepars.flux[thisncompline,*]
              if not tag_exist(initdat,'sorttype') then $
                 isort = sort(sigtmp[igd]) $
              else if initdat.sorttype eq 'wave' then $
                 isort = sort(linepars.wave[line,igd]) $
              else if initdat.sorttype eq 'reversewave' then $
                 isort = reverse(sort(linepars.wave[line,igd]))
              if tag_exist(initdat,'flipsort') then $
                 if flipsort[i,j] then isort = reverse(isort)
           endif
           if thisncomp gt 0 then begin
             foreach line,lines_with_doublets do begin
                kcomp = 1
                foreach sindex,isort do begin
                   cstr='c'+string(kcomp,format='(I0)')
                   emlwav[cstr,line,i,j]=linepars.wave[line,sindex]
                   emlwaverr[cstr,line,i,j]=linepars.waveerr[line,sindex]
                   emlsig[cstr,line,i,j]=linepars.sigma[line,sindex]
                   emlsigerr[cstr,line,i,j]=linepars.sigmaerr[line,sindex]
                   emlweq['f'+cstr,line,i,j]=lineweqs.comp[line,sindex]
                   emlflx['f'+cstr,line,i,j]=linepars.flux[line,sindex]
                   emlflxerr['f'+cstr,line,i,j]=linepars.fluxerr[line,sindex]
                   emlflx['f'+cstr+'pk',line,i,j]=linepars.fluxpk[line,sindex]
                   emlflxerr['f'+cstr+'pk',line,i,j]=linepars.fluxpkerr[line,sindex]
                   kcomp++
                 endforeach
              endforeach
;             Print line fluxes to a text file
              if not linepars.nolines then $ 
                 ifsf_printlinpar,lines_with_doublets,linlun,i+1,j+1,$
                                  initdat.maxncomp,linepars
           endif
        endif
        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Process and plot continuum data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;       Create / populate output data cubes                            
        if firstcontproc then begin
           hostcube = $
              {dat: dblarr(cube.ncols,cube.nrows,cube.nz),$
               err: dblarr(cube.ncols,cube.nrows,cube.nz),$
               dq:  dblarr(cube.ncols,cube.nrows,cube.nz), $
               norm_div: dblarr(cube.ncols,cube.nrows,cube.nz), $
               norm_sub: dblarr(cube.ncols,cube.nrows,cube.nz) $
              }
           if tag_exist(initdat,'decompose_ppxf_fit') then begin
              contcube = $
                 {wave: struct.wave,$
                  all_mod: dblarr(cube.ncols,cube.nrows,cube.nz),$
                  stel_mod: dblarr(cube.ncols,cube.nrows,cube.nz),$
                  poly_mod: dblarr(cube.ncols,cube.nrows,cube.nz),$
                  stel_mod_tot: dblarr(cube.ncols,cube.nrows)+bad,$
                  poly_mod_tot: dblarr(cube.ncols,cube.nrows)+bad,$
                  poly_mod_tot_pct: dblarr(cube.ncols,cube.nrows)+bad,$
                  stel_sigma: dblarr(cube.ncols,cube.nrows)+bad,$
                  stel_sigma_err: dblarr(cube.ncols,cube.nrows,2)+bad,$
                  stel_z: dblarr(cube.ncols,cube.nrows)+bad,$
                  stel_z_err: dblarr(cube.ncols,cube.nrows,2)+bad,$
                  stel_rchisq: dblarr(cube.ncols,cube.nrows)+bad,$
                  stel_ebv: dblarr(cube.ncols,cube.nrows)+bad,$
                  stel_ebv_err: dblarr(cube.ncols,cube.nrows,2)+bad $
                 }
           endif else if tag_exist(initdat,'decompose_qso_fit') then begin
              contcube = $
                 {wave: struct.wave,$
                  qso_mod: dblarr(cube.ncols,cube.nrows,cube.nz),$
                  qso_poly_mod: dblarr(cube.ncols,cube.nrows,cube.nz),$
                  host_mod: dblarr(cube.ncols,cube.nrows,cube.nz),$
                  poly_mod: dblarr(cube.ncols,cube.nrows,cube.nz),$
                  npts: dblarr(cube.ncols,cube.nrows)+bad,$
                  stel_sigma: dblarr(cube.ncols,cube.nrows)+bad,$
                  stel_sigma_err: dblarr(cube.ncols,cube.nrows,2)+bad,$
                  stel_z: dblarr(cube.ncols,cube.nrows)+bad,$
                  stel_z_err: dblarr(cube.ncols,cube.nrows,2)+bad,$
                  stel_rchisq: dblarr(cube.ncols,cube.nrows)+bad,$
                  stel_ebv: dblarr(cube.ncols,cube.nrows)+bad, $
                  stel_ebv_err: dblarr(cube.ncols,cube.nrows,2)+bad $
                 }
           endif else begin
              contcube = $
                 {all_mod: dblarr(cube.ncols,cube.nrows,cube.nz),$
                  stel_z: dblarr(cube.ncols,cube.nrows)+bad,$
                  stel_z_err: dblarr(cube.ncols,cube.nrows,2)+bad,$
                  stel_rchisq: dblarr(cube.ncols,cube.nrows)+bad,$
                  stel_ebv: dblarr(cube.ncols,cube.nrows)+bad, $
                  stel_ebv_err: dblarr(cube.ncols,cube.nrows,2)+bad $
                 }
           endelse
           firstcontproc=0
        endif
        hostcube.dat[i,j,struct.fitran_indx] = struct.cont_dat
        hostcube.err[i,j,struct.fitran_indx] = err[struct.fitran_indx]
        hostcube.dq[i,j,struct.fitran_indx] = dq[struct.fitran_indx]
        hostcube.norm_div[i,j,struct.fitran_indx] = struct.cont_dat / struct.cont_fit
        hostcube.norm_sub[i,j,struct.fitran_indx] = struct.cont_dat - struct.cont_fit
        if tag_exist(initdat,'decompose_ppxf_fit') then begin
           add_poly_degree = 4d ; this must be the same as in IFSF_FITSPEC
           if tag_exist(initdat,'argscontfit') then $
              if tag_exist(initdat.argscontfit,'add_poly_degree') then $
                 add_poly_degree = initdat.argscontfit.add_poly_degree
;          Compute polynomial
           log_rebin,[struct.wave[0],$
                      struct.wave[n_elements(struct.wave)-1]],$
                      struct.spec,$
                      dumy_log,wave_log
           xnorm = cap_range(-1d,1d,n_elements(wave_log))
           cont_fit_poly_log = 0d ; Additive polynomial
           for k=0,add_poly_degree do $
              cont_fit_poly_log += legendre(xnorm,k)*struct.ct_add_poly_weights[k]
           cont_fit_poly=interpol(cont_fit_poly_log,wave_log,ALOG(struct.wave))           
;          Compute stellar continuum
           cont_fit_stel = struct.cont_fit - cont_fit_poly
;          Total flux from different components
           cont_fit_tot = total(struct.cont_fit)
           contcube.all_mod[i,j,struct.fitran_indx] = struct.cont_fit
           contcube.stel_mod[i,j,struct.fitran_indx] = cont_fit_stel
           contcube.poly_mod[i,j,struct.fitran_indx] = cont_fit_poly
           contcube.stel_mod_tot[i,j] = total(cont_fit_stel)
           contcube.poly_mod_tot[i,j] = total(cont_fit_poly)
           contcube.poly_mod_tot_pct[i,j] = $
              contcube.poly_mod_tot[i,j] / cont_fit_tot
           contcube.stel_sigma[i,j] = struct.ct_ppxf_sigma
           contcube.stel_z[i,j] = struct.zstar
           if tag_exist(struct,'ct_errors') then $
              contcube.stel_sigma_err[i,j,*] = struct.ct_errors['ct_ppxf_sigma'] $
           else contcube.stel_sigma_err[i,j,*] = $
              [struct.ct_ppxf_sigma_err,struct.ct_ppxf_sigma_err]
           if tag_exist(struct,'ct_errors') then $
              contcube.stel_z_err[i,j,*] = struct.ct_errors['zstar'] $
           else contcube.stel_z_err[i,j,*] = $
              [struct.zstar_err,struct.zstar_err]
;;          Total flux near NaD in different components
;           ilow = value_locate(struct.wave,5850d*(1d + struct.zstar))
;           ihigh = value_locate(struct.wave,5950d*(1d + struct.zstar))
;           if ilow ne -1 AND ihigh ne -1 then begin
;              if tag_exist(initdat,'dored') then $
;                 cont_fit_nad = total(struct.cont_fit[ilow:ihigh]*$
;                                      cont_fit_redfact[ilow:ihigh]) $
;              else cont_fit_nad = total(struct.cont_fit[ilow:ihigh])
;              contcube.cont_fit_stel_nad[i,j] = total(cont_fit_stel[ilow:ihigh])
;              contcube.cont_fit_poly_nad[i,j] = total(cont_fit_poly[ilow:ihigh])
;              contcube.cont_fit_poly_nad_pct[i,j] = $
;                 contcube.cont_fit_poly_nad[i,j] / cont_fit_nad
;           endif
        endif else if tag_exist(initdat,'decompose_qso_fit') then begin
           if initdat.fcncontfit eq 'ifsf_fitqsohost' then begin
;              if tag_exist(initdat.argscontfit,'fitord') then $
;                 fitord=initdat.argscontfit.fitord else fitord=0b
              if tag_exist(initdat.argscontfit,'qsoord') then $
                 qsoord=initdat.argscontfit.qsoord else qsoord=0b
              if tag_exist(initdat.argscontfit,'hostord') then $
                 hostord=initdat.argscontfit.hostord else hostord=0b
              if tag_exist(initdat.argscontfit,'blrpar') then $
                 blrterms=n_elements(initdat.argscontfit.blrpar) else blrterms=0b
;              default here must be same as in IFSF_FITQSOHOST
              if tag_exist(initdat.argscontfit,'add_poly_degree') then $
                 add_poly_degree=initdat.argscontfit.add_poly_degree $
              else add_poly_degree=30
;             These lines mirror ones in IFSF_FITQSOHOST
              struct_tmp = struct
;             Get and renormalize template
              restore,file=initdat.argscontfit.qsoxdr
              qsowave = qsotemplate.wave
              qsoflux_full = qsotemplate.flux
              iqsoflux = where(qsowave ge struct_tmp.fitran[0]*0.99999d AND $
                               qsowave le struct_tmp.fitran[1]*1.00001d)
              qsoflux = qsoflux_full[iqsoflux]              
              qsoflux /= median(qsoflux)
              struct = struct_tmp
;             If polynomial residual is re-fit with PPXF, separate out best-fit
;             parameter structure created in IFSF_FITQSOHOST and compute polynomial
;             and stellar components
              if tag_exist(initdat.argscontfit,'refit') then begin
                 par_qsohost = struct.ct_coeff.qso_host
                 par_stel = struct.ct_coeff.stel
                 log_rebin,[struct.wave[0],$
                            struct.wave[n_elements(struct.wave)-1]],$
                            struct.spec,$
                            dumy_log,wave_log
                 xnorm = cap_range(-1d,1d,n_elements(wave_log))
                 if add_poly_degree ge 0 then begin
                   par_poly = struct.ct_coeff.poly
                   polymod_log = 0d ; Additive polynomial
                   for k=0,add_poly_degree do $
                      polymod_log += legendre(xnorm,k)*par_poly[k]
                   polymod_refit=interpol(polymod_log,wave_log,ALOG(struct.wave))
                 endif else begin
                    polymod_refit = dblarr(n_elements(struct.wave))
                 endelse
                 contcube.stel_sigma[i,j] = struct.ct_coeff.ppxf_sigma
                 contcube.stel_z[i,j] = struct.zstar
                 if tag_exist(struct,'ct_errors') then $
                    contcube.stel_sigma_err[i,j,*] = struct.ct_errors['ct_ppxf_sigma'] $
                 else contcube.stel_sigma_err[i,j,*] = $
                    [struct.ct_ppxf_sigma_err,struct.ct_ppxf_sigma_err]
                 if tag_exist(struct,'ct_errors') then $
                    contcube.stel_z_err[i,j,*] = struct.ct_errors['zstar'] $
                 else contcube.stel_z_err[i,j,*] = $
                    [struct.zstar_err,struct.zstar_err]
              endif else begin
                 par_qsohost = struct.ct_coeff
                 polymod_refit = 0d
              endelse
;             Produce fit with template only and with template + host. Also
;             output QSO multiplicative polynomial
              qsomod_polynorm=0b
              ifsf_qsohostfcn,struct.wave,par_qsohost,qsomod,qsoflux=qsoflux,$
                              /qsoonly,blrterms=blrterms,qsoscl=qsomod_polynorm,$
                              qsoord=qsoord,hostord=hostord
              hostmod = struct.cont_fit_pretweak - qsomod
;             If continuum is tweaked in any region, subide the resulting
;             residual proportionally (at each wavelength) between the QSO
;             and host components.
              qsomod_notweak = qsomod
              if tag_exist(initdat,'tweakcntfit') then begin
                 modresid = struct.cont_fit - struct.cont_fit_pretweak
                 inz = where(qsomod ne 0 AND hostmod ne 0)
                 qsofrac = dblarr(n_elements(qsomod))
                 qsofrac[inz] = qsomod[inz]/(qsomod[inz]+hostmod[inz])
                 qsomod += modresid*qsofrac
                 hostmod += modresid*(1d - qsofrac)
              endif
;             Components of QSO fit for plotting
;              qsomod_normonly = qsoflux*par_qsohost[fitord+1]
;             TODO: fix this
              qsomod_normonly = qsoflux
              if tag_exist(initdat.argscontfit,'blrpar') then $
                 ifsf_qsohostfcn,struct.wave,par_qsohost,qsomod_blronly,$
                                 qsoflux=qsoflux,/blronly,blrterms=blrterms,$
                                 qsoord=qsoord,hostord=hostord
;              qsomod_polynorm *= median(qsomod)
           endif else if initdat.fcncontfit eq 'ppxf' AND $
                         tag_exist(initdat,'qsotempfile') then begin
              struct_star = struct
              restore,file=initdat.qsotempfile
              struct_qso = struct
              struct = struct_star
              qsomod = struct_qso.cont_fit * $
                       struct.ct_coeff[n_elements(struct.ct_coeff)-1]
              hostmod = struct.cont_fit - qsomod
           endif
        endif else begin
           contcube.all_mod[i,j,struct.fitran_indx] = struct.cont_fit
           contcube.stel_z[i,j] = struct.zstar
           if tag_exist(struct,'ct_errors') then $
              contcube.stel_z_err[i,j,*] = struct.ct_errors['zstar'] $
;          for backwards compatibility
           else if tag_exist(struct,'zstar_err') then $
              contcube.stel_z_err[i,j,*] = [struct.zstar_err,struct.zstar_err] $
           else contcube.stel_z_err[i,j,*] = [0,0]
        endelse
        contcube.stel_ebv[i,j] = struct.ct_ebv
        if tag_exist(struct,'ct_errors') then $
           contcube.stel_ebv_err[i,j,*] = struct.ct_errors['ct_ebv']
;       for backwards compatibility
        if tag_exist(struct,'stel_rchisq') then $
           contcube.stel_rchisq[i,j] = struct.ct_rchisq $
        else contcube.stel_rchisq[i,j] = 0d


;       Print PPXF results to STDOUT
        if tag_exist(initdat,'decompose_ppxf_fit') OR $
           tag_exist(initdat,'decompose_qso_fit') then begin            
           if tag_exist(initdat,'argscontfit') then begin
              if tag_exist(initdat.argscontfit,'print_output') then begin
                 print,'PPXF results:'
                 if tag_exist(initdat,'decompose_ppxf_fit') then begin
                    ct_coeff_tmp = struct.ct_coeff
                    poly_tmp_pct = contcube.poly_mod_tot_pct[i,j]
                 endif else begin
                    ct_coeff_tmp = struct.ct_coeff.stel
                    poly_tmp_pct = total(polymod_refit) / total(hostmod)
                 endelse
                 inz = where(ct_coeff_tmp ne 0d,ctnz)
                 if ctnz gt 0 then begin
                    coeffgd = ct_coeff_tmp[inz]
;                   normalize coefficients to % of total stellar coeffs.
                    totcoeffgd = total(coeffgd)
                    coeffgd /= totcoeffgd
;                   re-normalize to % of total flux
                    coeffgd *= (1d - poly_tmp_pct)
                    restore,initdat.startempfile
                    agesgd = template.ages[inz]
;                   sum coefficients over age ranges
                    iyoung = where(agesgd le 1d7,ctyoung)
                    iinter1 = where(agesgd gt 1d7 AND agesgd le 1d8,ctinter1)
                    iinter2 = where(agesgd gt 1d8 AND agesgd le 1d9,ctinter2)
                    iold = where(agesgd gt 1d9,ctold)
                    if ctyoung gt 0 then $
                       coeffyoung = total(coeffgd[iyoung])*100d $
                    else coeffyoung=0d
                    if ctinter1 gt 0 then $
                       coeffinter1 = total(coeffgd[iinter1])*100d $
                    else coeffinter1=0d
                    if ctinter2 gt 0 then $
                       coeffinter2 = total(coeffgd[iinter2])*100d $
                    else coeffinter2=0d
                    if ctold gt 0 then $
                       coeffold = total(coeffgd[iold])*100d $
                    else coeffold=0d
                    print,'   ',string(round(coeffyoung),format='(I0)')+$
                       '% contribution from ages <= 10 Myr.'
                    print,'   ',string(round(coeffinter1),format='(I0)')+$
                       '% contribution from 10 Myr < age <= 100 Myr.'
                    print,'   ',string(round(coeffinter2),format='(I0)')+$
                       '% contribution from 100 Myr < age <= 1 Gyr.'
                    print,'   ',string(round(coeffold),format='(I0)')+$
                       '% contribution from ages > 1 Gyr.'
                 endif
                 print,'   ','Stellar template convolved with sigma = '+$
                       string(struct.ct_ppxf_sigma,format='(I0)')+' km/s'
              endif
           endif
        endif
        
;       Plot QSO- and host-only continuum fit
        if tag_exist(initdat,'decompose_qso_fit') then begin
           struct_host = struct
           struct_host.spec -= qsomod
           struct_host.cont_dat -= qsomod
           struct_host.cont_fit -= qsomod
           struct_qso = struct
           struct_qso.spec -= hostmod
           struct_qso.cont_dat -= hostmod
           struct_qso.cont_fit -= hostmod
           contcube.qso_mod[i,j,struct.fitran_indx] = qsomod
           contcube.qso_poly_mod[i,j,struct.fitran_indx] = qsomod_polynorm
           contcube.host_mod[i,j,struct.fitran_indx] = hostmod
           contcube.poly_mod[i,j,struct.fitran_indx] = polymod_refit
           contcube.npts[i,j] = n_elements(struct.fitran_indx)
           if tag_exist(initdat,'remove_scattered') then $
              contcube.host_mod[i,j,struct.fitran_indx] -= polymod_refit
;          Data minus (emission line model + QSO model)
;           contcube.host_dat[i,j,*] = struct.cont_dat - qsomod
;          Update hostcube.dat to remove tweakcnt mods
;          Data minus (emission line model + QSO model, tweakcnt mods not 
;          included in QSO model)
           hostcube.dat[i,j,struct.fitran_indx] = struct.cont_dat - qsomod_notweak
           if ~keyword_set(noplots) AND $
              total(struct_host.cont_fit) ne 0d then begin
              if tag_exist(initdat.argscontfit,'refit') then begin
                 compspec = [[polymod_refit],[hostmod-polymod_refit]]
                 comptit = ['ord. '+string(add_poly_degree,format='(I0)')+$
                            ' Leg. poly.','stel. temp.']
              endif else begin
                 compspec = hostmod
                 comptit = ['exponential terms']
              endelse
              if tag_exist(initdat,'argspltcont') then $
                 call_procedure,fcnpltcont,struct_host,outfile+'_cnt_host',$
                                compspec=compspec,comptit=comptit,title='Host',$
                                fitran=initdat.fitran,_extra=initdat.argspltcont $
              else $
                 call_procedure,fcnpltcont,struct_host,outfile+'_cnt_host',$
                                compspec=compspec,comptit=comptit,title='Host',$
                                fitran=initdat.fitran
              if tag_exist(initdat.argscontfit,'blrpar') then begin
                 qsomod_blrnorm = $
                    median(qsomod)/max(qsomod_blronly)
                 compspec = [[qsomod_normonly],[qsomod_blronly*qsomod_blrnorm]]
                 comptit = ['raw template','scattered$\times$'+$
                            string(qsomod_blrnorm,format='(D0.2)')]
              endif else begin
                 compspec = [[qsomod_normonly]]
                 comptit = ['raw template']
              endelse
              if tag_exist(initdat,'argspltcont') then $
                 call_procedure,fcnpltcont,struct_qso,outfile+'_cnt_qso',$
                                compspec=compspec,comptit=comptit,title='QSO',$
                                fitran=initdat.fitran,_extra=initdat.argspltcont $
              else $
                 call_procedure,fcnpltcont,struct_qso,outfile+'_cnt_qso',$
                                compspec=compspec,comptit=comptit,title='QSO',$
                                fitran=initdat.fitran
           endif
        endif

;       Plot continuum
;       Make sure fit doesn't indicate no continuum; avoids
;       plot range error in continuum fitting routine, as well as a blank
;       plot!
        if ~keyword_set(noplots) AND total(struct.cont_fit) ne 0d then $
           if tag_exist(initdat,'decompose_qso_fit') then $
              if tag_exist(initdat,'argspltcont') then $
                 call_procedure,fcnpltcont,struct,outfile+'_cnt',$
                                compspec=[[qsomod],[hostmod]],$
                                title='Total',comptit=['QSO','host'],$
                                fitran=initdat.fitran,_extra=initdat.argspltcont $
               else $
                 call_procedure,fcnpltcont,struct,outfile+'_cnt',$
                                compspec=[[qsomod],[hostmod]],$
                                title='Total',comptit=['QSO','host'],$
                                fitran=initdat.fitran $
           else if tag_exist(initdat,'decompose_ppxf_fit') then $
              if tag_exist(initdat,'argspltcont') then $
                 call_procedure,fcnpltcont,struct,outfile+'_cnt',$
                                compspec=[[cont_fit_stel],[cont_fit_poly]],$
                                title='Total',$
                                comptit=['stel. temp.','ord. '+$
                                         string(add_poly_degree,format='(I0)')+$
                                         ' Leg. poly'],$
                                fitran=initdat.fitran,_extra=initdat.argspltcont $
              else $
                 call_procedure,fcnpltcont,struct,outfile+'_cnt',$
                                compspec=[[cont_fit_stel],[cont_fit_poly]],$
                                title='Total',$
                                comptit=['stel. temp.','ord. '+$
                                         string(add_poly_degree,format='(I0)')+$
                                         ' Leg. poly'],$
                                fitran=initdat.fitran $
           else $
              if tag_exist(initdat,'argspltcont') then $
                 call_procedure,fcnpltcont,struct,outfile+'_cnt',$
                                fitran=initdat.fitran,_extra=initdat.argspltcont $
              else $
                 call_procedure,fcnpltcont,struct,outfile+'_cnt',$
                                fitran=initdat.fitran
                                
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Process NaD (normalize, compute quantities and save, plot)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        if tag_exist(initdat,'donad') then begin

           if tag_exist(initdat,'decompose_qso_fit') then begin
              if tag_exist(initdat,'remove_scattered') then begin
                 hostmod_tmp = hostmod - polymod_refit
                 qsomod_tmp = qsomod + polymod_refit
              endif else begin
                 hostmod_tmp = hostmod
                 qsomod_tmp = qsomod
              endelse
              nadnormcont = (struct.cont_dat - qsomod_tmp)/hostmod_tmp
              nadnormconterr = struct.spec_err/hostmod_tmp
              nadnormstel = hostmod_tmp
           endif else begin
              nadnormcont = struct.cont_dat/struct.cont_fit
              nadnormconterr = struct.spec_err/struct.cont_fit
              nadnormstel = struct.cont_fit
           endelse
         
           if tag_exist(initnad,'argsnormnad') then begin
              normnad = ifsf_normnad(struct.wave,$
                                     nadnormcont,$
                                     nadnormconterr,$
                                     struct.zstar,fitpars_normnad,$
                                     _extra=initnad.argsnormnad)
              normnadem = ifsf_normnad(struct.wave,$
                                       struct.emlin_dat,$
                                       struct.spec_err,$
                                       struct.zstar,fitpars_normnadem,$
                                       /nosncut,/subtract,$
                                       _extra=initnad.argsnormnad)
              normnadstel = ifsf_normnad(struct.wave,$
                                         nadnormstel,$
                                         struct.spec_err,$
                                         struct.zstar,fitpars_normnadstel,$
                                         _extra=initnad.argsnormnad)
              if ~ tag_exist(initnad.argsnormnad,'fitranlo') then $
                 fitranlo = (1d +struct.zstar)*[5810d,5865d] $
              else fitranlo = initnad.argsnormnad.fitranlo
              if ~ tag_exist(initnad.argsnormnad,'fitranhi') then $
                 fitranhi = (1d +struct.zstar)*[5905d,5960d] $
              else fitranhi = initnad.argsnormnad.fitranhi
           endif else begin
              normnad = ifsf_normnad(struct.wave,$
                                     nadnormcont,$
                                     nadnormconterr,$
                                     struct.zstar,fitpars_normnad)
              normnadem = ifsf_normnad(struct.wave,$
                                       struct.emlin_dat,$
                                       struct.spec_err,$
                                       struct.zstar,fitpars_normnadem,$
                                       /nosncut,/subtract)
              normnadstel = ifsf_normnad(struct.wave,$
                                         nadnormstel,$
                                         struct.spec_err,$
                                         struct.zstar,fitpars_normnadstel)
              fitranlo = (1d +struct.zstar)*[5810d,5865d]
              fitranhi = (1d +struct.zstar)*[5905d,5960d]

           endelse
;          Check data quality
           if normnad ne !NULL then igd = where(normnad.nflux gt 0d,ctgd) $
           else ctgd = 0
;          Compute empirical equivalent widths and emission-line fluxes
           if ctgd gt 0 then begin
;          Create output data cube
              if firstnadnorm then begin
                 ilo = value_locate(cube.wave,fitranlo[0])+1
                 ihi = value_locate(cube.wave,fitranhi[1])
                 dat_normnadwave = cube.wave[ilo:ihi]
                 nz = ihi-ilo+1
                 nadcube = $
                    {wave: dblarr(cube.ncols,cube.nrows,nz),$
                     cont: dblarr(cube.ncols,cube.nrows,nz),$
                     dat: dblarr(cube.ncols,cube.nrows,nz),$
                     err: dblarr(cube.ncols,cube.nrows,nz)+bad,$
                     weq: dblarr(cube.ncols,cube.nrows,4)+bad,$
                     stelweq: dblarr(cube.ncols,cube.nrows,2)+bad,$
                     iweq: dblarr(cube.ncols,cube.nrows,4)+bad,$
                     emflux: dblarr(cube.ncols,cube.nrows,2)+bad,$
                     emul: dblarr(cube.ncols,cube.nrows,4)+bad,$
                     vel: dblarr(cube.ncols,cube.nrows,6)+bad}
                 firstnadnorm = 0
              endif
;             Defaults
              emflux=dblarr(2)
              emul=dblarr(4)+bad
              vel = dblarr(6)+bad
              if tag_exist(initnad,'argsnadweq') then begin
                 weq = ifsf_cmpnadweq(normnad.wave,normnad.nflux,normnad.nerr,$
                                      snflux=normnadem.nflux,unerr=normnadem.nerr,$
                                      emflux=emflux,emul=emul,$
                                      _extra=initnad.argsnadweq)

;                These need to be compatible with the IFSF_CMPNADWEQ defaults
                 if tag_exist(initnad.argsnadweq,'emwid') then $
                    emwid=initnad.argsnadweq.emwid else emwid=20d
                 if tag_exist(initnad.argsnadweq,'iabsoff') then $
                    iabsoff=initnad.argsnadweq.iabsoff else iabsoff=4l
              endif else begin
                 weq = ifsf_cmpnadweq(normnad.wave,normnad.nflux,normnad.nerr,$
                                      snflux=normnadem.nflux,unerr=normnadem.nerr,$
                                      emflux=emflux,emul=emul)
;                These need to be compatible with the IFSF_CMPNADWEQ defaults
                 emwid=20d
                 iabsoff=4l
              endelse
;             Compute stellar continuum NaD equivalent widths from fit
              stelweq = ifsf_cmpnadweq(normnadstel.wave,normnadstel.nflux,$
                                       normnadstel.nerr,$
                                       wavelim=[5883d*(1d +initdat.zsys_gas),$
                                                6003d*(1d +initdat.zsys_gas),$
                                                0d,0d])

;             Compute empirical velocities
              size_weq = size(weq)
              if size_weq[0] eq 2 then begin
                 if tag_exist(initnad,'argsnadvel') then $
                    vel = ifsf_cmpnadvel(normnad.wave,normnad.nflux,normnad.nerr,$
                                         weq[*,1],initdat.zsys_gas,$
                                         _extra=initnad.argsnadvel) $
                 else vel = ifsf_cmpnadvel(normnad.wave,normnad.nflux,normnad.nerr,$
                                           weq[*,1],initdat.zsys_gas)
              endif
              
              ilo = where(normnad.wave[0] eq dat_normnadwave)
              ihi = where(normnad.wave[n_elements(normnad.wave)-1] $
                    eq dat_normnadwave)

;             Assume that stellar fit is a good model but that the error spectrum
;             may not be perfect. Correct using stellar reduced chi squared
              if tag_exist(initnad,'errcorr_ctrchisq') then begin
                 normnad.nerr *= struct.ct_rchisq
                 weq[1,0] *= struct.ct_rchisq
                 weq[3,0] *= struct.ct_rchisq
              endif
              nadcube.wave[i,j,*] = dat_normnadwave
              nadcube.cont[i,j,ilo:ihi] = struct.cont_fit[normnad.ind]
              nadcube.dat[i,j,ilo:ihi] = normnad.nflux
              nadcube.err[i,j,ilo:ihi] = normnad.nerr
              nadcube.weq[i,j,*] = weq[*,0]
              nadcube.iweq[i,j,*] = weq[*,1]
              nadcube.stelweq[i,j,*] = stelweq[0:1]
              nadcube.emflux[i,j,*] = emflux
              nadcube.emul[i,j,*] = emul
              nadcube.vel[i,j,*] = vel
;             Plot data
              if not keyword_set(noplots) then $
                 if tag_exist(initnad,'argspltnormnad') then $
                    ifsf_pltnaddat,normnad,fitpars_normnad,struct.zstar,$
                                   outfile+'_nad_norm',autoindices=weq[*,1],$
                                   emwid=emwid,iabsoff=iabsoff,$
                                   _extra=initnad.argspltnormnad else $
                    ifsf_pltnaddat,normnad,fitpars_normnad,struct.zstar,$
                                   outfile+'_nad_norm',autoindices=weq[*,1],$
                                   emwid=emwid,iabsoff=iabsoff
           endif

        endif

        endelse

     endfor

  endfor

  free_lun,fitlun

  if ~ tag_exist(initdat,'noemlinfit') then begin
     free_lun,linlun
     emlkeys=!NULL
;    Apply a sigma cut on a line-by-line basis
;    Apply only to total flux summed over all components.
;    Previous sigma cuts may have detected a line but not significantly in 
;    every species.
;    Presently this means that sums of doublet fluxes may have different cuts
;    than the components of the doublet
     if tag_exist(initdat,'emlsigcut') then begin
        foreach line,lines_with_doublets do begin
           ibd = where(emlflx['ftot',line] gt 0d AND $
                       emlflx['ftot',line] lt $
                       emlflxerr['ftot',line]*initdat.emlsigcut,ctbd)
           if ctbd gt 0 then begin
              emlweq['ftot',line,ibd] = bad
              emlflx['ftot',line,ibd] = bad
              if ~ tag_exist(initdat,'emlkeeperr') then $
                 emlflxerr['ftot',line,ibd] = bad
              for k=0,initdat.maxncomp-1 do begin
                 cstr='c'+string(k+1,format='(I0)')
                 emlweq['f'+cstr,line,ibd] = bad
                 emlflx['f'+cstr,line,ibd] = bad
                 emlflx['f'+cstr+'pk',line,ibd] = bad
                 if ~ tag_exist(initdat,'emlkeeperr') then begin
                    emlflxerr['f'+cstr,line,ibd] = bad
                    emlflxerr['f'+cstr+'pk',line,ibd] = bad
                 endif
                 emlwav[cstr,line,ibd] = bad
                 emlwaverr[cstr,line,ibd] = bad
                 emlsig[cstr,line,ibd] = bad
                 emlsigerr[cstr,line,ibd] = bad
              endfor
           endif
        endforeach
     endif
;    Compute CVDFs
     if not tag_exist(initdat,'nocvdf') then begin
        if tag_exist(initdat,'cvdf_vlimits') then vlimits=initdat.cvdf_vlimits $
        else vlimits=0
        if tag_exist(initdat,'cvdf_vstep') then vstep=initdat.cvdf_vstep $
        else vstep=0
        emlcvdf = ifsf_cmpcvdf(emlwav,emlwaverr,emlsig,emlsigerr,emlflx,emlflxerr,$
                               initdat.maxncomp,linelist_with_doublets,$
                               initdat.zsys_gas,vlimits=vlimits,vstep=vstep)
        save,emlwav,emlwaverr,emlsig,emlsigerr,emlweq,emlflx,emlflxerr,emlcvdf,$
             file=initdat.outdir+initdat.label+'.lin.xdr'
     endif else begin
        save,emlwav,emlwaverr,emlsig,emlsigerr,emlweq,emlflx,emlflxerr,$
             file=initdat.outdir+initdat.label+'.lin.xdr'
     endelse
  endif

  ;  if tag_exist(initdat,'decompose_ppxf_fit') OR $
  ;     tag_exist(initdat,'decompose_qso_fit') then $
  save,contcube,file=initdat.outdir+initdat.label+'.cont.xdr'

   if tag_exist(initdat,'host') then begin
;     Change initial wavelength -- use last restored fit structure. For some 
;     reason SXADDPAR doesn't work when applied directly to structures, so
;     create new variables.
      newheader_dat = header.dat
      newheader_var = header.var
      newheader_dq = header.dq
;     These lines are for producing a host spectrum in one dimension.
      if tag_exist(initdat.host,'singlespec') then begin
         sxaddpar,newheader_dat,'CRVAL1',cube.wave[0]
         sxaddpar,newheader_var,'CRVAL1',cube.wave[0]
         sxaddpar,newheader_dq,'CRVAL1',cube.wave[0]
         sxaddpar,newheader_dat,'CRPIX1',1
         sxaddpar,newheader_var,'CRPIX1',1
         sxaddpar,newheader_dq,'CRPIX1',1
         cdelt1 = sxpar(header.dat,'CDELT1')
         sxaddpar,newheader_dat,'CDELT1',cdelt1
         sxaddpar,newheader_var,'CDELT1',cdelt1
         sxaddpar,newheader_dq,'CDELT1',cdelt1
         writefits,initdat.host.dat_fits,cube.phu,header.phu
         writefits,initdat.host.dat_fits,$
                   reform(hostcube.dat,n_elements(hostcube.dat)),$
                   newheader_dat,/append
         if dqext eq 2 AND varext eq 3 then begin
            writefits,initdat.host.dat_fits,$
                      reform(hostcube.dq,n_elements(hostcube.dat)),$
                      newheader_dq,/append
            writefits,initdat.host.dat_fits,$
                      reform(hostcube.err,n_elements(hostcube.dat))^2d,$
                      newheader_var,/append
         endif else begin
            writefits,initdat.host.dat_fits,$
                      reform(hostcube.err,n_elements(hostcube.dat))^2d,$
                      newheader_var,/append
            writefits,initdat.host.dat_fits,$
                      reform(hostcube.dq,n_elements(hostcube.dat)),$
                      newheader_dq,/append
         endelse
      endif else begin
;     These lines are for producing a host spectrum in the third dimension.
;        Check to see if wave extension needs to be added. Note that it is added
;        as the last extension, regardless of what wave extension is specified
;        in the input file.
         writewaveext=0b
         if tag_exist(initdat,'argsreadcube') then $
            if tag_exist(initdat.argsreadcube,'waveext') then $
               writewaveext=1b
         sxaddpar,newheader_dat,'CRVAL3',cube.wave[0]
         sxaddpar,newheader_var,'CRVAL3',cube.wave[0]
         sxaddpar,newheader_dq,'CRVAL3',cube.wave[0]
         sxaddpar,newheader_dat,'CRPIX3',1
         sxaddpar,newheader_var,'CRPIX3',1
         sxaddpar,newheader_dq,'CRPIX3',1
;        case of no PHU
         if datext lt 0 then begin
            writefits,initdat.host.dat_fits,hostcube.dat,newheader_dat
;        case of PHU
         endif else begin
            writefits,initdat.host.dat_fits,cube.phu,header.phu
            writefits,initdat.host.dat_fits,hostcube.dat,newheader_dat,/append
         endelse
         if dqext eq 2 AND varext eq 3 then begin
            writefits,initdat.host.dat_fits,hostcube.dq,newheader_dq,/append
            writefits,initdat.host.dat_fits,hostcube.err^2d,newheader_var,/append
         endif else begin
            writefits,initdat.host.dat_fits,hostcube.err^2d,newheader_var,/append
            writefits,initdat.host.dat_fits,hostcube.dq,newheader_dq,/append
         endelse
         if writewaveext then $
            writefits,initdat.host.dat_fits,cube.wave,header.wave,/append
         if tag_exist(initdat.host,'norm_div') then begin
            if datext lt 0 then begin
               writefits,initdat.host.norm_div,hostcube.norm_div,newheader_dat
            endif else begin
               writefits,initdat.host.norm_div,cube.phu,header.phu
               writefits,initdat.host.norm_div,hostcube.norm_div,newheader_dat,/append
            endelse
            if writewaveext then $
               writefits,initdat.host.norm_div,cube.wave,header.wave,/append
         endif
         if tag_exist(initdat.host,'norm_sub') then begin
            if datext lt 0 then begin
               writefits,initdat.host.norm_sub,hostcube.norm_sub,newheader_dat
            endif else begin
               writefits,initdat.host.norm_sub,cube.phu,header.phu
               writefits,initdat.host.norm_sub,hostcube.norm_sub,newheader_dat,/append
            endelse
            if writewaveext then $
               writefits,initdat.host.norm_sub,cube.wave,header.wave,/append
         endif
      endelse
   endif

  if tag_exist(initdat,'donad') then $
     save,nadcube,file=initdat.outdir+initdat.label+'.nadspec.xdr'

end
