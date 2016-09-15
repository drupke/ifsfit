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
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Various plots and data files.
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
;
; :Copyright:
;    Copyright (C) 2013--2016 David S. N. Rupke
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
          verbose=verbose,_extra=ex

  bad = 1d99
  fwhmtosig = 2d*sqrt(2d*alog(2d))
  if keyword_set(verbose) then quiet=0 else quiet=1
  if keyword_set(oned) then oned=1 else oned=0

; Get fit initialization
  if keyword_set(_extra) then $
     initdat = call_function(initproc,_extra=ex) $
  else $
     initdat = call_function(initproc)   
  if tag_exist(initdat,'donad') then begin
    initnad={dumy: 1}
    if keyword_set(_extra) then $
      initdat = call_function(initproc,initnad=initnad,_extra=ex) $
    else $
      initdat = call_function(initproc,initnad=initnad)
  endif

  if ~ tag_exist(initdat,'noemlinfit') then begin
;   Get linelist
    linelist = ifsf_linelist(initdat.lines)
    nlines = linelist.count()
  endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if not tag_exist(initdat,'datext') then datext=1 else datext=initdat.datext
  if not tag_exist(initdat,'varext') then varext=2 else varext=initdat.varext
  if not tag_exist(initdat,'dqext') then dqext=3 else dqext=initdat.dqext

  cube = ifsf_readcube(initdat.infile,quiet=quiet,oned=oned,$
                       datext=datext,varext=varext,dqext=dqext)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Initialize output files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if ~ tag_exist(initdat,'noemlinfit') then begin
     if not tag_exist(initdat,'outlines') then $
        outlines = (linelist->keys())->toarray() $
     else outlines = initdat.outlines
     ifsf_printlinpar,outlines,linlun,$
                      outfile=initdat.outdir+initdat.label+'.lin.dat'
  endif
  ifsf_printfitpar,fitlun,$
                   outfile=initdat.outdir+initdat.label+'.fit.dat'
                   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Initialize line hash
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if ~ tag_exist(initdat,'noemlinfit') then begin
      linmaps = hash()
      tlinmaps = hash()
      foreach line,outlines do begin
         linmaps[line] = dblarr(cube.ncols,cube.nrows,initdat.maxncomp,5) + bad
         tlinmaps[line] = dblarr(cube.ncols,cube.nrows,1,2) + bad
      endforeach
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

        if oned then begin
           flux = cube.dat[*,i]
           err = sqrt(abs(cube.var[*,i]))
           lab = string(i+1,format='(I04)')
        endif else begin
           if keyword_set(verbose) then $
              print,'  Row ',j+1,' of ',cube.nrows,format='(A,I0,A,I0)'
           flux = reform(cube.dat[i,j,*],cube.nz)
           err = reform(sqrt(abs(cube.var[i,j,*])),cube.nz)
           lab = string(i+1,'_',j+1,format='(I04,A,I04)')
        endelse

;       Restore fit, after a couple of sanity checks
;       TODO: Mark zero-ed spectra as bad earlier!
        infile = initdat.outdir+initdat.label+'_'+lab+'.xdr'
        outfile = initdat.outdir+initdat.label+'_'+lab
        if ~ file_test(infile) then begin
           print,'IFSFA: No XDR file for ',i+1,', ',j+1,'.',$
           format='(A0,I4,A0,I4,A0)'
           goto,nofit
        endif
        nodata = where(flux ne 0d,ct)
        if ct le 0 then begin
           print,'IFSFA: Data is everywhere 0 for ',i+1,', ',j+1,'.',$
           format='(A0,I4,A0,I4,A0)'
           goto,nofit
        endif
        restore,file=infile

;       Restore original error.
;       TODO: Is this necessary?
        struct.spec_err = err[struct.gd_indx]

        if ~ struct.noemlinfit then begin
;          Get line fit parameters
           tflux=1b
           linepars = ifsf_sepfitpars(linelist,struct.param,struct.perror,$
                                      struct.parinfo,tflux=tflux)
        endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plot emission-line data and print data to a file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        if not keyword_set(noplots) then begin

;          Plot continuum
           if tag_exist(initdat,'fcnpltcont') then $
              fcnpltcont=initdat.fcnpltcont $
           else fcnpltcont='ifsf_pltcont'
;          Make sure fit doesn't indicate no continuum; avoids
;          plot range error in continuum fitting routine, as well as a blank
;          plot!
           if total(struct.cont_fit) ne 0d then $
              call_procedure,fcnpltcont,struct,outfile+'_cnt'
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
;          Save fit parameters to an array
           foreach line,outlines do begin
              linmaps[line,i,j,*,*] = [[linepars.flux[line,*]],$
                                       [linepars.fluxerr[line,*]],$
                                       [linepars.wave[line,*]],$
                                       [linepars.sigma[line,*]],$
                                       [linepars.fluxpk[line,*]]]
              tlinmaps[line,i,j,0,*] = [[tflux.tflux[line]],$
                                        [tflux.tfluxerr[line]]]
           endforeach

;          Print line fluxes and Halpha Weq to a text file
           if not linepars.nolines then $ 
              ifsf_printlinpar,outlines,linlun,i+1,j+1,initdat.maxncomp,linepars

        endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Process continuum data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;       Create / populate output data cube                            
        if firstcontproc then begin
           if tag_exist(initdat,'decompose_ppxf_fit') then begin
              contcube = $
                 {cont_fit_stel_tot: dblarr(cube.ncols,cube.nrows)+bad,$
                  cont_fit_poly_tot: dblarr(cube.ncols,cube.nrows)+bad,$
                  cont_fit_poly_tot_pct: dblarr(cube.ncols,cube.nrows)+bad,$
                  cont_fit_stel_nad: dblarr(cube.ncols,cube.nrows)+bad,$
                  cont_fit_poly_nad: dblarr(cube.ncols,cube.nrows)+bad,$
                  cont_fit_poly_nad_pct: dblarr(cube.ncols,cube.nrows)+bad,$
                  cont_fit_ebv: dblarr(cube.ncols,cube.nrows)+bad,$
                  cont_fit_ages: dblarr(cube.ncols,cube.nrows,3)+bad,$
                  cont_fit_norms: dblarr(cube.ncols,cube.nrows,3)+bad}
           endif else if tag_exist(initdat,'decompose_qso_fit') then begin
              contcube = $
                 {qso: dblarr(cube.ncols,cube.nrows,$
                              n_elements(struct.gd_indx)),$
                  host: dblarr(cube.ncols,cube.nrows,$
                               n_elements(struct.gd_indx))}
           endif
           firstcontproc=0
        endif
        if tag_exist(initdat,'decompose_ppxf_fit') then begin
;          Get stellar templates
           restore,initdat.startempfile
;          Redshift stellar templates
           templatelambdaz = $
              reform(template.lambda,n_elements(template.lambda)) * $
              (1d + struct.zstar)
           if tag_exist(initdat,'vacuum') then airtovac,templatelambdaz
;          Interpolate template to same grid as data
           temp = ifsf_interptemp(struct.wave,templatelambdaz,template.flux)
;          Compute stellar continuum
           cont_fit_stel = temp # struct.ct_coeff                                
           if tag_exist(initdat,'dored') then begin
              contcube.cont_fit_ebv[i,j] = struct.ct_ebv
              cont_fit_redfact = ppxf_reddening_curve(struct.wave,struct.ct_ebv)
              cont_fit_stel *= cont_fit_redfact
           endif
           cont_fit_poly = struct.cont_fit - cont_fit_stel
;          Total flux from different components
           if tag_exist(initdat,'dored') then $
              cont_fit_tot = total(struct.cont_fit * cont_fit_redfact) $
           else cont_fit_tot = total(struct.cont_fit)
           contcube.cont_fit_stel_tot[i,j] = total(cont_fit_stel)
           contcube.cont_fit_poly_tot[i,j] = total(cont_fit_poly)
           contcube.cont_fit_poly_tot_pct[i,j] = $
              contcube.cont_fit_poly_tot[i,j] / cont_fit_tot
;          Total flux near NaD in different components
           ilow = value_locate(struct.wave,5850d*(1d + struct.zstar))
           ihigh = value_locate(struct.wave,5950d*(1d + struct.zstar))
           if ilow ne -1 AND ihigh ne -1 then begin
              if tag_exist(initdat,'dored') then $
                 cont_fit_nad = total(struct.cont_fit[ilow:ihigh]*$
                                      cont_fit_redfact[ilow:ihigh]) $
              else cont_fit_nad = total(struct.cont_fit[ilow:ihigh])
              contcube.cont_fit_stel_nad[i,j] = total(cont_fit_stel[ilow:ihigh])
              contcube.cont_fit_poly_nad[i,j] = total(cont_fit_poly[ilow:ihigh])
              contcube.cont_fit_poly_nad_pct[i,j] = $
                 contcube.cont_fit_poly_nad[i,j] / cont_fit_nad
           endif
;          Find best fit stellar ages and coefficients. Include only three 
;          templates with largest coefficients, ordered in decreasing importance.
           icoeffgd = where(struct.ct_coeff ne 0,countgd)
           if countgd gt 0 then begin
              coeffgd = struct.ct_coeff[icoeffgd]
              totcoeffgd = total(coeffgd)
              coeffgd /= totcoeffgd
              agesgd = template.ages[icoeffgd]
              sortcoeffgd = reverse(sort(coeffgd))
              if countgd gt 3 then countgd=3
              contcube.cont_fit_norms[i,j,0:countgd-1] = $
                 coeffgd[sortcoeffgd[0:countgd-1]]
              contcube.cont_fit_ages[i,j,0:countgd-1] = $
                 agesgd[sortcoeffgd[0:countgd-1]]
           endif
        endif
        if tag_exist(initdat,'decompose_qso_fit') then begin
           if initdat.fcncontfit eq 'ifsf_fitqsohost' then begin
              if tag_exist(initdat.argscontfit,'fitord') then $
                 fitord=initdat.argscontfit.fitord else fitord=0b
              if tag_exist(initdat.argscontfit,'qsoord') then $
                 qsoord=initdat.argscontfit.qsoord else qsoord=0b
              if tag_exist(initdat.argscontfit,'expterms') then $
                 expterms=initdat.argscontfit.expterms else expterms=0b
;             These lines mirror ones in IFSF_FITQSOHOST
              struct_tmp = struct
;             Get and renormalize template
              restore,file=initdat.argscontfit.qsoxdr
              qsowave = struct.wave
              qsoflux = struct.cont_fit
              qsoflux /= median(qsoflux)
              struct = struct_tmp
;             Produce fit with template only and with template + host
              ifsf_qsohostfcn,struct.wave,struct.ct_coeff,qsomod,fitord=fitord,$
                              qsoord=qsoord,expterms=expterms,/qsoonly,$
                              qsoflux=qsoflux
              ifsf_qsohostfcn,struct.wave,struct.ct_coeff,hostmod,fitord=fitord,$
                              qsoord=qsoord,expterms=expterms,/hostonly,$
                              qsoflux=qsoflux
;             If continuum is tweaked in any region, divide the resulting
;             residual proportionally (at each wavelength) between the QSO
;             and host components.
              if tag_exist(initdat,'tweakcntfit') then begin
                 modresid = struct.cont_fit - (qsomod+hostmod)
                 qsofrac = qsomod/(qsomod+hostmod)
                 qsomod += modresid*qsofrac
                 hostmod += modresid*(1d - qsofrac)
              endif
              contcube.qso[i,j,*] = qsomod
;              hostmod = totmod - qsomod
;              contcube.host[i,j,*] = hostmod
              contcube.host[i,j,*] = hostmod
           endif else if initdat.fcncontfit eq 'ppxf' AND $
                         tag_exist(initdat,'qsotempfile') then begin
              struct_star = struct
              restore,file=initdat.qsotempfile
              struct_qso = struct
              struct = struct_star
              qsomod = struct_qso.cont_fit * $
                       struct.ct_coeff[n_elements(struct.ct_coeff)-1]
              contcube.qso[i,j,*] = qsomod
              hostmod = struct.cont_fit - qsomod
              contcube.host[i,j,*] = hostmod
           endif
        endif
        
;       Plot host-only continuum fit
        if tag_exist(initdat,'decompose_qso_fit') then begin
           struct_host = struct
           struct_host.spec -= qsomod
           struct_host.cont_dat -= qsomod
           struct_host.cont_fit -= qsomod
;          Make sure fit to host doesn't indicate no continuum; avoids
;          plot range error in continuum fitting routine, as well as a blank
;          plot!
           if total(struct_host.cont_fit) ne 0d then $
              call_procedure,fcnpltcont,struct_host,outfile+'_cnt_host'
        endif

        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Process NaD (normalize, compute quantities and save, plot)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        if tag_exist(initdat,'donad') then begin

           if tag_exist(initdat,'decompose_qso_fit') then begin
              nadnormcont = (struct.cont_dat - qsomod)/hostmod
              nadnormconterr = struct.spec_err/hostmod
              nadnormstel = hostmod
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
                                       struct.zstar,fitpars_normnadem,/subtract,$
                                       _extra=initnad.argsnormnad)
              normnadstel = ifsf_normnad(struct.wave,$
                                         nadnormstel,$
                                         struct.spec_err,$
                                         struct.zstar,fitpars_normnadstel,$
                                         _extra=initnad.argsnormnad)
           endif else begin
              normnad = ifsf_normnad(struct.wave,$
                                     nadnormcont,$
                                     nadnormconterr,$
                                     struct.zstar,fitpars_normnad)
              normnadem = ifsf_normnad(struct.wave,$
                                       struct.emlin_dat,$
                                       struct.spec_err,$
                                       struct.zstar,fitpars_normnadem,/subtract)
              normnadstel = ifsf_normnad(struct.wave,$
                                         nadnormstel,$
                                         struct.spec_err,$
                                         struct.zstar,fitpars_normnadstel)
           endelse
;          Compute empirical equivalent widths and emission-line fluxes
           emflux=dblarr(2)
           emul=dblarr(4)+bad
           if tag_exist(initnad,'argsnadweq') then begin
              weq = ifsf_cmpnadweq(normnad.wave,normnad.nflux,normnad.nerr,$
                                      snflux=normnadem.nflux,unerr=normnadem.nerr,$
                                      emflux=emflux,emul=emul,$
                                      _extra=initnad.argsnadweq)

;             These need to be compatible with the IFSF_CMPNADWEQ defaults
              if tag_exist(initnad.argsnadweq,'emwid') then $
                 emwid=initnad.argsnadweq.emwid else emwid=20d
              if tag_exist(initnad.argsnadweq,'iabsoff') then $
                 iabsoff=initnad.argsnadweq.iabsoff else iabsoff=4l
           endif else begin
              weq = ifsf_cmpnadweq(normnad.wave,normnad.nflux,normnad.nerr,$
                                   snflux=normnadem.nflux,unerr=normnadem.nerr,$
                                   emflux=emflux,emul=emul)
;             These need to be compatible with the IFSF_CMPNADWEQ defaults
              emwid=20d
              iabsoff=4l
           endelse
;          Compute stellar continuum NaD equivalent widths from fit
           stelweq = ifsf_cmpnadweq(normnadstel.wave,normnadstel.nflux,$
                                    normnadstel.nerr,$
                                    wavelim=[5883d*(1d +initdat.zsys_gas),$
                                             6003d*(1d +initdat.zsys_gas),$
                                             0d,0d])

;          Compute empirical velocities
           size_weq = size(weq)
           if size_weq[0] eq 2 then begin
              if tag_exist(initnad,'argsnadvel') then $
                 vel = ifsf_cmpnadvel(normnad.wave,normnad.nflux,normnad.nerr,$
                                      weq[*,1],initdat.zsys_gas,$
                                      _extra=initnad.argsnadvel) $
              else vel = ifsf_cmpnadvel(normnad.wave,normnad.nflux,normnad.nerr,$
                                        weq[*,1],initdat.zsys_gas)
           endif else vel = dblarr(6)+bad   
;          Create / populate output data cube                            
           if firstnadnorm then begin
              nadcube = $
                 {wave: dblarr(cube.ncols,cube.nrows,n_elements(normnad.wave)),$
                  cont: dblarr(cube.ncols,cube.nrows,n_elements(normnad.wave)),$
                  dat: dblarr(cube.ncols,cube.nrows,n_elements(normnad.wave)),$
                  err: dblarr(cube.ncols,cube.nrows,n_elements(normnad.wave)),$
                  weq: dblarr(cube.ncols,cube.nrows,4)+bad,$
                  stelweq: dblarr(cube.ncols,cube.nrows,2)+bad,$
                  iweq: dblarr(cube.ncols,cube.nrows,4)+bad,$
                  emflux: dblarr(cube.ncols,cube.nrows,2)+bad,$
                  emul: dblarr(cube.ncols,cube.nrows,4)+bad,$
                  vel: dblarr(cube.ncols,cube.nrows,6)+bad}
              firstnadnorm = 0
           endif
           nadcube.wave[i,j,*] = normnad.wave
           nadcube.cont[i,j,*] = struct.cont_fit[normnad.ind]
           nadcube.dat[i,j,*] = normnad.nflux
           nadcube.err[i,j,*] = normnad.nerr
           nadcube.weq[i,j,*] = weq[*,0]
           nadcube.iweq[i,j,*] = weq[*,1]
           nadcube.stelweq[i,j,*] = stelweq[0:1]
           nadcube.emflux[i,j,*] = emflux
           nadcube.emul[i,j,*] = emul
           nadcube.vel[i,j,*] = vel
;          Plot data
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

nofit:

     endfor

  endfor

  free_lun,fitlun

  if ~ tag_exist(initdat,'noemlinfit') then begin
     free_lun,linlun
     save,linmaps,file=initdat.outdir+initdat.label+'.lin.xdr'
     save,tlinmaps,file=initdat.outdir+initdat.label+'.tlin.xdr'
  endif

  if tag_exist(initdat,'decompose_ppxf_fit') OR $
     tag_exist(initdat,'decompose_qso_fit') then $
     save,contcube,file=initdat.outdir+initdat.label+'.cont.xdr'
  if tag_exist(initdat,'donad') then $
     save,nadcube,file=initdat.outdir+initdat.label+'.nadspec.xdr'

end
