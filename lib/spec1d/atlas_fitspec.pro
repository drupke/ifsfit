;+
; NAME:
;       ATLAS_FITSPEC()
;
; PURPOSE:
;       Fit all the galaxies in the K/M spectral atlas. 
;
; CALLING SEQUENCE:
;       specdata = atlas_fitspec(speclist,[datapath=,eigenspec=,$
;          eigendir=,elinefile=,alinefile=,linepath=,snrcut=,psname=,$
;          suffix=,specfit=],/doplot,/debug,/postscript,/write,_extra=extra)
;
; INPUTS:
;       speclist - FITS file list of galaxy spectra to fit
;
; OPTIONAL INPUTS:
;       datapath  - I/O path
;       eigenspec - spectral templates to use in the fitting
;          0 - jacoby_templates.fits (default)
;          1 - bc03_020_templates.fits (BC03 solar metallicity)
;          2 - bc03_008_templates.fits (BC03 LMC metallicity)
;       eigendir  - path to EIGENSPEC (default ${ISPEC_DIR}/templates)
;       elinefile - text file listing which emission lines to fit
;                   (default 'elinelist.dat')
;       alinefile - text file listing which absorption lines to fit
;                   (default 'alinelist.dat')
;       linepath  - path to ELINEFILE and ALINEFILE (default
;                   ${ISPEC_DIR}/etc) 
;       snrcut    - compute upper limits on lines with S/N < SNRCUT
;                   (default 3.0)
;       psname    - postscript output file name if POSTSCRIPT=1
;                   (default MJDSTR+'_atlas_specfit.ps', where MJDSTR
;                   is the modified Julian date)
;       suffix    - unique output file name suffix (default 'atlas') 
;
; KEYWORD PARAMETERS:
;       doplot     - generate a three windows (pages) of plots of the
;                    fitting results (if POSTSCRIPT=1 then DOPLOT=0
;                    always)
;       debug      - enable additional diagnostic plots in IFITSPEC
;       postscript - generate postscript output of the fitting results
;                    (see PSNAME)
;       write      - write the results of the spectral line and
;                    continuum fitting as two binary FITS files (see
;                    COMMENTS) 
;
; OUTPUTS:
;       specdata   - line-fitting results (fluxes, EWs, widths, etc.) 
;
; OPTIONAL OUTPUTS:
;       specfit    - continuum and continuum error spectra and
;                    emission-line spectrum [NPIX,3,NSPEC array] 
;
; COMMENTS:
;
; DATA FILES:
;       ${ISPEC_DIR}/etc/elinelist.dat
;       ${ISPEC_DIR}/etc/alinelist.dat
;       ${ISPEC_DIR}/etc/jacoby_templates.fits
;       ${ISPEC_DIR}/etc/bc03_020_templates.fits
;       ${ISPEC_DIR}/etc/bc03_008_templates.fits
;
; PROCEDURES USED:
;       CWD(), SPLOG, MRDFITS(), MAKE_WAVE(), SXPAR(), MWRFITS,
;       IM_OPENCLOSE, READ_LINEPARS(), GET_JULDATE, RD1DSPEC(),
;       STRUCT_TRIMTAGS(), IFITSPEC, STRUCT_ADDTAGS(),
;       PARSE_LINEFIT(), MKHDR, SXDELPAR, SXADDPAR,
;       ATLAS_FIGURE_FITSPEC, ICLEANUP 
;
; EXAMPLE:
;       To fit a single galaxy spectrum with no extinction curve and
;       three polynomial background terms:
;
;          IDL>  specdata = atlas_fitspec('galaxy.fits',/nodust,$
;                nback=3,/postscript,/write)
;
; MODIFICATION HISTORY:
;       J. Moustakas, Spring 2003, U of A
;-

function atlas_fitspec, speclist, datapath=datapath, eigenspec=eigenspec, $
  eigendir=eigendir, elinefile=elinefile, alinefile=alinefile, $
  linepath=linepath, snrcut=snrcut, psname=psname, suffix=suffix, specfit=specfit, $
  doplot=doplot, debug=debug, postscript=postscript, write=write, _extra=extra

    nspec = n_elements(speclist)
    if nspec eq 0L then begin
       print, 'Syntax - specdata = atlas_fitspec(speclist,[datapath=,eigenspec=,$'
       print, '   eigendir=,elinefile=,alinefile=,linepath=,snrcut=,psname=,$'
       print, '   suffix=,specfit=],/doplot,/debug,/postscript,/write,_extra=extra)'
       return, -1 
    endif

    light = 2.99792458D5        ; speed of light [km/s]

    if n_elements(datapath) eq 0L then datapath = cwd()
    
; specify the eigentemplates

    if n_elements(eigendir) eq 0L then eigendir = filepath('',$
      root_dir=getenv('ISPEC_DIR'),subdirectory='templates')

    if n_elements(eigenspec) eq 0L then eigenspec = 0L
    case eigenspec of
       0L: begin
          eigenfile = 'jacoby_templates.fits'
          templates = 'Jacoby Stellar Templates'
       end
       1L: begin
          eigenfile = 'bc03_020_templates.fits' ; BC03 solar metallicity
          templates = 'Bruzual & Charlot 2003 (Solar Metallicity)'
       end
       2L: begin
          eigenfile = 'bc03_008_templates.fits' ; BC03 LMC metallicity
          templates = 'Bruzual & Charlot 2003 (LMC Metallicity)'
       end
       else: begin
          eigenfile = 'jacoby_templates.fits'
          templates = 'Jacoby Stellar Templates'
       endelse
    endcase

    if file_test(eigendir+eigenfile,/regular) eq 0L then begin
       splog, 'Eigentemplate file '+eigendir+eigenfile+' not found.'
       return, -1L
    endif

; open the log file

    if n_elements(suffix) eq 0L then suffix = 'atlas'
       
    get_juldate, jd
    mjdstr = string(long(jd-2400000L),format='(I5)')

    if keyword_set(write) then begin

       logfile = mjdstr+'_'+suffix+'.log'
       specfitfile = mjdstr+'_'+suffix+'_specfit.fits'
       specdatafile = mjdstr+'_'+suffix+'_specdata.fits'

       splog, filename=logfile
       splog, 'Log file '+logfile+' opened '+systime()

    endif

    splog, 'Current datapath is ', datapath
    
; read the eigen-templates

    splog, 'Reading '+eigendir+eigenfile+'.'
    eigenflux = mrdfits(eigendir+eigenfile,0,eigenhead,/silent)
    eigenivar = eigenflux*0.0+1.0
    eigenwave = make_wave(eigenhead)
    eigenres = sxpar(eigenhead,'EIGENRES') ; [Angstrom]
    
    nstar = sxpar(eigenhead,'NAXIS2')
    splog, 'Fitting with '+strn(nstar)+' templates.'

; read the emission and absorption line parameters.  for the emission
; lines: on the first iteration constrain the widths and redshifts of
; all the lines.  on the second iteration use the constraints
; specified in ELINELIST

    if n_elements(linepath) eq 0L then $  
      linepath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='etc')
    if n_elements(elinefile) eq 0L then elinefile = 'elinelist.dat'
    if n_elements(alinefile) eq 0L then alinefile = 'alinelist.dat'

    linepars = read_linepars(linepath=linepath,linefile=elinefile)
    nline = n_elements(linepars)
    
    zindex = linepars.zindex
    windex = linepars.windex
    
    linepars.zindex = replicate('zemission',nline)
    linepars.windex = replicate('wemission',nline)

    abslines = read_linepars(linepath=linepath,linefile=alinefile)

; we will store the results of the continuum fitting for each galaxy
; in SPECLIST in a single extension of ?????_specfit.fits.
; consequently in the zeroth extension we generate a table of contents
; (toc) structure describing the contents of each extension

; table of contents structure for ?????_SPECFIT.FITS
    
    toc = {$
      date:      im_today(), $
      mjdstr:    mjdstr, $
      nspec:     nspec, $
      eigenfile: eigenfile, $
      eigendir:  eigendir, $
      templates: templates, $
      ap00:      'continuum', $
      ap01:      'continuum error', $
      ap02:      'emission line spectrum fit'}

    if keyword_set(write) then begin

       splog, 'Creating '+specfitfile+'.'
       mwrfits, toc, specfitfile, /create

    endif

    if n_elements(psname) eq 0L then psname = mjdstr+'_'+suffix+'_specfit.ps'

    if keyword_set(postscript) then begin
       doplot = 0
       splog, 'Opening plot file '+psname
       im_openclose, psname, /postscript, /silent
    endif
 
; ---------------------------------------------------------------------------
; loop on each spectrum in SPECLIST
; ---------------------------------------------------------------------------
    
    splog, 'IDL version: ' + string(!version,format='(99(a," "))')
    spawn, 'uname -a', uname
    splog, 'UNAME: ' + uname[0]
    
    stime0 = systime(1)

    for iobj = 0L, nspec-1L do begin

       splog, format='("Fitting object #",I0,"/",I0,".")', iobj+1, nspec
       
; read the spectrum and the redshift
       
       scube = rd1dspec(speclist[iobj],datapath=datapath)
       flux = scube.spec
       ferr = scube.sigspec
       wave = scube.wave
       invvar = 1.0/ferr^2.0
       npix = scube.npix
       header = scube.header

; data spectral resolution
    
       specres = replicate(7.6,npix)

; extract redshift information on GALAXY
       
       galaxy = strcompress(sxpar(header,'GALAXY'),/remove)
       z_obj = float(sxpar(header,'Z',count=zcount))
       if zcount eq 0L then begin
          splog, 'WARNING: No redshift information for '+speclist[iobj]+'.'
          z_obj = 0.0
       endif

       z_obj_err = float(sxpar(header,'Z_ERR',count=zcount))
       if (zcount eq 0L) then z_err = float(-999.0)
       
; information structure to prepend to SPECDATA
       
       basicinfo = {galaxy: galaxy, specfile: speclist[iobj], $
         specres: float(djs_mean(specres)), fit_id: iobj, $
         z_obj: z_obj, z_obj_err: z_obj_err}

; fit the whole spectrum

       ifitspec, flux, wave, eigenflux, eigenwave, invvar=invvar, $
         zobj=z_obj, specres=specres, eigenres=eigenres, snrcut=snrcut, $
         linepars=linepars, abslines=abslines, zindex=zindex, windex=windex, $
         findex=findex, fvalue=fvalue, backfit=backfit, linefit=linefit, $
         absfit=absfit, indices=indices, speclinefit=speclinefit, $
         speclineewfit=speclineewfit, starflux=starflux, doplot=debug, $
         /combine_blends, _extra=extra

; parse the results: append the absorption-line measurements, Lick
; indices, and continuum-fitting results and emission-line
; measurements

       if n_elements(abslines) ne 0L then begin

          abschi2 = absfit.linechi2
          absgood = where(abschi2 gt 0.0,nabsgood)
          if nabsgood ne 0L then abschi2[absgood] = abschi2[absgood]/(absfit[absgood].linedof-1)
       
          abs = create_struct($
           'ABS_'+absfit[0].linename+'_EW',[absfit[0].lineew_area,absfit[0].lineew_area_err],$
           'ABS_'+absfit[0].linename+'_CHI2',abschi2[0])
          for k = 1L, n_elements(absfit)-1L do abs = create_struct($
           abs,'ABS_'+absfit[k].linename+'_EW',[absfit[k].lineew_area,absfit[k].lineew_area_err],$
           'ABS_'+absfit[k].linename+'_CHI2',abschi2[k])

          abs = struct_addtags({abs_linename: absfit.linename},abs)
          temp = struct_addtags(abs,indices)
          
       endif else temp = indices
          
       backinfo = struct_trimtags(backfit,select=['CONTINUUM_CHI2','EBV','STARCOEFF','Z_ABS','Z_ABS_ERR'])
       plinefit = parse_linefit(linefit,specres=basicinfo.specres,snrcut=snrcut)
       
       specdata1 = struct_addtags(basicinfo,struct_addtags(struct_addtags(temp,backinfo),plinefit))
          
       if (iobj eq 0L) then specdata = specdata1 else specdata = [ [specdata], [specdata1] ]
       
; initialize the output FITS vector and header

       specfit = float([ [backfit.continuum], [backfit.continuum_sigma], [speclinefit] ])
;      specfit = float([ [backfit.continuum], [backfit.continuum_sigma], [speclinefit], [speclineewfit] ])
       
       mkhdr, spechead, specfit, /extend
       sxdelpar, spechead, 'COMMENT'
       sxdelpar, spechead, 'DATE'
       sxaddpar, spechead, 'GALAXY', galaxy
       sxaddpar, spechead, 'FILE', speclist[iobj]

; generate a figure of the fit

       if keyword_set(doplot) or keyword_set(postscript) then begin

          if keyword_set(postscript) then splog, 'Generating postscript for '+strn(speclist[iobj])+'.'

          atlas_figure_fitspec, specdata1, specfit, wave=wave, flux=flux, toc=toc, postscript=postscript

          if (keyword_set(doplot) and (nspec gt 1L)) then begin
             splog, 'Press any key to continue.'
             cc = get_kbrd(1)
          endif

       endif
          
       if keyword_set(write) then begin

          splog, 'Updating '+specfitfile+'.'
          mwrfits, specfit, specfitfile, spechead          

       endif 

; clean up memory

       abs = 0 & absfit = 0 & backfit = 0 & backfit = 0 & backinfo = 0 
       basicinfo = 0 & flux = 0 & ferr = 0 & invvar = 0 & linefit = 0 
       plinefit = 0 & specdata1 = 0 & indices = 0 & speclinefit = 0
       speclineewfit = 0 & specres = 0 & wave = 0 

       icleanup, scube

    endfor 

    splog, format='("Total time for ATLAS_FITSPEC = ",G0," minutes.")', (systime(1)-stime0)/60.0
    
; ---------------------------------------------------------------------------
; end main loop
; ---------------------------------------------------------------------------

; close and GZIP files
    
    specdata = reform(specdata)
    if keyword_set(write) then begin
       mwrfits, specdata, specdatafile, /create
       spawn, ['gzip -f '+specdatafile], /sh
       spawn, ['gzip -f '+specfitfile], /sh
    endif

    if keyword_set(postscript) then begin
       im_openclose, postscript=postscript, /close
       spawn, ['gzip -f '+psname], /sh
    endif

return, specdata
end
