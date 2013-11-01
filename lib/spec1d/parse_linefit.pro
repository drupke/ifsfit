function parse_linefit, linefit, specres=specres, snrcut=snrcut, $
  combineblends=combineblends
; LINEFIT is one structure per emission line
; initialize the output structure with the line fluxes and add a
; galaxy, specfile, and a redshift tag 

    light = 2.99792458D5        ; speed of light [km/s]

    if n_elements(specres) eq 0L then specres = 5.0 ; [Angstrom]
    if n_elements(snrcut) eq 0L then snrcut = 3.0
    nline = n_elements(linefit.linename)

    linename = strcompress(strupcase(linefit.linename),/remove)
    linewave = linefit.linewave

    for i = 0L, nline-1L do begin

; desired output tags.  also listed in COMBINE_BLENDS()
       
       tags = linename[i]+['','_BOX','_EW','_CHI2','_SIGMA','_CONTINUUM','_WAVE']

; initialize the structure fields for each emission line fitted
       
       out1 = create_struct($
         tags[0], [linefit[i].linearea,linefit[i].linearea_err], $
         tags[1], [linefit[i].linebox,linefit[i].linebox_err], $
         tags[2], [linefit[i].lineew_area,linefit[i].lineew_area_err], $
         tags[3], -1.0, $
         tags[4], [linefit[i].linesigma,linefit[i].linesigma_err], $
         tags[5], [linefit[i].linecontlevel,linefit[i].linecontlevel_err], $
         tags[6], linewave[i])

; compute the reduced chi2.  N.B. the out1.(3) is not general if TAGS
; changes 
       
       if (linefit[i].linedof gt 0.0) then out1.(3) = linefit[i].linechi2/(linefit[i].linedof-1)

; remove the instrumental resolution from the SIGMA line width for
; good lines

       if (linefit[i].linearea_err gt 0.0) then begin
          sigmafit = (out1.(4))[0]
          sigmains = light*specres/out1.(6)
          if sigmains lt sigmafit then out1.(4) = [sqrt(sigmafit^2.0 - sigmains^2.0),(out1.(4))[1]]
       endif
       
; grow the output structure       
       
       if (i eq 0L) then out = out1 else out = struct_addtags(out,out1)
       
    endfor

; add the mean redshift

    good = where(linefit.linearea_err gt 0.0,ngood)
    if ngood ne 0L then begin
       z_line = djs_mean(linefit[good].linez)   ; mean redshift
       z_line_err = djsig(linefit[good].linez) ; standard error
    endif else begin
       z_line = 0.0
       z_line_err = -1.0
    endelse

    parsed = struct_addtags({z_line: float(z_line), z_line_err: $
      float(z_line_err), linename: linename},out)

    if keyword_set(combineblends) then parsed = combine_blends(parsed)
    
return, parsed
end
