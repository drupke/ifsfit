;+
; NAME:
;       FILL_EW()
;
; PURPOSE:
;       Compute the equivalent width from a measure of the continuum
;       and total line flux. 
;
; INPUTS:
;       linefit - output structure from ILINEFIT() or IABSLINEFIT() 
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       linefit - (modified)
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       IM_COMPUTE_ERROR()
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 June 12, U of A, excised from IFITSPEC 
;-

function fill_ew, linefit

    nline = n_elements(linefit)
    if nline eq 0L then begin
       print, 'Syntax - linefit = fill_ew(linefit)'
       return, -1L
    endif

    for k = 0L, nline-1L do begin

       case linefit[k].linearea_err of

          -2.0: ; not measured

          -3.0: begin ; upper limit

             linefit[k].lineew_area = linefit[k].linearea / linefit[k].linecontlevel
             linefit[k].lineew_area_err = -3.0

             linefit[k].lineew_box = 0.0
             linefit[k].lineew_box_err = -3.0

          end

          else: begin
       
             linefit[k].lineew_area = linefit[k].linearea / linefit[k].linecontlevel
             linefit[k].lineew_area_err = im_compute_error(linefit[k].linearea,$
               linefit[k].linearea_err,linefit[k].linecontlevel,linefit[k].linecontlevel_err,/quotient)

             linefit[k].lineew_box = linefit[k].linebox / linefit[k].linecontlevel
             linefit[k].lineew_box_err = im_compute_error(linefit[k].linebox,$
               linefit[k].linebox_err,linefit[k].linecontlevel,linefit[k].linecontlevel_err,/quotient)

          end
             
       endcase

    endfor
    
return, linefit
end

