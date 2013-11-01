;+
; NAME:
;       COMBINE_BLENDS
;
; PURPOSE:
;       Combine doublet line flux measurements (e.g., [OII]) into a
;       single line measurement. 
;
; CALLING SEQUENCE:
;       sout = combine_blends(s)
;
; INPUTS:
;       s    - linefit structure from IFITSPEC
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       sout - s (modified)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       So far only supports the [O II] doublet.
;
; PROCEDURES USED:
;       STRUCT_TRIMTAGS(), STRUCT_ADDTAGS()
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 Feb 27, U of A, written
;       jm03jun12uofa - major rewrite
;-

function combine_blends, s
    
    known = ['OII_3727']
    nknown = n_elements(known)

    idtags = ['','_BOX','_EW','_CHI2','_SIGMA','_CONTINUUM','_WAVE']
    ntags = n_elements(idtags)

    linename = s.linename
    nline = n_elements(linename)

    headstruct = struct_trimtags(s,select=['Z_LINE*','LINENAME'])
    headstruct = struct_trimtags(headstruct,except=['LINENAME'])
    
    for k = 0L, nknown-1L do begin
    
       tags   = known[k]+idtags      ; desired output tags
       tags_1 = known[k]+'_1'+idtags ; first line
       tags_2 = known[k]+'_2'+idtags ; second line

       s1 = struct_trimtags(s,select=tags_1)
       s2 = struct_trimtags(s,select=tags_2)
       sclean = struct_trimtags(s,except=['Z_LINE*','LINENAME',tags_1,tags_2])
       
       init = create_struct($
         tags[0], fltarr(2), $
         tags[1], fltarr(2), $
         tags[2], fltarr(2), $
         tags[3], -99.0,     $
         tags[4], fltarr(2), $
         tags[5], fltarr(2), $
         tags[6], 0.0)

; mean wavelength

       init.(6) = djs_mean([s1.(6),s2.(6)])
       
; if one or the either line was "not measured", then make the combined
; line "not measured".  if one or both lines are computed "upper
; limits", then make the combined line an "upper limit"
       
       errflag = 0.0                                                         ; measured well
       if ((s1.(0))[1] eq -2.0) or ((s2.(0))[1] eq -2.0) then errflag = -2.0 ; not measured
       if ((s1.(0))[1] eq -3.0) or ((s2.(0))[1] eq -3.0) then errflag = -3.0 ; upper limit
       
       case errflag of

          0.0: begin ; both lines are well-measured

; linearly add the Gaussian flux and add the error in quadrature
             init.(0) = [(s1.(0))[0]+(s2.(0))[0],sqrt((s1.(0))[1]^2+(s2.(0))[1]^2)]             ; Gaussian flux
; linearly add the box flux and add the error in quadrature
             init.(1) = [(s1.(1))[0]+(s2.(1))[0],sqrt((s1.(1))[1]^2+(s2.(1))[1]^2)]             ; box flux
; linearly add the Gaussian EW and add the error in quadrature
             init.(2) = [(s1.(2))[0]+(s2.(2))[0],sqrt((s1.(2))[1]^2+(s2.(2))[1]^2)]             ; EW
; assign the chi2 to be the mean chi2
             init.(3) = djs_mean([s1.(3),s2.(3)])                                               ; chi2
; add the sigma line widths and errors in quadrature
             init.(4) = [sqrt((s1.(2))[0]^2+(s2.(2))[0]^2),sqrt((s1.(2))[1]^2+(s2.(2))[1]^2)]   ; sigma
; compute the mean continuum and add the errors in quadrature
             init.(5) = [djs_mean([(s1.(5))[0],(s2.(5))[0]]),sqrt((s1.(5))[1]^2+(s2.(5))[1]^2)] ; continuum
             
          end
          
          -2.0: begin ; set the combined line as not measured as well

             init.(0) = [0.0,-2.0] ; Gaussian flux
             init.(1) = [0.0,-2.0] ; box flux
             init.(2) = [0.0,-2.0] ; EW
             init.(3) = -2.0       ; chi2
             init.(4) = [0.0,-2.0] ; sigma
             init.(5) = [0.0,-2.0] ; continuum
             
          end

          -3.0: begin ; set the combined Gaussian flux and EW as upper limits

; linearly add the Gaussian flux and set the error to -3.0
             init.(0) = [(s1.(0))[0]+(s2.(0))[0],-3.0]              ; Gaussian flux
; the box flux is not meaningful
             init.(1) = [0.0,-3.0]                                 ; box flux
; linearly add the Gaussian EW and set the error to -3.0
             init.(2) = [(s1.(2))[0]+(s2.(2))[0],-3.0]             ; EW
; set the chi2 to -3.0
             init.(3) = -3.0                                       ; chi2
; add the sigma line widths and set the error to -3.0
             init.(4) = [sqrt((s1.(2))[0]^2+(s2.(2))[0]^2),-3.0]   ; sigma
; compute the mean continuum and set the error to -3.0
             init.(5) = [djs_mean([(s1.(5))[0],(s2.(5))[0]]),-3.0] ; continuum
             
          end

          else: message, 'Problem here!'

       endcase

       keep = lindgen(nline)
       remove, where(known[k]+['_1','_2'] eq linename), keep
       newlinename = [known[k],linename[keep]]

       sout = struct_addtags(struct_addtags(struct_addtags(headstruct,{linename: newlinename}),init),sclean)
       
    endfor

return, sout
end
