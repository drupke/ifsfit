function calculate_gaussian, param, wave, area = area

ggg = param[0]*exp( -(wave - param[1])^2/(2. * param[2]^2) )

area = param[2]*param[0]*sqrt(!dpi*2)

return, ggg
end

;-----------------------------------------------------------------------------

function measure_equivalent_width, str, nsig = nsig

if not keyword_set(nsig) then nsig = 2

gg   = where(str.spec ne 0)
ll   = str.wave[gg]

;tfile = '/home/group10/jabran/metal/bc_model/templates/template_structure.genx'
;restgen, file= tfile, template
;new_temp = interpol_template(ll, template.lambda, template.flux)
;new_temp = [[new_temp], [dblarr(n_elements(ll))+1]]

;spec = str.spec + str.starcoeff ## new_temp
;sp   = spec[gg]/str.norm

sp = str.spec

nlines = n_elements(str.linename)
ew     = dblarr(nlines)

for i = 0, nlines - 1 do begin

  if str.param[i*3+1] gt min(str.wave) and $
     str.param[i*3+1] lt max(str.wave) and $
     str.param[i*3+2] gt 0 then begin  

    ii = where(str.linename ne str.linename[i], count)
    mg = fltarr(n_elements(str.wave))
    gg = calculate_gaussian(str.param[i*3:i*3+2], str.wave, area = area)
    for j = 0, count - 1 do $
      mg = mg + calculate_gaussian(str.param[ii[j]*3:ii[j]*3+2], str.wave)
    ss  = sp - mg
    cen = str.param[i*3+1]
    sig = str.param[i*3+2]
    ff = where(str.wave ge cen-nsig*sig and str.wave le cen+nsig*sig, count)
;    if count gt 0 then ew[i] = total(gg[ff])/total( (ss - gg)[ff] ) $
;                  else ew[i] = 0
    if count gt 0 then ew[i] = area/total( (ss - gg)[ff] ) $
                  else ew[i] = 0

  endif

endfor

;nstr = add_tag(str, ew, 'EW')
return, ew

end
