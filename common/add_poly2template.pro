;
; History
;  2009xxxYY  Jabran Zahid  created
;  2010apr13  DSNR          added exponentials
;
; Add polynomials (and possibly exponentials) to stellar continuum
; libraries for spectral fitting.
;

function add_poly2template, template, nterms, addexp=addexp

  if nterms eq -1 then return, template

  el = double(n_elements(template[*,0]))
  new_temp = template

  for i=0, nterms do $
     new_temp = [[new_temp], [(dindgen(el)/el)^i] ]

  if keyword_set(addexp) then begin
     new_temp = [[new_temp], [exp(-dindgen(el)/el)] ]
     new_temp = [[new_temp], [exp(-2d*dindgen(el)/el)] ]
     new_temp = [[new_temp], [exp(-4d*dindgen(el)/el)] ]
     new_temp = [[new_temp], [exp(-8d*dindgen(el)/el)] ]
     new_temp = [[new_temp], [1d -exp(-dindgen(el)/el)] ]
     new_temp = [[new_temp], [1d -exp(-2d*dindgen(el)/el)] ]
     new_temp = [[new_temp], [1d -exp(-4d*dindgen(el)/el)] ]
     new_temp = [[new_temp], [1d -exp(-8d*dindgen(el)/el)] ]
  endif

  return, new_temp

end
