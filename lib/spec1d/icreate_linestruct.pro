function icreate_linestruct, nline

; initialize the output structure, and set default return values 

   linestruct = create_struct($
    'linename'          ,   ' ', $
    'linewave'          ,   0.0, $
    'linez'             ,   0.0, $
    'linez_err'         ,  -1.0, $
    'linesigma'         ,   0.0, $
    'linesigma_err'     ,  -99.0,$
    'linearea'          ,   0.0, $
    'linearea_err'      ,  -99.0,$
    'linebox'           ,   0.0, $
    'linebox_err'       ,  -99.0,$
    'linecontlevel'     ,   0.0, $
    'linecontlevel_err' ,  -99.0,$
    'lineew_area'       ,   0.0, $
    'lineew_area_err'   ,  -99.0,$
    'lineew_box'        ,   0.0, $
    'lineew_box_err'    ,  -99.0,$
    'line_blend'        ,   '',  $
    'linenpix'          ,   0,   $
    'linedof'           ,   0.0, $
    'linechi2'          ,  -99.0)
   if (keyword_set(nline)) then linestruct = replicate(linestruct,nline)

return, linestruct
end
