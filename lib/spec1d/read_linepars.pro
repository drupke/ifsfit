;+
; NAME:
;       READ_LINEPARS()
;
; PURPOSE:
;
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;
;
; PROCEDURES USED:
;
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 October 31, U of A, written
;       jm03jul3uofa - fully documented
;-

function read_linepars, linepath=linepath, linefile=linefile

    if n_elements(linepath) eq 0L then $  
      linepath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='etc')
    if n_elements(linefile) eq 0L then linefile = 'elinelist.dat'
    if file_test(linepath+linefile,/regular) eq 0L then message, 'Line list file '+$
      linepath+linefile+' not found.'

    readcol, linepath+linefile, line, center, zindex, windex, findex, fvalue, blend, $
      format = 'A,F,A,A,A,F,A', /silent, comment='#'

    linepars = {line: ' ', wave: 0.0, zindex: '', windex: '', findex: '', fvalue: 0.0, blend: ''}
    nline = n_elements(line)
    linepars = replicate(linepars,nline)

    linepars.line = line
    linepars.wave = center
    linepars.zindex = zindex
    linepars.windex = windex
    linepars.findex = findex
    linepars.fvalue = fvalue
    linepars.blend = blend

return, linepars
end    
