;
; History
;  10jul21  DSNR  created
;

pro uhsf_readnadfit,parfile,abspars,empars,opars

  openr,lun,parfile,/get_lun
  line = ''
  readf,lun,line
  readf,lun,line
  linespl = strsplit(line,' ',/extract)
  nabs = fix(linespl[0])
  nem = fix(linespl[1])
  llo = double(linespl[2])
  lhi = double(linespl[3])
  readf,lun,line
  readf,lun,line
  readf,lun,line
  if nabs gt 0 then begin
     abspars = dblarr(4,nabs)
     for i=0,nabs-1 do begin
        for j=0,3 do begin
           readf,lun,line
           linespl = strsplit(line,' ',/extract)
           abspars[j,i] = double(linespl[4])
        endfor
     endfor
  endif
  if nem gt 0 then begin
     empars = dblarr(3,nem)
     for i=0,nem-1 do begin
        for j=0,2 do begin
           readf,lun,line
           linespl = strsplit(line,' ',/extract)
           empars[j,i] = double(linespl[4])
        endfor
     endfor
  endif
  free_lun,lun

  opars = {nabs: nabs,$
           nem:  nem,$
           llo:  llo,$
           lhi:  lhi $
          }
           

end
