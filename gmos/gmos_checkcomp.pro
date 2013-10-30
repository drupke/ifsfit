;
; History
;  10jul19  DSNR  modified for new GMOS data
;

function gmos_checkcomp,struct,z,sigcut=sigcut,siglim=siglim

  if ~ keyword_set(siglim) then siglim=[0.7d,10d]
  if ~ keyword_set(sigcut) then sigcut = 3d

  ncomp = struct.param[1]

  linepars = sepfitpars(struct.param,struct.perror)
  nlines = n_elements(linepars.flux[*,0])

  haind = where(struct.linelabel eq 'Halpha')
  n2ind = where(struct.linelabel eq '[NII]6583')
  goodcomp = -1

  maxhaflux = max(linepars.fluxpk[haind,*])
  maxn2flux = max(linepars.fluxpk[n2ind,*])

  for i=0,ncomp-1 do begin 
     hasig = linepars.sigma[haind,i]
     haflux = linepars.fluxpk[haind,i]
     hafluxerr = linepars.fluxpkerr[haind,i]
     n2sig = linepars.sigma[n2ind,i]
     n2flux = linepars.fluxpk[n2ind,i]
     n2fluxerr = linepars.fluxpkerr[n2ind,i]
     keep = 0
     if ((haflux gt sigcut*hafluxerr AND $
;          haflux gt comppctthresh*maxhaflux AND $
          hasig ge siglim[0] AND $
          hasig le siglim[1]) OR $
         (n2flux gt sigcut*n2fluxerr AND $
;          n2flux gt comppctthresh*maxn2flux AND $
          n2sig ge siglim[0] AND $
          n2sig le siglim[1])) $
     then keep = 1
;     if i gt 0 then begin
;        if (abs(z.gas[i]-z.gas[0]) lt 0.0002 AND $
;        abs(hasig - linepars.sigma[haind,0]) le 1.5 $
;        then keep = 0
;     endif
     if keep then begin
        if goodcomp[0] eq -1 then goodcomp = i $
        else goodcomp = [goodcomp,i]
     endif

  endfor

  return,goodcomp

end
