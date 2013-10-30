;
; History
;  13mar05  DSNR  created
;
; Set (flux errors eq 0) as a criterion b/c if line is specified to be
; fit but out of the fit range, it can yield a nonzero flux and 0
; error.
;

function nifs_checkcomp,struct,z,dblsigh2=dblsigh2,singsigh2=singsigh2,$
                        paasig=paasig,siglim=siglim,doubleline=doubleline,$
                        use_h2_10_s2=use_h2_10_s2
                        
; Sigma = 1.6 A corresponds to R = 5290 (from NIFS website) at 2um
; Sigma = 85 A corresponds to 3000 km/s FWHM at 2um
  if ~ keyword_set(siglim) then siglim_use=[1.6d,85d] else siglim_use = siglim
  allowsigpeg = [1,1]
  if siglim_use[0] lt 0 then begin
     siglim_use[0] = abs(siglim_use[0])
     allowsigpeg[0] = 0
  endif
  if siglim_use[1] lt 0 then begin
     siglim_use[1] = abs(siglim_use[1])
     allowsigpeg[1] = 0
  endif
  if ~ keyword_set(paasig) then paasig = 3d
  if ~ keyword_set(dblsigh2) then dblsigh2 = 3d
  if ~ keyword_set(singsigh2) then singsigh2 = 5d

  ncomp = struct.param[1]

  linepars = sepfitpars(struct.param,struct.perror)
  nlines = n_elements(linepars.flux[*,0])

  goodcomp = -1

; H II cuts
  paind = where(struct.linelabel eq 'Paa')
  pasig = linepars.sigma[paind,0]
  paflux = linepars.fluxpk[paind,0]
  pafluxerr = linepars.fluxpkerr[paind,0]
  if allowsigpeg[0] then pasiglo_gd = pasig ge siglim_use[0] $
  else pasiglo_gd = pasig gt siglim_use[0]
  if allowsigpeg[1] then pasighi_gd = pasig le siglim_use[1] $
  else pasighi_gd = pasig lt siglim_use[1]
  if (paflux gt paasig*pafluxerr AND $
      pafluxerr gt 0 AND $
      pasiglo_gd AND pasighi_gd) then goodcomp = 0

; H_2 cuts. First if statement requires that line be detected in Paa.
  if goodcomp eq 0 AND ncomp gt 1 then begin
     h21ind = where(struct.linelabel eq 'H2_10_S1')
     h22ind = where(struct.linelabel eq 'H2_10_S2')
     h23ind = where(struct.linelabel eq 'H2_10_S3')
     i = 1
     h21sig = linepars.sigma[h21ind,i]
     h21flux = linepars.fluxpk[h21ind,i]
     h21fluxerr = linepars.fluxpkerr[h21ind,i]
     h22sig = linepars.sigma[h22ind,i]
     h22flux = linepars.fluxpk[h22ind,i]
     h22fluxerr = linepars.fluxpkerr[h22ind,i]
     h23sig = linepars.sigma[h23ind,i]
     h23flux = linepars.fluxpk[h23ind,i]
     h23fluxerr = linepars.fluxpkerr[h23ind,i]
     keep = 0
;       Require that each H_2 component appear significantly in at
;       least 2 lines.
     h21gd = 0
     h22gd = 0
     h23gd = 0
     h21vgd = 0
     h22vgd = 0
     h23vgd = 0
     if allowsigpeg[0] then h2siglo_gd = h21sig ge siglim_use[0] $
     else h2siglo_gd = h21sig gt siglim_use[0]
     if allowsigpeg[1] then h2sighi_gd = h21sig le siglim_use[1] $
     else h2sighi_gd = h21sig lt siglim_use[1]
     if (h21flux gt dblsigh2*h21fluxerr AND $
         h21fluxerr gt 0 AND $
         h2siglo_gd AND h2sighi_gd) then h21gd=1
     if keyword_set(use_h2_10_s2[0]) then $
        if (h22flux gt dblsigh2*h22fluxerr AND $
            h22fluxerr gt 0 AND $
            h2siglo_gd AND h2sighi_gd) then h22gd=1
     if (h23flux gt dblsigh2*h23fluxerr AND $
         h23fluxerr gt 0 AND $
         h2siglo_gd AND h2sighi_gd) then h23gd=1
     if (h21flux gt singsigh2*h21fluxerr AND $
         h21fluxerr gt 0 AND $
         h2siglo_gd AND h2sighi_gd) then h21vgd=1
     if keyword_set(use_h2_10_s2[0]) then $
        if (h22flux gt singsigh2*h22fluxerr AND $
            h22fluxerr gt 0 AND $
            h2siglo_gd AND h2sighi_gd) then h22vgd=1
     if (h23flux gt singsigh2*h23fluxerr AND $
         h23fluxerr gt 0 AND $
         h2siglo_gd AND h2sighi_gd) then h23vgd=1
     h2gd = h21gd+h22gd+h23gd
     h2vgd = h21vgd+h22vgd+h23vgd
     if h2gd gt 1 OR h2vgd ge 1 then goodcomp = [goodcomp,i]
     if doubleline eq 'h2' AND ncomp eq 3 then begin
        i = 2
        h21sig = linepars.sigma[h21ind,i]
        h21flux = linepars.fluxpk[h21ind,i]
        h21fluxerr = linepars.fluxpkerr[h21ind,i]
        h22sig = linepars.sigma[h22ind,i]
        h22flux = linepars.fluxpk[h22ind,i]
        h22fluxerr = linepars.fluxpkerr[h22ind,i]
        h23sig = linepars.sigma[h23ind,i]
        h23flux = linepars.fluxpk[h23ind,i]
        h23fluxerr = linepars.fluxpkerr[h23ind,i]
        keep = 0
;       Require that each H_2 component appear significantly in at
;       least 2 lines.
        h21gd = 0
        h22gd = 0
        h23gd = 0
        h21vgd = 0
        h22vgd = 0
        h23vgd = 0
        if allowsigpeg[0] then h2siglo_gd = h21sig ge siglim_use[0] $
        else h2siglo_gd = h21sig gt siglim_use[0]
        if allowsigpeg[1] then h2sighi_gd = h21sig le siglim_use[1] $
        else h2sighi_gd = h21sig lt siglim_use[1]
        if (h21flux gt dblsigh2*h21fluxerr AND $
            h21fluxerr gt 0 AND $
            h2siglo_gd AND h2sighi_gd) then h21gd=1
        if keyword_set(use_h2_10_s2[1]) then $
           if (h22flux gt dblsigh2*h22fluxerr AND $
               h22fluxerr gt 0 AND $
               h2siglo_gd AND h2sighi_gd) then h22gd=1
        if (h23flux gt dblsigh2*h23fluxerr AND $
            h23fluxerr gt 0 AND $
            h2siglo_gd AND h2sighi_gd) then h23gd=1
        if (h21flux gt singsigh2*h21fluxerr AND $
            h21fluxerr gt 0 AND $
            h2siglo_gd AND h2sighi_gd) then h21vgd=1
        if keyword_set(use_h2_10_s2[1]) then $
           if (h22flux gt singsigh2*h22fluxerr AND $
               h22fluxerr gt 0 AND $
               h2siglo_gd AND h2sighi_gd) then h22vgd=1
        if (h23flux gt singsigh2*h23fluxerr AND $
            h23fluxerr gt 0 AND $
            h2siglo_gd AND h2sighi_gd) then h23vgd=1
        h2gd = h21gd+h22gd+h23gd
        h2vgd = h21vgd+h22vgd+h23vgd
        if h2gd gt 1 OR h2vgd ge 1 then goodcomp = [goodcomp,i]
     endif
  endif
  
; Paa cuts if second component is HII. First if statement requires
; that first component be detected in Paa.
  if doubleline eq 'paa' AND ncomp eq 3 then begin
     pasig = linepars.sigma[paind,2]
     paflux = linepars.fluxpk[paind,2]
     pafluxerr = linepars.fluxpkerr[paind,2]
     if allowsigpeg[0] then pasiglo_gd = pasig ge siglim_use[0] $
     else pasiglo_gd = pasig gt siglim_use[0]
     if allowsigpeg[1] then pasighi_gd = pasig le siglim_use[1] $
     else pasighi_gd = pasig lt siglim_use[1]
     if (paflux gt paasig*pafluxerr AND $
         pafluxerr gt 0 AND $
         pasiglo_gd AND pasighi_gd) then goodcomp = [0,1,2] $
     else goodcomp = [0,1]
  endif
  
  return,goodcomp

end
