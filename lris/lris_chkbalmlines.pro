function lris_chkbalmlines,fluxes,rv=rv,ebmv=ebmv

;
; History
;   09may26 DSNR created
;   09jul   DSNR re-written
;   09aug25 DSNR added calculation for total flux over components
;   09sep08 DSNR changed from Calzetti 2000 to CCM89 reddening curve,
;                at Calzetti's advice (via Lisa): don't use C00 for
;                emission lines.
;

  if ~ keyword_set(rv) then rv = 3.1

  bad = 999d

; Sum over components to get Balmer correction
  fluxes_info = size(fluxes.value)
  if fluxes_info[0] gt 1 then begin
     ncomp = fluxes_info[2]
     fluxvals = [[total(fluxes.value,2)],[fluxes.value]]
     fluxerrs = [[sqrt(total((fluxes.error)^2,2))],[fluxes.error]]
  endif else begin
     ncomp=1
     fluxvals = [[fluxes.value],[fluxes.value]]
     fluxerrs = [[fluxes.error],[fluxes.error]]
  endelse

  nlines = n_elements(fluxes.value)/ncomp
  
; Unextincted Balmer ratios, from Brocklehurst 1971, MNRAS, 153, 471
; T = 10^4 K, N_e = 10^2 cm^-2

;  HxHb_th = [2.859d,1d,0.469d,0.259d,0.159d,0.105d,0.073d,$
;             0.053d,0.040d,0.031d]
;  lambda = [6563d,4871d,4341d,4104d,3970d,3889d,3835d,3798d,3771d,3750d]
  HxHb_th = [2.859d,1d,0.469d,0.259d,0.159d]
  lambda = [6563d,4861d,4341d,4104d,3970d]
  alamav = extcurve_ccm(lambda,rv=rv)

; Balmer line ratios
  HxHb = dblarr(nlines-1,2,ncomp+1)-bad
; E(B-V)
  ebmv = dblarr(nlines-1,2,ncomp+1)-bad
; Predicted Ha/Hb based on other Balmer lines
  HaHb_pred = dblarr(nlines-1,2,ncomp+1)-bad
; Predicted Ha/Hb based on other Balmer lines divided by actual Ha/Hb;
; Factor by which to multiply blue LRIS data to correct it.
  dHaHb = dblarr(nlines-1,2,ncomp+1)-bad

  for j=0,ncomp do begin

;    DO HA/HB FIRST
  
;    Compute observed Ha/Hb
     HxHb[0,0,j] = fluxvals[0,j]/fluxvals[1,j]
     HxHb[0,1,j] = HxHb[0,0,j]*sqrt((fluxerrs[0,j]/fluxvals[0,j])^2+$
                                    (fluxerrs[1,j]/fluxvals[1,j])^2)
;    Compute E(B-V) from Ha/Hb
     this_ebmv = ebv_ccm([lambda[0],lambda[1]],$
                         [fluxvals[0,j],fluxvals[1,j]],HxHb_th[0],$
                         rv=rv,fluxerr=[fluxerrs[0,j],fluxerrs[1,j]])
     ebmv[0,0,j] = this_ebmv[0]
     ebmv[0,1,j] = this_ebmv[1]

;    Predicted Ha/Hb based on this E(B-V) and corresponding blue/red
;    corr.
     HaHb_pred[0,0,j] = HxHb[0,0,j]
     HaHb_pred[0,1,j] = 1d
     dHaHb[0,0,j] = 1d
     dHaHb[0,1,j] = 0d
;    Set correction to less than one if measured extinction is less
;    than 0.
     if ebmv[0,0,j] lt 0d then $
        dHaHb[0,0,j] = HxHb[0,0,j] / HxHb_th[0]

;    NOW THE OTHER BALMER LINE RATIOS     
    
     for i=1,nlines-2 do begin
  
        if (fluxvals[i+1,j] gt 0) then begin

           HxHb[i,0,j] = fluxvals[i+1,j]/fluxvals[1,j]
           HxHb[i,1,j] = HxHb[i,0,j]*sqrt((fluxerrs[i+1,j]/fluxvals[i+1,j])^2+$
                                          (fluxerrs[i,j]  /fluxvals[i,j])^2)
           this_ebmv = ebv_ccm([lambda[i+1],lambda[1]],$
                               [fluxvals[i+1,j],fluxvals[1,j]],HxHb_th[i+1],$
                               rv=rv,fluxerr=[fluxerrs[i+1,j],fluxerrs[1,j]])
           ebmv[i,0,j] = this_ebmv[0]
           ebmv[i,1,j] = this_ebmv[1]

;          Set effective extinction to zero if it comes out less
           usebmv = ebmv[i,0,j]
           if usebmv lt 0d then usebmv = 0d

           HaHb_pred[i,0,j] = 10d^(-0.4d*(usebmv*rv*(alamav[0]-alamav[1])-$
                                          2.5d*alog10(HxHb_th[0])))
           HaHb_pred[i,1,j] = -HaHb_pred[i,0,j]*0.4d*rv*(alamav[0]-alamav[1])*$
                              alog(10d)*ebmv[i,1,j]
           dHaHb[i,0,j] = HxHb[0,0,j] / HaHb_pred[i,0,j]
           dHaHb[i,1,j] = dHaHb[i,0,j]*sqrt((HxHb[0,1,j]/HxHb[0,0,j])^2+$
                                            (HaHb_pred[i,1,j]/HaHb_pred[i,0,j])^2)
           
           if dHaHb[i,0,j] lt 0 then begin
              dHaHb[i,0,j] = -999d
              dHaHb[i,1,j] = -999d
           endif

        endif

     endfor

  endfor

  return,dHaHb
  
end
