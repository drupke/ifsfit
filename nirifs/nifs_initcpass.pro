;
; History
;  13jun20  DSNR  created
;
pro nifs_initcpass,gal

; Minimum number of nearest neighbors with the same number of
; components or higher 
  nearestneighbors = 2

  if gal eq 'mrk231' then begin
     rootindir = '/Users/drupke/winds/ao/red/'+gal+'/'
     rootindir += 'science/merged/'
     dx = 45
     dy = 45
     dz = 3
  endif
  
; Find number of components and best fitting redshift from first round
; of fits.
  ncomp = dblarr(dx,dy)
  zinit = dblarr(dx,dy,dz)
  readcol300,'/Users/drupke/winds/ao/specfits/'+gal+'/'+gal+'.fit.dat',$
             col_v,row_v,comp_v,fwhm,z,$
             /skip,/silent,format='(I,I,I,X,X,D,D)'        
  col_v--
  row_v--
  for nz=1,dz do begin
     em_v = where(comp_v eq nz,ct_v)
     for i=0,ct_v-1 do begin
        zinit[col_v[em_v[i]],row_v[em_v[i]],nz-1] = z[em_v[i]]
        ncomp[col_v[em_v[i]],row_v[em_v[i]]]++
     endfor
  endfor

; Array to hold number of nearest neighbors with # components >=
  ncomp_neighbors = dblarr(dx,dy)
; Reduce # components by one unless otherwise instructed
  ncomp_new = ncomp
  inz = where(ncomp gt 0)
  ncomp_new[inz]--
  ncomp_new = -ncomp_new

; Find number of neighbor spaxels (touching edge or corner) that have
; the same # of components or greater.
  for ic=1,dz do begin
     ic_xy = array_indices(ncomp,where(ncomp eq ic,ctc))
     for i=0,ctc-1 do begin
        for ix=-1,1 do begin
           for iy=-1,1 do begin
              if ncomp[ic_xy[0,i]+ix,ic_xy[1,i]+iy] ge ic then $
                 ncomp_neighbors[ic_xy[0,i],ic_xy[1,i]]++
           endfor
        endfor
     endfor
;    Maintain # of components in spaxels that have at least
;    NEARESTNEIGHBORS neighbors with same # of components or greater
;    (+1 factor b/c spaxel itself included in count).
     ic_gd = where(ncomp eq ic AND ncomp_neighbors ge nearestneighbors+1)
     ncomp_new[ic_gd] = -ic
  endfor
  
  initcpass = {zinit: zinit, $
               ncomp: ncomp_new}
  
  save,initcpass,file=rootindir+gal+'_initcpass.xdr'

end
