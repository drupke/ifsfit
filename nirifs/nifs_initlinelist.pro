function nifs_initlinelist,bre=bre,hei187=hei187,hei206=hei206,si6=si6,$
                           h2s4=h2s4,h2s5=h2s5
;
; History
;   13mar05
;
; Wavelength sources:
;   1. H_2: http://www.jach.hawaii.edu/UKIRT/astronomy/calib/spec_cal/h2_s.html
;   2. He: http://www.jach.hawaii.edu/UKIRT/astronomy/calib/spec_cal/he_4.html
;   3. H: http://www.jach.hawaii.edu/UKIRT/astronomy/calib/spec_cal/recomb_menu.html
;   4. [SiVI]: NIST
;   Checked He and H lines in NIST and there is good agreement.
; 

  

  label = ['H2_10_S3','H2_10_S2','H2_10_S1',$
           'Paa','Brd','Brg']
  if keyword_set(bre) then label = [label,'Bre']
  if keyword_set(hei187) then label = [label,'HeI187']
  if keyword_set(hei206) then label = [label,'HeI206']
  if keyword_set(si6) then label = [label,'Si6']
  if keyword_set(h2s4) then label = [label,'H2_10_S4']
  if keyword_set(h2s5) then label = [label,'H2_10_S5']

; Vacuum wavelengths!
  wave = double([19576,20338,21218,$
                 18756,19451,21661])
  if keyword_set(bre) then wave = [wave,double(18181)]
  if keyword_set(hei187) then wave = [wave,double(18691)]
  if keyword_set(hei206) then wave = [wave,double(20587)]
  if keyword_set(si6) then wave = [wave,double(19646)]
  if keyword_set(h2s4) then wave = [wave,double(18920)]
  if keyword_set(h2s5) then wave = [wave,double(18358)]

  linelist = {label:label, wave:wave}

  return,linelist

end
