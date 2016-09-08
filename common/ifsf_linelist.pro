; docformat = 'rst'
;
;+
;
; A lookup table for wavelengths. Presently only air wavelengths output.
;
; Wavelength sources:
;   1. NIST ASD -- For UV, observed if present, otherwise Ritz
;   2. HeI 4686 -- found random wavelength compilation. In NIST, shows
;                  up as a complex of many fine structure lines with
;                  range of several tenths of an A.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Hash of wavelengths for specified input lines.
;
; :Params:
;    inlines: in, required, type=strarr
;      Line labels for which to retrieve wavelengths.
;
; :Keywords:
;    all: in, optional, type=byte
;      Return all lines.
; 
; :Author:
;    David S. N. Rupke::
;    Anthony To::
;      Rhodes College
;      Department of Physics
;      2000 N. Parkway
;      Memphis, TN 38104
;      drupke@gmail.com
;
; :History:
;    ChangeHistory::
;      2009jul08, DSNR, created
;      2010nov04, DSNR, wavelengths corrected from Mappings III to
;                       NIST values
;      2013sep, DSNR, added more lines for 1-slit GMOS config.
;      2013nov22, DSNR, documented, renamed, added copyright and
;                       license, and turned output into a hash
;      2014jan14, DSNR, switched to ordered hashes
;      2014jan16, DSNR, fixed one wrong label
;      2014feb26, DSNR, replaced ordered hashes with hashes
;      2015jan06, DSNR, added Mg1b lines
;      2015xxxYY, AT, added UV lines
;      2016jul12, DSNR, added line labels
;      2016aug03, DSNR, edited UV line labels
;      2016sep07, DSNR, heavily edited UV lines; added ALL keyword
;    
; :Copyright:
;    Copyright (C) 2013--2016 David S. N. Rupke, Anthony To
;
;    This program is free software: you can redistribute it and/or
;    modify it under the terms of the GNU General Public License as
;    published by the Free Software Foundation, either version 3 of
;    the License or any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;    General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program.  If not, see
;    http://www.gnu.org/licenses/.
;
;-
function ifsf_linelist,inlines,linelab=linelab,all=all

;  Associated line labels:
   lines = hash()
;  Optical lines
   lines['Halpha'] = 6562.80d
   lines['Hbeta'] = 4861.32d
   lines['Hgamma'] = 4340.47d
   lines['HeI5876'] = 5875.661d
   lines['HeI6678'] = 6678.15d
   lines['HeI7065'] = 7065.19d
   lines['HeII4686'] = 4686.7d
   lines['[NeIII]3869'] = 3868.76d
   lines['[NI]5198'] = 5197.90d
   lines['[NI]5200'] = 5200.26d
   lines['[NII]5755'] = 5754.59
   lines['[NII]6548'] = 6548.05d
   lines['[NII]6583'] = 6583.45d
   lines['[OI]5577'] = 5577.34d
   lines['[OI]6300'] = 6300.30d
   lines['[OI]6364'] = 6363.78d
   lines['[OII]3726'] = 3726.032d
   lines['[OII]3729'] = 3728.815d
   lines['[OIII]4959'] = 4958.91d
   lines['[OIII]5007'] = 5006.84d
   lines['[SII]6716'] = 6716.44d
   lines['[SII]6731'] = 6730.82d
   lines['[SIII]6312'] = 6312.06d
   lines['Mg1b5167'] = 5167.3213
   lines['Mg1b5173'] = 5172.6844
   lines['Mg1b5184'] = 5183.6043  
   lines['NaD2'] = 5889.95d
   lines['NaD1'] = 5895.92d
   lines['OH8344'] = 8344.602d
   lines['OH8399'] = 8399.160d
   lines['OH8430'] = 8430.170d
   lines['[CaV]5309'] = 5309.11d
   lines['[FeVII]5159'] = 5158.89d
   lines['[FeVII]5276'] = 5276.38d
   lines['[FeVII]5721'] = 5720.7d
   lines['[FeVII]6087'] = 6087.0d
   lines['[FeX]6375'] = 6374.51d

;  UV lines
;  Sources:
;    T = Margaret Trippe's informal COS linelist
;    P01 = Prochaska+01, Table 2
;    (https://ui.adsabs.harvard.edu/#abs/2001ApJS..137...21P/abstract)
;    S = Savaglio (http://www.pha.jhu.edu/~savaglio/gdds/absorption.dat)
;    Lower state is ground unless labeled s (*) => lower state is excited
;
;  P01 includes a handful of lines from other elements (Ar, Co, Ga, Mn) that are
;  not listed here

;  P01 lists down to Ly19
   lines['Lyalpha'] = 1215.6701d
   lines['Lybeta'] = 1025.728d
   lines['Lygamma'] = 972.517d
   lines['Lydelta'] = 949.742d
   lines['Lyepsilon'] = 937.814d
   lines['Ly6'] = 930.751d
   lines['Ly7'] = 926.249d
   lines['Ly8'] = 923.148d
   lines['Ly9'] = 920.947d
   lines['Ly10'] = 919.342d
   lines['Ly11'] = 918.125d
   
;  Inlcludes all C lines from P01 up to lambda = 1600 A
   lines['CI1139'] = 1239.792d
   lines['CI1157'] = 1157.910d
   lines['CI1277'] = 1277.245d
   lines['CI1280'] = 1280.135d ; T only. f_jk comparable to other lines
   lines['CI1328'] = 1328.834d
   lines['CI1560'] = 1560.310d
   lines['CIs1560'] = 1560.683d ; another line at 1560.708, f_jk 3x less
   lines['CII1036'] = 1036.337d
   lines['CII1334'] = 1334.532d
;  Oscillator strength 10x lower than other CIIs transition
;   lines['CIIs1335'] = 1335.663d
   lines['CIIs1335'] = 1335.708d
   lines['CIII977'] = 977.03d
;  Actually a complex of several lines from an excited state; see, e.g.,
;  http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:astro-ph/9609161
;   lines['CIIIs1176'] = 1175.7217d ; T only
;  See Feibelman & Johansson 1995
;  (http://adsabs.harvard.edu/abs/1995ApJS..100..405F)
;   lines['CIV1168a'] = 1168.873d ; T only
;   lines['CIV1168b'] = 1168.990d ; T only
   lines['CIV1548'] = 1548.202d
   lines['CIV1550'] = 1550.774d

;  S lists one other ClI line. Incorrectly labeled CII in T.
;   lines['ClI1347'] = 1347.2396d ; T/S

;   lines['CuII1359'] = 1358.773d ; T/S

;  Includes most Fe lines from P01 up to lambda = 1600 A
;  Fe II lines are only those with f_jk > 0.01
   lines['FeII1063'] = 1063.170d
   lines['FeII1081'] = 1081.872d
   lines['FeII1096'] = 1096.871d
   lines['FeII1121'] = 1121.967d
   lines['FeII1125'] = 1125.442d
   lines['FeII1143'] = 1143.220d
   lines['FeII1145'] = 1144.926d
   lines['FeII1260'] = 1260.525d
   lines['FeIII1122'] = 1122.53d ; no energy structure listed in ASD

;  Includes all N lines from P01  up to lambda = 1600 A
   lines['NI953a'] = 953.415d ; added by DSNR; f_jk similar to 963/4/5
   lines['NI953b'] = 953.655d ; added by DSNR; f_jk similar to 963/4/5
   lines['NI953c'] = 953.970d ; added by DSNR; f_jk similar to 963/4/5
   lines['NI963'] = 963.990d
   lines['NI964'] = 964.626d
   lines['NI965'] = 965.041d
   lines['NI1134a'] = 1134.165d
   lines['NI1134b'] = 1134.415d
   lines['NI1134c'] = 1134.980d
   lines['NI1199'] = 1199.550d
   lines['NI1200a'] = 1200.223d
   lines['NI1200b']  = 1200.710d
   lines['NII915'] = 915.612d
   lines['NII1083'] = 1083.990d
   lines['NIII989'] = 989.790d
;   lines['NIIIs1183'] = 1183.031d ; T only; f_jk similar to NIII989
;   lines['NIIIs1184'] = 1184.550d ; T only; f_jk similar to NIII989
;   lines['NIVs1168'] = 1168.599d; T only; not sure of line ID
   lines['NV1238'] = 1238.821d
   lines['NV1242'] = 1242.804d

;  Includes all Ni lines from P01 up to lambda = 1600 A
;  No energy level info. on NI II lines in ASD
   lines['NiII1317'] = 1317.22d
   lines['NiII1370'] = 1370.14d
   lines['NiII1454'] = 1454.8420d ; lambda from P01 (not in ASD)
   lines['NiII1467a'] = 1467.2590d ; lambda from P01 (not in ASD)
   lines['NiII1467b'] = 1467.7560d ; lambda from P01 (not in ASD)

;  Includes all O lines from P01 up to lambda = 1600 A
;  Except OI 921, 922, 925, 930, which don't appear in ASD
   lines['OI924'] = 924.950d
   lines['OI929'] = 929.517d
   lines['OI936'] = 936.629d
   lines['OI948'] = 948.646d
   lines['OI950'] = 950.885d
   lines['OI971'] = 971.738d
   lines['OI976'] = 976.448d
   lines['OI988'] = 988.773d ; 2 other v. nearby lines, but 5-50x weaker
   lines['OI1025'] = 1025.762d
   lines['OI1039'] = 1039.234d
   lines['OI1302'] = 1302.168d
   lines['OI1355'] = 1355.598d
   lines['OVI1031'] = 1031.912d
   lines['OVI1037'] = 1037.613d
   
;  Includes all P lines from P01 up to lambda = 1600 A
   lines['PII963'] = 963.8010d ; lambda from P01 (not in ASD)
   lines['PII1152'] = 1152.818d ; lambda from P01 (not in ASD)
   lines['PII1301'] = 1301.87d
   lines['PIII917'] = 917.120d
   lines['PIV950'] = 950.655d ; added by DSNR; high f_jk
   lines['PV1117'] = 1117.98d
   lines['PV1128'] = 1128.01d
   
;  Includes all S lines from P01 up to lambda = 1600 A
   lines['SII1250'] = 1250.578d
   lines['SII1253'] = 1253.805d
   lines['SII1259'] = 1259.518d
   lines['SIII1012'] = 1012.504
   lines['SIII1190'] = 1190.206d
   lines['SIV1062'] = 1062.656
   lines['SVI933'] = 933.376d
   lines['SVI944'] = 944.525d
   
;  Includes all Si lines from P01 up to lambda = 1600 A
   lines['SiII989'] = 989.87d
   lines['SiII1020'] = 1020.70d
   lines['SiII1190'] = 1190.42d
   lines['SiII1193'] = 1193.28d
   lines['SiII1260'] = 1260.42d
   lines['SiII1304'] = 1304.37d
   lines['SiII1526'] = 1526.72d
   lines['SiIII1206'] = 1206.51d
   lines['SiIV1393'] = 1393.76d
   lines['SiIV1402'] = 1402.77d
   
;   lines['SnII1400'] = 1400.454d ; T only
 
   if keyword_set(linelab) then begin
      linelab = hash()
;     Optical lines
      linelab['Halpha'] = 'H$\alpha$'
      linelab['Hbeta'] = 'H$\beta$'
      linelab['Hgamma'] = 'H$\gamma$'
      linelab['HeI5876'] = 'HeI 5876'
      linelab['HeI6678'] = 'HeI 6678'
      linelab['HeI7065'] = 'HeI 7065'
      linelab['HeII4686'] = 'HeII 4686'
      linelab['[NeIII]3869'] = '[NeIII] 3869'
      linelab['[NI]5198'] = '[NI] 5198'
      linelab['[NI]5200'] = '[NI] 5200'
      linelab['[NII]5755'] = '[NII] 5755'
      linelab['[NII]6548'] = '[NII] 6548'
      linelab['[NII]6583'] = '[NII] 6583'
      linelab['[OI]5577'] = '[OI] 5577'
      linelab['[OI]6300'] = '[OI] 6300'
      linelab['[OI]6364'] = '[OI] 6364'
      linelab['[OII]3726'] = '[OII] 3726'
      linelab['[OII]3729'] = '[OII] 3729'
      linelab['[OIII]4959'] = '[OIII] 4959'
      linelab['[OIII]5007'] = '[OIII] 5007'
      linelab['[SII]6716'] = '[SII] 6716'
      linelab['[SII]6731'] = '[SII] 6731'
      linelab['[SIII]6312'] = '[SIII] 6312'
      linelab['Mg1b5167'] = 'MgI 5167'
      linelab['Mg1b5173'] = 'MgI 5173'
      linelab['Mg1b5184'] = 'MgI 5184'
      linelab['NaD2'] = 'NaI 5890'
      linelab['NaD1'] = 'NaI 5896'
      linelab['OH8344'] = ''
      linelab['OH8399'] = ''
      linelab['OH8430'] = ''
      linelab['[CaV]5309'] = '[CaV] 5309'
      linelab['[FeVII]5159'] = '[FeVII] 5159'
      linelab['[FeVII]5276'] = '[FeVII] 5276'
      linelab['[FeVII]5721'] = '[FeVII] 5721'
      linelab['[FeVII]6087'] = '[FeVII] 6087'
      linelab['[FeX]6375'] = '[FeX] 6375'

;     UV lines

      linelab['Lyalpha'] = 'Ly$\alpha$'
      linelab['Lybeta'] = 'Ly$\beta$'
      linelab['Lygamma'] = 'Ly$\gamma$'
      linelab['Lydelta'] = 'Ly$\delta$'
      linelab['Lyepsilon'] = 'Ly$\epsilon$'
      linelab['Ly6'] = 'Ly6'
      linelab['Ly7'] = 'Ly7'
      linelab['Ly8'] = 'Ly8'
      linelab['Ly9'] = 'Ly9'
      linelab['Ly10'] = 'Ly10'
      linelab['Ly11'] = 'Ly11'

      linelab['CI1139'] ='CI 1139'
      linelab['CI1157'] ='CI 1157'
      linelab['CI1277'] ='CI 1277'
      linelab['CI1280'] ='CI 1280'
      linelab['CI1328'] = 'CI 1328'
      linelab['CI1560'] = 'CI 1560'
      linelab['CIs1560'] = 'CI* 1560'
      linelab['CII1036'] = 'CII 1036'
      linelab['CII1334'] = 'CII 1334'
      linelab['CIIs1335'] ='CII* 1335'
      linelab['CIII977'] = 'CIII 977'
      linelab['CIV1548'] = 'CIV 1548'
      linelab['CIV1550'] = 'CIV 1550'

      linelab['FeIII1122'] = 'FeIII 1122'
      linelab['FeII1063'] = 'FeII 1063'
      linelab['FeII1081'] = 'FeII 1081'
      linelab['FeII1096'] = 'FeII 1096'
      linelab['FeII1121'] = 'FeII 1121'
      linelab['FeII1125'] = 'FeII 1125'
      linelab['FeII1143'] = 'FeII 1143'
      linelab['FeII1145'] = 'FeII 1145'
      linelab['FeII1260'] = 'FeII 1260'

      linelab['NI953a'] = 'NI 953'
      linelab['NI953b'] = 'NI 953'
      linelab['NI953c'] = 'NI 953'
      linelab['NI963'] = 'NI 963'
      linelab['NI964'] = 'NI 964'
      linelab['NI965'] = 'NI 965'
      linelab['NI1134a'] = 'NI 1134'
      linelab['NI1134b'] = 'NI 1134'
      linelab['NI1134c'] = 'NI 1134'
      linelab['NI1199'] = 'NI 1199'
      linelab['NI1200a'] = 'NI 1200'
      linelab['NI1200b']  = 'NI 1200'
      linelab['NII915'] = 'NII 915'
      linelab['NII1083'] = 'NII 1083'
      linelab['NIII989'] = 'NIII 989'
      linelab['NV1238'] = 'NV 1238'
      linelab['NV1242'] = 'NV 1242'

      linelab['NiII1317'] ='NiII 1317'
      linelab['NiII1370'] = 'NiII 1370'
      linelab['NiII1454'] = 'NiII 1454'
      linelab['NiII1467a'] = 'NiII 1467'
      linelab['NiII1467b'] = 'NiII 1467'

      linelab['OI924'] = 'OI 924'
      linelab['OI929'] = 'OI 929'
      linelab['OI936'] = 'OI 936'
      linelab['OI948'] = 'OI 948'
      linelab['OI950'] = 'OI 950'
      linelab['OI971'] = 'OI 971'
      linelab['OI976'] = 'OI 976'
      linelab['OI988'] = 'OI 988'
      linelab['OI1025'] = 'OI 1025'
      linelab['OI1039'] = 'OI 1039'
      linelab['OI1302'] = 'OI 1302'
      linelab['OI1355'] = 'OI 1355'
      linelab['OVI1031'] = 'OVI 1031'
      linelab['OVI1037'] = 'OVI 1037'

      linelab['PII963'] = 'PII 963'
      linelab['PII1152'] = 'PII 1152'
      linelab['PII1301'] = 'PII 1301'
      linelab['PIII917'] = 'PIII 917'
      linelab['PIV950'] = 'PIV 950'
      linelab['PV1117'] = 'PV 1117'
      linelab['PV1128'] = 'PV 1128'

      linelab['SII1250'] = 'SII 1250'
      linelab['SII1253'] ='SII 1253'
      linelab['SII1259'] = 'SII 1259'
      linelab['SIII1012'] = 'SIII 1012'
      linelab['SIII1190'] = 'SIII 1190'
      linelab['SIV1062'] = 'SIV 1062'
      linelab['SVI933'] = 'SVI 933'
      linelab['SVI944'] = 'SVI 944'

      linelab['SiII989'] ='SiII 989'
      linelab['SiII1020'] ='SiII 1020'
      linelab['SiII1190'] ='SiII 1190'
      linelab['SiII1193'] = 'SiII 1193'
      linelab['SiII1260'] = 'SiII 1260'
      linelab['SiII1304'] ='SiII 1304'
      linelab['SiII1526'] ='SiII 1526'
      linelab['SiIII1206'] = 'SiIII 1206'
      linelab['SiIV1393'] = 'SiIV 1393'
      linelab['SiIV1402'] = 'SiIV 1402'

;      lines['SnII1400'] = 'SnII 1400'

    endif

   if keyword_set(all) then outlines = lines $
   else begin
      outlines = hash()
      for i=0, n_elements(inlines)-1 do begin
         imatch = where(inlines[i] eq lines->keys(),ctmatch)
         if ctmatch eq 1 then outlines[inlines[i]] = lines[inlines[i]] $
         else print,'IFSF_LINELIST: ERROR: ',inlines[i],$
            ' not found in wavelength list.'
       endfor
   endelse
 
   return,outlines
 
end
