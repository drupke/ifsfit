; docformat = 'rst'
;
;+
;
; A lookup table for wavelengths.
;
; Wavelength sources:
;   1. NIST ASD -- For UV, observed if present, otherwise Ritz
;   2. HeI 4686 -- found random wavelength compilation. In NIST, shows
;                  up as a complex of many fine structure lines with
;                  range of several tenths of an A.
;   3. Some other random sources -- see inline comments.
;   
; For IR lines:
;   1. H_2: http://www.jach.hawaii.edu/UKIRT/astronomy/calib/spec_cal/h2_s.html
;   2. He: http://www.jach.hawaii.edu/UKIRT/astronomy/calib/spec_cal/he_4.html
;   3. H: http://www.jach.hawaii.edu/UKIRT/astronomy/calib/spec_cal/recomb_menu.html
;   4. [SiVI]: NIST
;   Checked He and H lines in NIST and there is good agreement.
;
; To get a hash of all lines with labels:
; 
; IDL> linelab = 1b
; IDL> lines = ifsf_linelist(/all,linelab=linelab)
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Hash of wavelengths for specified input lines, in Angstroms.
;
; :Params:
;    inlines: in, required, type=strarr
;      Line labels for which to retrieve wavelengths.
;
; :Keywords:
;    all: in, optional, type=byte
;      Return all lines.
;    linelab: in, optional, type=hash
;      Return line labels.
;    waveunit: in, optional, type=double
;      If set, multiply each wavelength by this number. Default is to Angstroms.
;    quiet: in, optional, type=byte
;      Don't output anything to STDOUT.
;    vacuum: in, optional, type=byte
;      Default is to air wavelengths. If set, output vacuum wavelengths from a
;      corresponding lookup table.
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
;      2016sep12, DSNR, added IR lines
;      2020may26, DSNR, improved some documentation, added VACUUM option
;      2020jun05, DSNR, added [NeV], [NeIII], MgII doublet cases;
;                       now check lines/labels only on output lines
;      2020jun08, DSNR, now returns only line labels corresponding to output lines
;    
; :Copyright:
;    Copyright (C) 2013--2020 David S. N. Rupke, Anthony To
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
function ifsf_linelist,inlines,linelab=linelab,all=all,waveunit=waveunit,$
                       vacuum=vacuum,quiet=quiet

;  Associated line labels:
   lines = hash()

   if ~ keyword_set(quiet) then $
      print,'NOTE: Calling sequence to IFSF_LINELIST has changed. If vacuum',$
            '      wavelengths are required, or other wavelength units, or ',$
            '      to suppress this message, use ARGSLINELIST structure in ',$
            '      initialization procedure. Default is now air waves.',$
            format='(A0/A0/A0/A0)'

   if keyword_set(vacuum) then begin

;  IR lines (vacuum)
      lines['Paa'] = 18756d
      lines['Pab'] = 12821.578d
      lines['Pag'] = 10941.17d
      lines['Pad'] = 10052.6d
      lines['Brg'] = 21661d
      lines['Brd'] = 19451d
      lines['Bre'] = 18181d
      lines['H2_10_S1'] = 21218d
      lines['H2_10_S2'] = 20338d
      lines['H2_10_S3'] = 19575d
      lines['H2_10_S4'] = 18920d
      lines['H2_10_S5'] = 18358d
      lines['HeI206'] = 20587d
      lines['HeI187'] = 18691d
      lines['HeI108'] = 10833d
      lines['HeI108a'] = 10832.057472d
      lines['HeI108b'] = 10833.216751d
      lines['HeI108c'] = 10833.306444d
      lines['[ArIII]7138'] = 7137.77d
      lines['[ArIII]7753'] = 7753.19d
      lines['[SIII]9071'] = 9071.1d
      lines['[SIII]9533'] = 9533.2d   
      lines['[SII]103'] = 10329d
      lines['[SII]103a'] = 10289.55d
      lines['[SII]103b'] = 10323.32d
      lines['[SII]103c'] = 10329.24d
      lines['[SII]103d'] = 10373.33d
      lines['[SiVI]196'] = 19646d
 
      
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
      lines['PV1117'] = 1117.9774d
      lines['PV1128'] = 1128.0078d

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

;  Lines seen in
;     Vestergaard & Wilkes 01
;     Veilleux 02 ASP Conf. Ser.
;     Vanden Berk+01
;  Wavelengths in vacuum from NIST ASD

;  add HeI lines ...
      lines['HeII1640'] = 1640.4d ; several lines
; (up to 7? from 1640.33 to 1640.53), according to NIST ASD
      lines['HeII2733'] = 2734.1d ; several lines ; in vacuum
;   lines['HeII2734'] = 2733.3d ; several lines ; in air
; (up to 13? from 2733.99 to 2734.16), according to NIST ASD

;  check whether allowed or intercombination or forbidden
      lines['NIII1747'] = 1747.848d
      lines['NIII1751a'] = 1751.218d
      lines['NIII1751b'] = 1751.657d

;  add CIII] + CII] lines ...
      lines['[CIII]1906'] = 1906.683d
      lines['CIII]1908'] = 1908.734d
      lines['CII]2326'] = 2326.11d
      lines['CII]2329'] = 2328.84d

;  add [NeIV] 1601 lines ...
      lines['[NeIV]2422'] = 2422.561d ; [M1] or [E2]
      lines['[NeIV]2425'] = 2425.139d ; [M1+E2]

      lines['[SiIII]1883'] = 1882.7d ; from Laor+97
      lines['SiIII]1892'] = 1892.03d ; from http://astronomy.nmsu.edu/drewski/tableofemissionlines.html

;  check wavelengths
      lines['[OII]2470'] = 2470.966d
      lines['[OII]2471'] = 2471.088d

;  add [OIII] lines

      lines['[FeIV]2568a'] = 2568.17d ; A = 3.8e-2 [M1]
      lines['[FeIV]2568b'] = 2568.38d ; A = 5.1e-2 [M1]
      lines['[FeIV]2830'] = 2830.19d ; A = 8.8e-1 [M1]
      lines['[FeIV]2836'] = 2836.57 ; A = 1.4 [M1]
      lines['[FeV]2708'] = 2708.2d ; A = 2.5e-1 [M1]
      lines['[FeVI]2145'] = 2145.75d ; A = 2.2e-1 [M1]; 4.5e-3 [E2]
      lines['[FeVI]2163'] = 2163.69d ; A = 6.6e-4 [E2]
      lines['[FeXI]2649'] = 2649.50d ; A = 92 [M1] ; 1.5e-1 [E2]
      lines['[FeXII]2405'] = 2405.8d ; A = 48 [M1]; 4.0e-2 [E2]

;  check wavelengths
      lines['[MgVI]1805'] = 1805.96d ; a guess!
      lines['[MgVII]2509'] = 2509.2d ; a guess! lower level close to but not at ground state
      lines['[MgVII]2629'] = 2629.1d ; a guess! lower level close to but not at ground state
      lines['[MgV]2782'] = 2782.7d
      lines['[MgV]2928'] = 2928.0d ; a guess! lower level close to but not at ground state

;  check wavelengths
      lines['[ArIV]2853'] = 2853.67d
      lines['[ArIV]2868'] = 2868.21d

;  check wavelengths
      lines['[CaVI]2214'] = 2214.51d
      lines['[CaVI]2242'] = 2242.13d

;  add NaV lines ...

;     Resonant lines
      lines['MgII2796'] = 2796.352 ; f = 0.6123
      lines['MgII2803'] = 2803.530 ; f = 0.3054
      lines['MgII2796+MgII2803'] = (2796.352d + 2803.530d)/2d ; f = 0.3054
      lines['MgI2852'] = 2852.965 ; f = 1.81
;;  UV5
;      lines['FeII2249'] = 2249.1795d
;;  UV4
;      lines['FeII2260'] = 2260.0809d
;      lines['FeII2233'] = 2233.7532d
;;  UV3
;      lines['FeII2343'] = 2343.4948d ; f = 0.1140 (~2 times smaller than 2599 line)
;;  UV2
;      lines['FeII2366'] = 2366.8674d
;      lines['FeII2373'] = 2373.73528d; f = 0.0313 (~8 times smaller than 2599 line)
;      lines['FeII2382'] = 2382.03733d; f = 0.320 (larger than 2599 line)
;;  UV1
;      lines['FeII2585'] = 2585.87560d ; f = 0.0691 (~3.5 times smaller than 2599 line)
;      lines['FeII2599'] = 2599.39515d ; f = 0.239
;;  non-resonant lines from Finley+17
;;  UV3
;      lines['FeII*2364'] = 2364.8278d ; goes with 2343
;      lines['FeII*2380'] = 2380.76131d ; goes with 2343
;;  UV2
;      lines['FeII*2395'] = 2395.62504d ; goes with 2373
;;  UV1
;      lines['FeII*2611'] = 2611.87336d ; goes with 2585
;      lines['FeII*2631'] = 2631.32292d ; goes with 2585
;      lines['FeII*2625'] = 2625.66685d ; goes with 2599

   endif else begin
   
;  Optical lines
      lines['Halpha'] = 6562.80d
      lines['Hbeta'] = 4861.35d
      lines['Hgamma'] = 4340.47d
      lines['Hdelta'] = 4101.73d
      lines['Hepsilon'] = 3970.075d
      lines['H8'] = 3889.064d
      lines['H9'] = 3835.397d
      lines['HeI5876'] = 5875.661d
      lines['HeI6678'] = 6678.15d
      lines['HeI7065'] = 7065.19d
      lines['HeII2733'] = 2733.3d ; a little uncertain; ~10 lines that could contribute
      lines['HeII3203'] = 3203.1d ; a little uncertain; ~10 lines that could contribute
      lines['HeII4686'] = 4686.7d
      lines['[NeIII]3869'] = 3868.76d
      lines['[NeIII]3967'] = 3967.47d
      lines['[NeIII]3869+[NeIII]3967'] = (3868.76d + 3967.47d)/2d
      lines['[NI]5198'] = 5197.90d
      lines['[NI]5200'] = 5200.26d
      lines['[NI]5198+[NI]5200'] = (5197.90d + 5200.26d)/2d
      lines['[NII]5755'] = 5754.59
      lines['[NII]6548'] = 6548.05d
      lines['[NII]6583'] = 6583.45d
      lines['[OI]5577'] = 5577.34d
      lines['[OI]6300'] = 6300.30d
      lines['[OI]6364'] = 6363.78d
      lines['[OII]3726'] = 3726.032d
      lines['[OII]3729'] = 3728.815d
      lines['[OII]3726+[OII]3729'] = (3726.032d + 3728.815d)/2d
      lines['[OIII]4363'] = 4363.209d
      lines['[OIII]4959'] = 4958.91d
      lines['[OIII]5007'] = 5006.84d
      lines['[SII]4068'] = 4068.60d
      lines['[SII]4076'] = 4076.35d
      lines['[SII]6716'] = 6716.44d
      lines['[SII]6731'] = 6730.82d
      lines['[SII]6716+[SII]6731'] = (6716.44d + 6730.82d)/2d
      lines['[SIII]6312'] = 6312.06d
      lines['Mg1b5167'] = 5167.3213
      lines['Mg1b5173'] = 5172.6844
      lines['Mg1b5184'] = 5183.6043  
      lines['NaD2'] = 5889.95d
      lines['NaD1'] = 5895.92d
      lines['OH6287a'] = 6287.407 ; 9-3 P1e(2.5), Osterbrock 1996
      lines['OH6287b'] = 6287.462 ; 9-3 P1f(2.5), Osterbrock 1996
      lines['OH6306a'] = 6306.869 ; 9-3 P1e(3.5)
      lines['OH6306b'] = 6306.981 ; 9-3 P1f(3.5)
      lines['OH6329a'] = 6329.747 ; 9-3 P1e(4.5)
      lines['OH6329b'] = 6329.933 ; 9-3 P1f(4.5)
      lines['OH6356a'] = 6356.167 ; 9-3 P1e(5.5)
      lines['OH6356b'] = 6356.441 ; 9-3 P1f(5.5)
      lines['OH8344'] = 8344.602d
      lines['OH8399'] = 8399.160d
      lines['OH8430'] = 8430.170d
      lines['[ArIV]4740'] = 4740.12d
      lines['[CaV]5309'] = 5309.11d
      lines['[FeVII]5159'] = 5158.89d
      lines['[FeVII]5276'] = 5276.38d
      lines['[FeVII]5721'] = 5720.7d
      lines['[FeVII]6087'] = 6087.0d
      lines['[FeX]6375'] = 6374.51d

;  resonant lines
;  air wavelengths
;  from P01 -- wavelengths from NIST ASD, f values from P01
      lines['MgI2025'] = 2025.824d ; f = 0.112 (~16 times smaller than 2852 line)
      lines['MgII2796'] = 2795.528d ; f = 0.6123
      lines['MgII2803'] = 2802.704d ; f = 0.3054
      lines['MgII2796+MgII2803'] = (2795.528d + 2802.704d)/2d ; f = 0.3054
      lines['MgI2852'] = 2852.127d ; f = 1.81
;  UV5
      lines['FeII2249'] = 2249.1795d ; from Morton 2003
;  UV4
      lines['FeII2260'] = 2260.0809d ; from Morton 2003
      lines['FeII2233'] = 2233.7532d ; from Morton 2003
;  UV3
      lines['FeII2343'] = 2343.4948d ; f = 0.1140 (~2 times smaller than 2599 line)
;  UV2
      lines['FeII2366'] = 2366.8674d ; from Morton 2003
      lines['FeII2373'] = 2373.73528d; f = 0.0313 (~8 times smaller than 2599 line)
      lines['FeII2382'] = 2382.03733d; f = 0.320 (larger than 2599 line)
;  UV1
      lines['FeII2585'] = 2585.87560d ; f = 0.0691 (~3.5 times smaller than 2599 line)
      lines['FeII2599'] = 2599.39515d ; f = 0.239
;  non-resonant lines from Finley+17
;  UV3
      lines['FeII*2364'] = 2364.8278d ; goes with 2343
      lines['FeII*2380'] = 2380.76131d ; goes with 2343
;  UV2
      lines['FeII*2395'] = 2395.62504d ; goes with 2373
;  UV1
      lines['FeII*2611'] = 2611.87336d ; goes with 2585
      lines['FeII*2631'] = 2631.32292d ; goes with 2585
      lines['FeII*2625'] = 2625.66685d ; goes with 2599
; From Morton 2003  
      lines['MnII2605'] = 2605.684d
      lines['MnII2593'] = 2593.724d
      lines['MnII2576'] = 2576.105d

;  Al lines from P01 ... check wavelengths and what type of line
;  add Al II at 2660/2669
      lines['AlII1670'] = 1670.7874d
      lines['AlIII1854'] = 1854.7164d
      lines['AlIII1862'] = 1862.7895d

;  other resonant lines and connected non-resonant lines
;  see Laha, Keenan, Ferland+16a,b for info. on SiII multiplet
      lines['SiII1808']  = 1808.00d ; A = 2.54e6
      lines['SiII*1816'] = 1816.92d ; A = 2.65e6
      lines['SiII*1817'] = 1817.45d ; A = 3.23e5
      lines['SiII2329']  = 2329.23d ; A = 2.35e1
      lines['SiII2335']  = 2335.12d ; A = 5.51e3
      lines['SiII*2335'] = 2335.33d ; A = 2.44e3
      lines['SiII*2344'] = 2344.92d ; A = 1.31e3
      lines['SiII*2350'] = 2350.89d ; A = 4.70e3

; Wavelengths from NIST ASD
      lines['[NeIV]2422'] = 2421.825d ; [M1] or [E2] resonant
      lines['[NeIV]2425'] = 2424.403d ; [M1+E2] resonant
      lines['[NeV]3345'] = 3345.83d ; non-resonant
      lines['[NeV]3426'] = 3425.87d ; non-resonant
      lines['[NeV]3345+[NeV]3426'] = (3345.83d + 3425.87d)/2d ; non-resonant

   endelse

   if keyword_set(linelab) then begin
      linelab = hash()

;     IR lines
      linelab['Paa'] = 'P$\alpha$'
      linelab['Pab'] = 'P$\beta$'
      linelab['Pag'] = 'P$\gamma$'
      linelab['Pad'] = 'P$\delta$'
      linelab['Brg'] = 'Br$\gamma$'
      linelab['Brd'] = 'Br$\delta$'
      linelab['Bre'] = 'Br$\epsilon$'
      linelab['H2_10_S1'] = 'H$\down2$(1-0) S(1)'
      linelab['H2_10_S2'] = 'H$\down2$(1-0) S(2)'
      linelab['H2_10_S3'] = 'H$\down2$(1-0) S(3)'
      linelab['H2_10_S4'] = 'H$\down2$(1-0) S(4)'
      linelab['H2_10_S5'] = 'H$\down2$(1-0) S(5)'
      linelab['HeI206'] = 'HeI 2.06$\mu$m'
      linelab['HeI187'] = 'HeI 1.87$\mu$m'
      linelab['HeI108'] = 'HeI 1.08$\mu$m'
      linelab['HeI108a'] = 'HeI 1.08320$\mu$m'
      linelab['HeI108b'] = 'HeI 1.08332$\mu$m'
      linelab['HeI108c'] = 'HeI 1.08333$\mu$m'
      linelab['[ArIII]7138'] = '[ArIII] 7138'
      linelab['[ArIII]7753'] = '[ArIII] 7753'
      linelab['[SIII]9071'] = '[SIII] 9071'
      linelab['[SIII]9533'] = '[SIII] 9533'
      linelab['[SII]103'] = '[SII] 1.03$\mu$m'
      linelab['[SII]103a'] = '[SII] 1.0289$\mu$m'
      linelab['[SII]103b'] = '[SII] 1.0323$\mu$m'
      linelab['[SII]103c'] = '[SII] 1.0329$\mu$m'
      linelab['[SII]103d'] = '[SII] 1.0373$\mu$m'
      linelab['[SiVI]196'] = '[SiVI] 1.96$\mu$m'

;     Optical lines
      linelab['Halpha'] = 'H$\alpha$'
      linelab['Hbeta'] = 'H$\beta$'
      linelab['Hgamma'] = 'H$\gamma$'
      linelab['Hdelta'] = 'H$\delta$'
      linelab['Hepsilon'] = 'H$\epsilon$'
      linelab['H8'] = 'H8'
      linelab['H9'] = 'H9'
      linelab['HeI5876'] = 'HeI 5876'
      linelab['HeI6678'] = 'HeI 6678'
      linelab['HeI7065'] = 'HeI 7065'
      linelab['HeII1640'] = 'HeII 1640'
      linelab['HeII2733'] = 'HeII 2733'
      linelab['HeII3203'] = 'HeII 3203'
      linelab['HeII4686'] = 'HeII 4686'

      linelab['[NI]5198'] = '[NI] 5198'
      linelab['[NI]5200'] = '[NI] 5200'
      linelab['[NI]5198+[NI]5200'] = '[NI] 5198, 5200'
      linelab['[NII]5755'] = '[NII] 5755'
      linelab['[NII]6548'] = '[NII] 6548'
      linelab['[NII]6583'] = '[NII] 6583'
      linelab['[OI]5577'] = '[OI] 5577'
      linelab['[OI]6300'] = '[OI] 6300'
      linelab['[OI]6364'] = '[OI] 6364'
      linelab['[OII]3726'] = '[OII] 3726'
      linelab['[OII]3729'] = '[OII] 3729'
      linelab['[OII]3726+[OII]3729'] = '[OII] 3726, 3729'
      linelab['[OIII]4363'] = '[OIII] 4363'
      linelab['[OIII]4959'] = '[OIII] 4959'
      linelab['[OIII]5007'] = '[OIII] 5007'
      linelab['[SII]4068'] = '[SII] 4068'
      linelab['[SII]4076'] = '[SII] 4076'
      linelab['[SII]6716'] = '[SII] 6716'
      linelab['[SII]6731'] = '[SII] 6731'
      linelab['[SII]6716+[SII]6731'] = '[SII] 6716, 6731'
      linelab['[SIII]6312'] = '[SIII] 6312'
      linelab['Mg1b5167'] = 'MgI 5167'
      linelab['Mg1b5173'] = 'MgI 5173'
      linelab['Mg1b5184'] = 'MgI 5184'
      linelab['NaD2'] = 'NaI 5890'
      linelab['NaD1'] = 'NaI 5896'
      linelab['OH6287a'] = ''
      linelab['OH6287b'] = ''
      linelab['OH6306a'] = ''
      linelab['OH6306b'] = ''
      linelab['OH6329a'] = ''
      linelab['OH6329b'] = ''
      linelab['OH6356a'] = ''
      linelab['OH6356b'] = ''
      linelab['OH8344'] = ''
      linelab['OH8399'] = ''
      linelab['OH8430'] = ''
      linelab['[ArIV]4740'] = '[ArIV] 4740'
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

      linelab['NIII1747'] = 'NIII] 1747'
      linelab['NIII1751a'] = 'NIII] 1751'
      linelab['NIII1751b'] = 'NIII] 1751'

      linelab['[CIII]1906'] = '[CIII] 1906'
      linelab['CIII]1908'] = 'CIII] 1908'
      linelab['CII]2326'] = 'CII] 2326'
      linelab['CII]2329'] = 'CII] 2329'

      linelab['[SiIII]1883'] = '[SiIII] 1883'
      linelab['SiIII]1892'] = 'SiIII] 1892'

      linelab['[NeIII]3869'] = '[NeIII] 3869'
      linelab['[NeIII]3967'] = '[NeIII] 3967'
      linelab['[NeIII]3869+[NeIII]3967'] = '[NeIII] 3869, 3967'
      linelab['[NeIV]2422'] = '[NeIV] 2422'
      linelab['[NeIV]2425'] = '[NeIV] 2425'
      linelab['[NeV]3345'] = '[NeV] 3345'
      linelab['[NeV]3426'] = '[NeV] 3426'
      linelab['[NeV]3345+[NeV]3426'] = '[NeV] 3345, 3426'

      linelab['[OII]2470'] = '[OII] 2470'
      linelab['[OII]2471'] = '[OII] 2471'

      linelab['[FeIV]2568a'] = '[FeIV] 2568a'
      linelab['[FeIV]2568b'] = '[FeIV] 2568b'
      linelab['[FeIV]2830'] = '[FeIV] 2830'
      linelab['[FeIV]2836'] = '[FeIV] 2836'
      linelab['[FeV]2708'] = '[FeV] 2708'
      linelab['[FeVI]2145'] = '[FeVI] 2145'
      linelab['[FeVI]2163'] = '[FeVI] 2163'
      linelab['[FeXI]2649'] = '[FeXI] 2649'
      linelab['[FeXII]2405'] = '[FeXII] 2405'

      linelab['[MgVI]1805'] = '[MgVI] 1805'
      linelab['[MgVII]2509'] = '[MgVII] 2509'
      linelab['[MgVII]2629'] = '[MgVII] 2629'
      linelab['[MgV]2782'] = '[MgV] 2782'
      linelab['[MgV]2928'] = '[MgV] 2928'

      linelab['[ArIV]2853'] = '[ArIV] 2853'
      linelab['[ArIV]2868'] = '[ArIV] 2868'

      linelab['[CaVI]2214'] = '[CaVI]2214'
      linelab['[CaVI]2242'] = '[CaVI]2242'

      linelab['MgI2025'] = 'MgI 2025'
      linelab['MgII2796'] = 'MgII 2796'
      linelab['MgII2803'] = 'MgII 2803'
      linelab['MgII2796+MgII2803'] = 'MgII 2796, 2803'
      linelab['MgI2852'] = 'MgI 2852'
      linelab['FeII2233'] = 'FeII 2233'
      linelab['FeII2249'] = 'FeII 2249'
      linelab['FeII2260'] = 'FeII 2260'
      linelab['FeII2343'] = 'FeII 2343'
      linelab['FeII2366'] = 'FeII 2366'
      linelab['FeII2373'] = 'FeII 2373'
      linelab['FeII2382'] = 'FeII 2382'
      linelab['FeII2585'] = 'FeII 2585'
      linelab['FeII2599'] = 'FeII 2599'
      linelab['FeII*2364'] = 'FeII* 2364'
      linelab['FeII*2380'] = 'FeII* 2380'
      linelab['FeII*2395'] = 'FeII* 2395'
      linelab['FeII*2611'] = 'FeII* 2611'
      linelab['FeII*2625'] = 'FeII* 2625'
      linelab['FeII*2631'] = 'FeII* 2631'
      linelab['MnII2605'] = 'MnII 2605'
      linelab['MnII2593'] = 'MnII 2593'
      linelab['MnII2576'] = 'MnII 2576'

      linelab['AlII1670'] = 'AlII 1670'
      linelab['AlIII1854'] = 'AlIII 1854'
      linelab['AlIII1862'] = 'AlIII 1862'

      linelab['SiII1808'] = 'SiII 1808'
      linelab['SiII*1816'] = 'SiII* 1816'
      linelab['SiII*1817'] = 'SiII* 1817'
      linelab['SiII2329'] = 'SiII* 2329'
      linelab['SiII2335'] = 'SiII 2335'
      linelab['SiII*2335'] = 'SiII* 2335'
      linelab['SiII*2344'] = 'SiII* 2344'
      linelab['SiII*2350'] = 'SiII* 2350'

    endif
   
   if keyword_set(waveunit) then $
      foreach key,lines.keys() do lines[key] *= waveunit

   if keyword_set(all) then outlines = lines $
   else begin
      outlines = hash()
      outlinelab = hash()
      for i=0, n_elements(inlines)-1 do begin
         imatch = where(inlines[i] eq lines->keys(),ctmatch)
         ctlabmatch = 1
         if keyword_set(linelab) then $
            ilabmatch = where(inlines[i] eq linelab->keys(),ctlabmatch)
         if ctmatch eq 1 AND ctlabmatch eq 1 then begin
            outlines[inlines[i]] = lines[inlines[i]]
            if keyword_set(linelab) then $
               outlinelab[inlines[i]] = linelab[inlines[i]]
         endif else print,'IFSF_LINELIST: ERROR: ',inlines[i],$
            ' not found in wavelength list OR line label missing.'
       endfor
   endelse

; No longer need this, previous logic now does this check. 
;;  Make sure every output line has a label, and vice versa!
;   if keyword_set(linelab) then begin
;      foreach key,outlines.keys() do $
;         if ~ linelab.haskey(key) then $
;            print,'IFSF_LINELIST: WARNING: Line ',key,' has no label.'
;      foreach key,linelab.keys() do $
;         if ~ outlines.haskey(key) then $
;            print,'IFSF_LINELIST: WARNING: Line label ',key,' has no wavelength.'
;   endif
 
   linelab = outlinelab
   return,outlines
 
end
