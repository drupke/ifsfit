; docformat = 'rst'
;
;+
;
; Produces plots showing the normalized continuum with absorption profiles
; as well as the continuum fit.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Postscript file with plots.
;
; :Params:
;    initfile: in, required, type=strarr
;      File that holds a list of galaxies as well as their redshifts.
;    directoryname: in, required, type=strarr
;      Location of the data files.
;    galaxyname: in, required, type=strarr
;      Galaxy being plotted, stored in a directory with the same name.
;
; :Keywords:
;    (AT LEAST ONE KEYWORD IS REQUIRED)
;    NV: in, optional, type=string
;      Set if galaxy has been processed and NV doublets were found.
;    OVI: in, optional, type=string
;      Set if galaxy has been processed and OVI doublets were found.
;
; :Author:
;    Anthony Dinh To::
;      Rhodes College
;      Department of Physics
;      2000 N. Parkway
;      Memphis, TN 38104
;      andto94@gmail.com
;
; :History:
;    ChangeHistory::
;      2016jul06, ADT, created
;
; :Copyright:
;    Copyright (C) 2016 Anthony To
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
PRO ifsf_smallerplots, initfile, directoryname, galaxyname, NV=NV, OVI=OVI


;List of emission/absorption lines and their corresponding wavelengths.

LineWavelength = $
    [977.0201,$
    989.799,$
    1025.7223,$
    1031.912,$
    1037.613,$
    1122.5240,$
    1125.4477,$
    1128.0078,$
    1134.1653,$
    1134.4149,$
    1134.9803,$
    1143.2260,$
    1144.9379,$
    1152.8180,$
    1168.599,$
    1168.873,$
    1168.990,$
    1175.7217,$
    1183.030,$
    1184.544,$
    1190.416,$
    1193.2890,$
    1199.5496,$
    1200.2233,$
    1200.7098,$
    1206.500,$
    1215.67,$
    1238.8210,$
    1242.804,$
    1250.578,$
    1253.8051,$
    1259.518,$
    1260.4221,$
    1277.245,$
    1280.14,$
    1302.1685,$
    1304.3702,$
    1317.21,$
    1328.83,$
    1334.532,$
    1335.6627,$
    1335.7077,$
    1347.2396,$
    1355.5977,$
    1358.77,$
    1370.132,$
    1393.76,$
    1400.450,$
    1402.7729,$
    1548.202,$
    1550.774]

  LineLabel = $
    ['C III $\lambda$977',$
    'N III $\lambda$990',$
    'Ly $\beta$1026',$
    'O VI $\lambda$1032',$
    'O VI $\lambda$1038',$
    'Fe III $\lambda$1123',$
    'Fe II $\lambda$1025',$
    'P V $\lambda$1128',$
    'N I $\lambda$1134.2',$
    'N I $\lambda$1134.4',$
    'N I $\lambda$1134.9',$
    'Fe II $\lambda$1143',$
    'Fe II $\lambda$1145',$
    'P II $\lambda$1153',$
    'N IV $\lambda$1169',$
    'C IV $\lambda$1169',$
    'C IV $\lambda$1169',$
    'C III $\lambda$1175.7',$
    'N III $\lambda$1183',$
    'N III $\lambda$1185',$
    'Si II $\lambda$1190',$
    'Si II $\lambda$1193',$
    'N I $\lambda$1199.5',$
    'N I $\lambda$1200.2',$
    'N I $\lambda$1200.71',$
    'Si III $\lambda$1207',$
    'Ly $\alpha$1216',$
    'N V $\lambda$1239',$
    'N V $\lambda$1243',$
    'S II $\lambda$1251',$
    'S II $\lambda$1254',$
    'S II $\lambda$1260', $
    'Si II $\lambda$1260',$
    'C I $\lambda$1277',$
    'C I $\lambda$1280',$
    'O I $\lambda$1302',$
    'Si II $\lambda$1304',$
    'Ni II $\lambda$1317',$
    'C I $\lambda$1329',$
    'C II $\lambda$1335',$
    'C II* $\lambda$1335.6',$
    'C II* $\lambda$1335.7',$
    'C II $\lambda$1347',$
    'O I $\lambda$1356',$
    'Cu II $\lambda$1359',$
    'Ni II $\lambda$1370',$
    'Si IV $\lambda$1394',$
    'Sn II $\lambda$1400',$
    'Si IV $\lambda$1403',$
    'C IV $\lambda$1548',$
    'C IV $\lambda$1551']


;Picks out the redshift of the specified galaxy.
readcol, directoryname+'/'+initfile, galaxynamelist,redshiftlist, SKIPLINE=1,format = '(A,D)' 
selectionparameter=WHERE(galaxynamelist eq galaxyname)
galaxyname=galaxynamelist[selectionparameter[0]]
zsys=redshiftlist[selectionparameter[0]]


;Shifts the emission/absorption wavelengths by the galaxy's redshift
ShiftedLines= LineWavelength*(1+zsys)

;Opens a Postscript file to draw plots on.
CGPS_OPEN,directoryname+'/'+galaxyname+'/'+galaxyname+'smallplot',/ENCAPSULATED
!P.Background='WHITE'
!P.CharThick=.8

print, 'Plotting...'    


;O VI only plots
 
;Read params
IF (keyword_set(OVI)) AND (~ keyword_set(NV)) THEN BEGIN
  readcol, directoryname+'/'+galaxyname+'/'+galaxyname+'OVI_fit_data.txt', wave, modflux, continuum, flux, normalizedflux
  readcol, directoryname+'/'+galaxyname+'/'+galaxyname+'OVI_fit_dataparam.txt', xran_1,xran_2
  readcol, directoryname+'/'+galaxyname+'/'+galaxyname+'OVI_fit_datamodabs.txt', nuvabs,elementsize,NUMLINE=1

;Initializing and fixing variables  
  nuvabs=nuvabs[0]
  unity=make_array(elementsize,Value=1)
  elementsize=elementsize[0]
  moduvabs=Make_array(nuvabs,elementsize)

;Rescaling flux values  
  continuum=continuum/1E-14
  flux = flux/1E-14
  moduvabs=moduvabs/1E-14

;Defining plot ranges  
  xran=[xran_1,xran_2]
  plotregion=[value_locate(wave,xran[0]):value_locate(wave,xran[1])]
  
;Y-range for the normalized spectra plot
  yrannormalized=[0,1.05*Max(normalizedflux[plotregion])]
  
;Y-range for the continuum plot
  yran=[0,1.05*Max(flux[plotregion])]
  
;Normalized plot and legend
  cgplot,wave,normalizedflux,xran=xran,yran=yrannormalized,xstyle=1,ystyle=1,$
    backg='White',axiscolor='Black',color='Black',$
    xtit='Wavelength ($\Angstrom$)',ytit='Normalized F!I$\lambda$',$
    position= [.1,.5,.5,.9], CHARSIZE=.6,thick=1, Title= 'O VI'
  AL_LEGEND,['Intrinsic','ISM','Continuum','Components','Continuum w/Absorption'], $
    Color=['Red','Blue','Orange','Sky Blue','Purple'], charsize=.3, $
    Linestyle=[0,0,0,0,0], Position = [.12,.98],/norm,box=0
    
;Plots absorption features, unity, and outlines spectra
  FOR M=0, nuvabs-1 DO BEGIN
    readcol, directoryname+'/'+galaxyname+'/'+galaxyname+'OVI_fit_datamodabs.txt',moduvabs,skipline=1+elementsize*M, NUMLINE=elementsize
    cgoplot,wave,moduvabs,color='Sky Blue',thick=4
  ENDFOR
  cgoplot,wave,unity,color='Orange',thick=4
  cgoplot,wave, modflux, color = 'Purple', thick = 4
  
;Plots absorption and emission labels
  FOR M = 0, N_ELEMENTS(LineWavelength)-1 DO BEGIN
    cgoplot, [LineWavelength[M],LineWavelength[M]], yrannormalized, color = 'Blue', thick=4
    cgoplot, [ShiftedLines[M],ShiftedLines[M]], yrannormalized, color = 'Red', thick=4
    index=Where(LineWavelength lt Max(wave) AND LineWavelength gt Min(wave))
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,LineWavelength[I], $
        .4*yrannormalized[1],  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .4
    ENDFOR
    index=Where(ShiftedLines lt Max(wave) AND ShiftedLines gt Min(wave))
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,ShiftedLines[I], $
        .4*yrannormalized[1],  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .4
    ENDFOR
  ENDFOR
  
;Plots continuum 
  cgplot,wave,flux,xran=xran,yran=yran,xstyle=1,ystyle=1,$
         backg='White',axiscolor='Black',color='Black',$
         xtickformat="(A1)",ytit='F!I$\lambda$!N/10!E-14!N (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)',$
         position= [.1,.03,.5,.43], /NoErase, CHARSIZE=.6
  cgoplot,wave, continuum, color = 'Orange', thick=4
  
;Plots absorption and emission labels
  FOR M = 0, N_ELEMENTS(LineWavelength)-1 DO BEGIN
    cgoplot, [LineWavelength[M],LineWavelength[M]], yran, color = 'Blue', thick=4
    cgoplot, [ShiftedLines[M],ShiftedLines[M]], yran, color = 'Red', thick=4
    index=Where(LineWavelength lt Max(wave) AND LineWavelength gt Min(wave))
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,LineWavelength[I], $
        .4*yran[1],  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .4
    ENDFOR
    index=Where(ShiftedLines lt Max(wave) AND ShiftedLines gt Min(wave))
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,ShiftedLines[I], $
        .4*yran[1],  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .4, $
        _REF_EXTRA=extra
    ENDFOR
  ENDFOR
  
ENDIF 


;N V only plots

;Read params
IF (keyword_set(NV)) AND (~ keyword_set(OVI)) THEN BEGIN
  readcol, directoryname+'/'+galaxyname+'/'+galaxyname+'NV_fit_data.txt', wave, modflux, continuum, flux, normalizedflux
  readcol, directoryname+'/'+galaxyname+'/'+galaxyname+'NV_fit_dataparam.txt', xran_1,xran_2
  readcol, directoryname+'/'+galaxyname+'/'+galaxyname+'NV_fit_datamodabs.txt', nuvabs,elementsize,NUMLINE=1

;Initializing and fixing params  
  nuvabs=nuvabs[0]
  unity=make_array(elementsize,Value=1)
  elementsize=elementsize[0]
  moduvabs=Make_array(nuvabs,elementsize)

;Rescaling flux values  
  continuum=continuum/1E-14
  flux = flux/1E-14
  moduvabs=moduvabs/1E-14

;Defining plot ranges  
  xran=[xran_1,xran_2]
  plotregion=[value_locate(wave,xran[0]):value_locate(wave,xran[1])]

;Y-range for the normalized spectra plot
  yrannormalized=[0,1.05*Max(normalizedflux[plotregion])]
  
;Y-range for the continuum plot
  yran=[0,1.05*Max(flux[plotregion])]
  
;Normalized plot and legend
  cgplot,wave,normalizedflux,xran=xran,yran=yrannormalized,xstyle=1,ystyle=1,$
    backg='White',axiscolor='Black',color='Black',$
    xtit='Wavelength ($\Angstrom$)',ytit='Normalized F!I$\lambda$',$
    position= [.1,.5,.5,.9], CHARSIZE=.6,thick=1, Title= 'N V'  
  AL_LEGEND,['Intrinsic','ISM','Continuum','Components','Continuum w/Absorption'], $
    Color=['Red','Blue','Orange','Sky Blue','Purple'], charsize=.3, $
    Linestyle=[0,0,0,0,0], Position = [.12,.98],/norm,box=0
    
;Plots absorption features, unity, and outlines spectra
  FOR M=0, nuvabs-1 DO BEGIN
    readcol, directoryname+'/'+galaxyname+'/'+galaxyname+'NV_fit_datamodabs.txt',moduvabs,skipline=1+elementsize*M, NUMLINE=elementsize
    cgoplot,wave,moduvabs,color='Sky Blue',thick=4
  ENDFOR
  cgoplot,wave,unity,color='Orange',thick=4
  cgoplot,wave, modflux, color = 'Purple', thick = 4
  
;Plots absorption and emission labels
  FOR M = 0, N_ELEMENTS(LineWavelength)-1 DO BEGIN
    cgoplot, [LineWavelength[M],LineWavelength[M]], yrannormalized, color = 'Blue', thick=4
    cgoplot, [ShiftedLines[M],ShiftedLines[M]], yrannormalized, color = 'Red', thick=4
    index=Where(LineWavelength lt Max(wave) AND LineWavelength gt Min(wave))
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,LineWavelength[I], $
        .4*yrannormalized[1],  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .4
    ENDFOR
    index=Where(ShiftedLines lt Max(wave) AND ShiftedLines gt Min(wave))
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,ShiftedLines[I], $
        .4*yrannormalized[1],  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .4
    ENDFOR
  ENDFOR

;Plots continuum  
  cgplot,wave,flux,xran=xran,yran=yran,xstyle=1,ystyle=1,$
    backg='White',axiscolor='Black',color='Black',$
    xtickformat="(A1)",ytit='F!I$\lambda$!N/10!E-14!N (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)',$
    position= [.1,.03,.5,.43], /NoErase, CHARSIZE=.6,thick=1
  cgoplot,wave, continuum, color = 'Orange', thick=4
  
;Plots absorption and emission labels
  FOR M = 0, N_ELEMENTS(LineWavelength)-1 DO BEGIN
    cgoplot, [LineWavelength[M],LineWavelength[M]], yran, color = 'Blue', thick=4
    cgoplot, [ShiftedLines[M],ShiftedLines[M]], yran, color = 'Red', thick=4
    index=Where(LineWavelength lt Max(wave) AND LineWavelength gt Min(wave))
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,LineWavelength[I], $
        .4*yran[1],  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .4
    ENDFOR
    index=Where(ShiftedLines lt Max(wave) AND ShiftedLines gt Min(wave))
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,ShiftedLines[I], $
        .4*yran[1],  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .4
    ENDFOR
  ENDFOR

  
ENDIF


;O VI and N V plots

;Read O VI params
IF (keyword_set(OVI)) AND (keyword_set(NV)) THEN BEGIN
  readcol, directoryname+'/'+galaxyname+'/'+galaxyname+'OVI_fit_data.txt', wave, modflux, continuum, flux, normalizedflux
  readcol, directoryname+'/'+galaxyname+'/'+galaxyname+'OVI_fit_dataparam.txt', xran_1,xran_2
  readcol, directoryname+'/'+galaxyname+'/'+galaxyname+'OVI_fit_datamodabs.txt', nuvabs,elementsize,NUMLINE=1

;Initializing and fixing variables  
  nuvabs=nuvabs[0]
  unity=make_array(elementsize,Value=1)
  elementsize=elementsize[0]
  moduvabs=Make_array(nuvabs,elementsize)
  
;Rescaling flux values  
  continuum=continuum/1E-14
  flux = flux/1E-14
  moduvabs=moduvabs/1E-14

;Defining plot ranges  
  xran=[xran_1,xran_2]
  plotregion=[value_locate(wave,xran[0]):value_locate(wave,xran[1])]
  
;Y-range for the normalized spectra plot
  yrannormalized=[0,1.05*Max(normalizedflux[plotregion])]

;Y-range for the continuum plot
  yran=[0,1.05*Max(flux[plotregion])]
  
;Plots normalized spectra
  cgplot,wave,normalizedflux,xran=xran,yran=yrannormalized,xstyle=1,ystyle=1,$
    backg='White',axiscolor='Black',color='Black',$
    xtit='Wavelength ($\Angstrom$)',ytit='Normalized F!I$\lambda$',$
    position= [.1,.5,.5,.9], CHARSIZE=.6,thick=1, Title= 'O VI'
  AL_LEGEND,['Intrinsic','ISM','Continuum','Components','Continuum w/Absorption'], $
    Color=['Red','Blue','Orange','Sky Blue','Purple'], charsize=.3, $
    Linestyle=[0,0,0,0,0], Position = [.12,.98],/norm,box=0
    
;Plots absorption features, unity, and outlines spectra
  FOR M=0, nuvabs-1 DO BEGIN
    readcol, directoryname+'/'+galaxyname+'/'+galaxyname+'OVI_fit_datamodabs.txt',moduvabs,skipline=1+elementsize*M, NUMLINE=elementsize
    cgoplot,wave,moduvabs,color='Sky Blue',thick=4
  ENDFOR
  cgoplot,wave,unity,color='Orange',thick=4
  cgoplot,wave, modflux, color = 'Purple', thick = 4
  
;Plots absorption and emission labels
  FOR M = 0, N_ELEMENTS(LineWavelength)-1 DO BEGIN
    cgoplot, [LineWavelength[M],LineWavelength[M]], yrannormalized, color = 'Blue', thick=4
    cgoplot, [ShiftedLines[M],ShiftedLines[M]], yrannormalized, color = 'Red', thick=4
    index=Where(LineWavelength lt Max(wave) AND LineWavelength gt Min(wave))
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,LineWavelength[I], $
        .4*yrannormalized[1],  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .4
    ENDFOR
    index=Where(ShiftedLines lt Max(wave) AND ShiftedLines gt Min(wave))
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,ShiftedLines[I], $
        .4*yrannormalized[1],  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .4
    ENDFOR
  ENDFOR
  
;Plots continuum
  cgplot,wave,flux,xran=xran,yran=yran,xstyle=1,ystyle=1,$
    backg='White',axiscolor='Black',color='Black',$
    xtickformat="(A1)",ytit='F!I$\lambda$!N/10!E-14!N (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)',$
    position= [.1,.03,.5,.43], /NoErase, CHARSIZE=.6,thick=1
  cgoplot,wave, continuum, color = 'Orange', thick=4, /NoErase

;Plots absorption and emission labels 
  FOR M = 0, N_ELEMENTS(LineWavelength)-1 DO BEGIN
    cgoplot, [LineWavelength[M],LineWavelength[M]], yran, color = 'Blue', thick=4
    cgoplot, [ShiftedLines[M],ShiftedLines[M]], yran, color = 'Red', thick=4
    index=Where(LineWavelength lt Max(wave) AND LineWavelength gt Min(wave))
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,LineWavelength[I], $
        .4*yran[1],  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .4
    ENDFOR
    index=Where(ShiftedLines lt Max(wave) AND ShiftedLines gt Min(wave))
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,ShiftedLines[I], $
        .4*yran[1],  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .4
    ENDFOR
  ENDFOR

;Read N V params  
  readcol, directoryname+'/'+galaxyname+'/'+galaxyname+'NV_fit_data.txt', wave, modflux, continuum, flux, normalizedflux
  readcol, directoryname+'/'+galaxyname+'/'+galaxyname+'NV_fit_dataparam.txt', xran_1,xran_2
  readcol, directoryname+'/'+galaxyname+'/'+galaxyname+'NV_fit_datamodabs.txt', nuvabs,elementsize,NUMLINE=1

;Initializing and fixing variables  
  nuvabs=nuvabs[0]
  unity=make_array(elementsize,Value=1)
  elementsize=elementsize[0]
  moduvabs=Make_array(nuvabs,elementsize)

;Rescaling flux values 
  continuum=continuum/1E-14
  flux = flux/1E-14
  moduvabs=moduvabs/1E-14

;Defining plot ranges  
  xran=[xran_1,xran_2]
  plotregion=[value_locate(wave,xran[0]):value_locate(wave,xran[1])]
  
;Y-range for the normalized spectra plot
  yrannormalized=[0,1.05*Max(normalizedflux[plotregion])]
  
;Y-range for the continuum plot
  yran=[0,1.05*Max(flux[plotregion])]
  
;Plots normalized spectra
  cgplot,wave,normalizedflux,xran=xran,yran=yrannormalized,xstyle=1,ystyle=1,$
    backg='White',axiscolor='Black',color='Black',$
    xtit='Wavelength ($\Angstrom$)',ytit='Normalized F!I$\lambda$',$
    position= [.55,.5,.95,.9], /NoErase, CHARSIZE=.6, Title= 'N V'
    
;Plots absorption features, unity, and outlines spectra
  FOR M=0, nuvabs-1 DO BEGIN
    readcol, directoryname+'/'+galaxyname+'/'+galaxyname+'NV_fit_datamodabs.txt',moduvabs,skipline=1+elementsize*M, NUMLINE=elementsize
    cgoplot,wave,moduvabs,color='Sky Blue',thick=4
  ENDFOR
  cgoplot,wave,unity,color='Orange',thick=4
  cgoplot,wave, modflux, color = 'Purple', thick = 4
  
;Plots absorption and emission labels
  FOR M = 0, N_ELEMENTS(LineWavelength)-1 DO BEGIN
    cgoplot, [LineWavelength[M],LineWavelength[M]], yrannormalized, color = 'Blue', thick=4
    cgoplot, [ShiftedLines[M],ShiftedLines[M]], yrannormalized, color = 'Red', thick=4
    index=Where(LineWavelength lt Max(wave) AND LineWavelength gt Min(wave))
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,LineWavelength[I], $
        .4*yrannormalized[1],  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .4
    ENDFOR
    index=Where(ShiftedLines lt Max(wave) AND ShiftedLines gt Min(wave))
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,ShiftedLines[I], $
        .4*yrannormalized[1],  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .4
    ENDFOR
  ENDFOR

;Plots continuum  
  cgplot,wave,flux,xran=xran,yran=yran,xstyle=1,ystyle=1,$
    backg='White',axiscolor='Black',color='Black',$
    xtickformat="(A1)",ytit='F!I$\lambda$!N/10!E-14!N (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)',$
    position= [.55,.03,.95,.43], /NoErase,Charsize=.6
  cgoplot,wave, continuum, color = 'Orange', thick=2

;Plots absorption and emission labels  
  FOR M = 0, N_ELEMENTS(LineWavelength)-1 DO BEGIN
    cgoplot, [LineWavelength[M],LineWavelength[M]], yran, color = 'Blue', thick=4
    cgoplot, [ShiftedLines[M],ShiftedLines[M]], yran, color = 'Red', thick=4
    index=Where(LineWavelength lt Max(wave) AND LineWavelength gt Min(wave))
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,LineWavelength[I], $
        .4*yran[1],  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .4
    ENDFOR
    index=Where(ShiftedLines lt Max(wave) AND ShiftedLines gt Min(wave))
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,ShiftedLines[I], $
        .4*yran[1],  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .4
    ENDFOR
  ENDFOR
  
ENDIF
CGPS_CLOSE

print, 'Fin.'
END