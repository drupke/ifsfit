; docformat = 'rst'
;
;+
;
; Produces plots showing the entire data spectra, as well as zoomed-in
; portions.
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
PRO ifsf_spectraplots, initfile, directoryname, galaxyname


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
    
;Read galaxy full names   
  readcol,directoryname+'/'+'Redshifts', galaxyfullnamelist, format='(A)'

;Read wavelength and flux
  readcol, directoryname+'/'+galaxyname+'/'+galaxyname+'.txt', wave, flux

;Read redshift
  readcol, directoryname+'/'+initfile, galaxynamelist,redshiftlist, SKIPLINE=1,format = '(A,D)'
  selectionparameter=WHERE(galaxynamelist eq galaxyname)
  galaxyname=galaxynamelist[selectionparameter[0]]
  galaxyfullname=galaxyfullnamelist[selectionparameter[0]]
  zsys=redshiftlist[selectionparameter[0]]

;Set ASPECT RATIO for plots. Equivalent to (y-range)/(x-range) in data coords.
  aratio=(.11)/(.9)
  
;Shifts the absorption/emission waveelngths by the galaxy's redshift
  ShiftedLines= LineWavelength*(1+zsys)
  
;Rescaling flux values
  relativeflux=flux/1E-14
  
;Set the y-range of the big spectra
  yran=[0,Max(relativeflux)]

;Acquires the range of the wavelength values
  baserange=Max(wave)-Min(wave)
  
;Sets window size in inches
  aspectRatio = FLOAT(!D.Y_VSIZE) / !D.X_VSIZE
  xsize=7.5
  ysize=8.5 
  IF ysize GT 10.5 THEN BEGIN
    ysize = 10.5
    xsize = ysize / aspectRatio
  ENDIF  
  
  xoffset=(8.5-xsize)/2.0d
  yoffset=(11.0-ysize)/2.0d
  
  
;Opens a Postscript value to draw plots on
  cgPS_OPEN,directoryname+'/'+galaxyname+'/'+galaxyname+'spectraplot',/Encapsulated,scale_factor=1, $
    Charsize=.15,font=-1, /NOMATCH
  DEVICE, xsize=xsize, ysize=ysize, xoffset=xoffset, yoffset=yoffset, /Inches
  !Y.MINOR=1
  !Y.THICK=.8
  !X.THICK=.8
  !Y.TICKFORMAT='(F5.1)'
  
  cgText, .5,.98,galaxyfullname +','+ 'z='+String(zsys), alignment=.5, Charsize = 1
  
;Full range plot
  cgplot, wave, relativeflux, xstyle=1, ystyle=0, yran=yran, $
    axiscolor='Black',color='Black',$
    xtit='Observe Wavelength ($\Angstrom$)' ,ytit='F!I$\lambda$!N/10!E-14!N (ergs s!E-1 !Ncm!E-2 !N$\Angstrom$!E-1!N)', $
    Position = [.05,.86,.95,.97], CHARSIZE=.35,thick=1, aspect=aratio,charthick=.4,/NoErase


  
;Plot legend
  AL_LEGEND,['Intrinsic','ISM'], $
    Color=['Red','Blue'], charsize=.3, charthick=.4, $
    Linestyle=[0,0], Position = [Min(wave)+.04*(Max(wave)-Min(wave)),.9*(Max(relativeflux)-Min(relativeflux))], $
    bthick=.6,clear=1,background_color='White'


;Plots 6 zoomed-in regions, one after another, along with absorption/emission lines and labels

  cgplot, wave, relativeflux, xstyle=1, ystyle=1, xran=[Min(wave),Min(wave)+(1d/6)*baserange], $
    axiscolor='Black',color='Black', $
    yran=[0, $
    3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(0d/6)*baserange): $
    value_locate(wave,Min(wave)+(1d/6)*baserange)])^2))], $
    Position = [.05,.715,.95,.825], /NoErase, CHARSIZE=.35,thick=1, aspect=aratio,charthick=.4
  FOR M = 0, N_ELEMENTS(LineWavelength)-1 DO BEGIN
    cgoplot, [LineWavelength[M],LineWavelength[M]], 2*yran, color = 'Blue', thick=1
    cgoplot, [ShiftedLines[M],ShiftedLines[M]], 2*yran, color = 'Red', thick=1
    xran=[Min(wave),Min(wave)+(1d/6)*baserange]
    index=Where(LineWavelength lt xran[1] AND LineWavelength gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,LineWavelength[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(0d/6)*baserange): $
        value_locate(wave,Min(wave)+(1d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .25, $
        charthick=.4
    ENDFOR
    index=Where(ShiftedLines lt xran[1] AND ShiftedLines gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,ShiftedLines[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(0d/6)*baserange): $
        value_locate(wave,Min(wave)+(1d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .25, $
        charthick=.4
    ENDFOR
  ENDFOR

  cgplot, wave, relativeflux, xstyle=1, ystyle=1, xran=[Min(wave)+(1d/6)*baserange,Min(wave)+(2d/6)*baserange], $
    axiscolor='Black',color='Black', $
    yran=[0, $
    3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(1d/6)*baserange): $
    value_locate(wave,Min(wave)+(2d/6)*baserange)])^2))], $
    Position = [.05,.575,.95,.685], /NoErase, CHARSIZE=.35,thick=1, aspect=aratio,charthick=.4
  FOR M = 0, N_ELEMENTS(LineWavelength)-1 DO BEGIN
    cgoplot, [LineWavelength[M],LineWavelength[M]], 2*yran, color = 'Blue', thick=1
    cgoplot, [ShiftedLines[M],ShiftedLines[M]], 2*yran, color = 'Red', thick=1
    xran=[Min(wave)+(1d/6)*baserange,Min(wave)+(2d/6)*baserange]
    index=Where(LineWavelength lt xran[1] AND LineWavelength gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,LineWavelength[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(1d/6)*baserange): $
        value_locate(wave,Min(wave)+(2d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .25, $
        charthick=.4
    ENDFOR
    index=Where(ShiftedLines lt xran[1] AND ShiftedLines gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,ShiftedLines[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(1d/6)*baserange): $
        value_locate(wave,Min(wave)+(2d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .25, $
        charthick=.4
    ENDFOR
  ENDFOR

  cgplot, wave, relativeflux, xstyle=1, ystyle=1, xran=[Min(wave)+(2d/6)*baserange,Min(wave)+(3d/6)*baserange], $
    axiscolor='Black',color='Black', $
    yran=[0, $
    3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(2d/6)*baserange): $
    value_locate(wave,Min(wave)+(3d/6)*baserange)])^2))], $
    Position = [.05,.435,.95,.545], /NoErase, CHARSIZE=.35,thick=1, aspect=aratio,charthick=.4
  FOR M = 0, N_ELEMENTS(LineWavelength)-1 DO BEGIN
    cgoplot, [LineWavelength[M],LineWavelength[M]], 2*yran, color = 'Blue', thick=1
    cgoplot, [ShiftedLines[M],ShiftedLines[M]], 2*yran, color = 'Red', thick=1
    xran=[Min(wave)+(2d/6)*baserange,Min(wave)+(3d/6)*baserange]
    index=Where(LineWavelength lt xran[1] AND LineWavelength gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,LineWavelength[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(2d/6)*baserange): $
        value_locate(wave,Min(wave)+(3d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .25, $
        charthick=.4
    ENDFOR
    index=Where(ShiftedLines lt xran[1] AND ShiftedLines gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,ShiftedLines[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(2d/6)*baserange): $
        value_locate(wave,Min(wave)+(3d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .25, $
        charthick=.4
    ENDFOR
  ENDFOR

  cgplot, wave, relativeflux, xstyle=1, ystyle=1, xran=[Min(wave)+(3d/6)*baserange,Min(wave)+(4d/6)*baserange], $
    axiscolor='Black',color='Black', $
    yran=[0, $
    3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(3d/6)*baserange): $
    value_locate(wave,Min(wave)+(4d/6)*baserange)])^2))], $
    Position = [.05,.295,.95,.405], /NoErase, CHARSIZE=.35,thick=1, aspect=aratio,charthick=.4
  FOR M = 0, N_ELEMENTS(LineWavelength)-1 DO BEGIN
    cgoplot, [LineWavelength[M],LineWavelength[M]], 2*yran, color = 'Blue', thick=1
    cgoplot, [ShiftedLines[M],ShiftedLines[M]], 2*yran, color = 'Red', thick=1
    xran=[Min(wave)+(3d/6)*baserange,Min(wave)+(4d/6)*baserange]
    index=Where(LineWavelength lt xran[1] AND LineWavelength gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,LineWavelength[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(3d/6)*baserange): $
        value_locate(wave,Min(wave)+(4d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .25, $
        charthick=.4
    ENDFOR
    index=Where(ShiftedLines lt xran[1] AND ShiftedLines gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,ShiftedLines[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(3d/6)*baserange): $
        value_locate(wave,Min(wave)+(4d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .25, $
        charthick=.4
    ENDFOR
  ENDFOR

  cgplot, wave, relativeflux, xstyle=1, ystyle=1, xran=[Min(wave)+(4d/6)*baserange,Min(wave)+(5d/6)*baserange], $
    axiscolor='Black',color='Black', $
    yran=[0, $
    3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(4d/6)*baserange): $
    value_locate(wave,Min(wave)+(5d/6)*baserange)])^2))], $
    Position = [.05,.155,.95,.265], /NoErase, CHARSIZE=.35,thick=1, aspect=aratio,charthick=.4
  FOR M = 0, N_ELEMENTS(LineWavelength)-1 DO BEGIN
    cgoplot, [LineWavelength[M],LineWavelength[M]], 2*yran, color = 'Blue', thick=1
    cgoplot, [ShiftedLines[M],ShiftedLines[M]], 2*yran, color = 'Red', thick=1
    xran=[Min(wave)+(4d/6)*baserange,Min(wave)+(5d/6)*baserange]
    index=Where(LineWavelength lt xran[1] AND LineWavelength gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,LineWavelength[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(4d/6)*baserange): $
        value_locate(wave,Min(wave)+(5d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .25,$
        charthick=.4
    ENDFOR
    index=Where(ShiftedLines lt xran[1] AND ShiftedLines gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,ShiftedLines[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(4d/6)*baserange): $
        value_locate(wave,Min(wave)+(5d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .25,$
        charthick=.4
    ENDFOR
  ENDFOR
  
  cgplot, wave, relativeflux, xstyle=1, ystyle=1, xran=[Min(wave)+(5d/6)*baserange,Min(wave)+(6d/6)*baserange], $
    axiscolor='Black',color='Black', $
    yran=[0, $
    3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(5d/6)*baserange): $
    value_locate(wave,Min(wave)+(6d/6)*baserange)])^2))], $
    Position = [.05,.02,.95,.13], /NoErase, CHARSIZE=.35,thick=1, aspect=aratio,charthick=.4
  FOR M = 0, N_ELEMENTS(LineWavelength)-1 DO BEGIN
    cgoplot, [LineWavelength[M],LineWavelength[M]], 2*yran, color = 'Blue', thick=1
    cgoplot, [ShiftedLines[M],ShiftedLines[M]], 2*yran, color = 'Red', thick=1
    xran=[Min(wave)+(5d/6)*baserange,Min(wave)+(6d/6)*baserange]
    index=Where(LineWavelength lt xran[1] AND LineWavelength gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,LineWavelength[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(5d/6)*baserange): $
        value_locate(wave,Min(wave)+(6d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .25,$
        charthick=.4
    ENDFOR
    index=Where(ShiftedLines lt xran[1] AND ShiftedLines gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,ShiftedLines[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(5d/6)*baserange): $
        value_locate(wave,Min(wave)+(6d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = .25,$
        charthick=.4
    ENDFOR
  ENDFOR
  CGPS_CLOSE,Allow_Transparent=1
END