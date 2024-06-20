; docformat = 'rst'
;
;+
;
; Load BPASS population synthesis models into a
; structure of fluxlength and flux arrays. Population synthesis models
; from Byrne & Stanway 2023.  See
; https://bpass.auckland.ac.nz/
; for more details.
;
; From p 14 of the v2.3 manual:
; These files contain the primary output of BPASS, which is the stellar spectral
;  energy distribution (SED). Flux values are given every 1 Angstrom in the 
;  range 1 - 100,000 A. Most users will wish to resample to lower resolution, 
;  depending on their use case. We caution that some of the stellar atmospheres 
;  we use also have slightly lower spectral resolution.
;  Each file has 52 columns and 105 rows. The first column lists a wavelength 
;  in angstroms, and each remaining column n (n>1) holds the model flux for the 
;  population at an age of 10^(6+0.1*(n-2)) years at that wavelength.
;  The units of flux are Solar Luminosities per Angstrom, normalised for a 
;  cluster of 1e6 Msun formed in a single instantaneous burst. The total 
;  luminosity of the SED can be simply calculated by summing all the rows 
;  together. Or the flux over the wavelength range from, for example, 2000 to 
;  3000 Angstroms can be calculated by summing the 2000th to 3000th rows.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    IDL save file. The file consists of a structure (named
;    'template') containing tags: lambda, flux, ages, zs. Lambda
;    is an n-element array, while flux is and nxm array of the format
;    [n_waves,nages*nz]. Ages and zs are the arrays of age and metallicity values.
;
; :Params:
;    indir: in, required, type=string
;      Path of input model FITS files.
;    outfile: in, required, type=string
;      Filename and path of output file.
;    wavelo: in, required, type=double
;    wavehi: in, required, type=double
;      Lower and upper wavelengths for output templates, in integer Angstroms
;
; :Keywords:
;    binary: in, optional, type=bool
;      Single models assumed; switch this on to include binaries
;    zlist: in, optional, type=dblarr
;      metallicities to include
; 
; :Author:
;    David S. N. Rupke::
;      Rhodes College
;      Department of Physics
;      2000 N. Parkway
;      Memphis, TN 38104
;      drupke@gmail.com
;
; :History:
;    ChangeHistory::
;      2024mar08, DSNR, created
;    
; :Copyright:
;    Copyright (C) 2024 David S. N. Rupke
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
pro ifsf_bpasstemp,indir,outfile,wavelo,wavehi,binary=binary,zlist=zlist

   bad = 1d99
   ; metallicities
   if ~ keyword_set(zlist) then $
      zlist = [0.001,0.002,0.003,0.004,0.006,0.008,0.010,0.014,0.020,0.030,0.040]
   nz = n_elements(zlist)

   ; ages for one metallicity
   nages = 51
   tonez = dblarr(nages)
   ; log age in yr
   for i=0,nages-1 do tonez[i] = 10d^(6d + 0.1d*i)

   tall = reform(transpose(rebin(transpose(tonez),nz,nages)),nages*nz)
   zall = rebin(zlist,nages*nz,/sample)

   ; wavelength spacing is 1 A; so number of wavlengths is:
   nwave = wavehi-wavelo+1
   waveall = dindgen(nwave)+wavelo

   ; output flux array
   fluxall = dblarr(nwave,nz*nages)
   
   ; to hold temporary SED points at each wavelength and all ages
   spectmp = dblarr(nages)
   ; cycle through files
   foreach zval,zlist,iz do begin
      ; single or binary
      sinorbin='sin'
      if keyword_set(binary) then sinorbin='bin'
      ; convert z val to string
      zvalstr = strsplit(string(zval,format='(D0.3)'),'.',/extract)
      filename = indir+$
         'spectra-'+sinorbin+'-imf135_300.a+00.z'+zvalstr[1]+'.dat'
      openr,readlun,filename,/get_lun
      wavetmp = 0d
      iflux = 0
      while wavetmp le wavehi do begin
         READF, readlun, format='(E16.7,'+string(nages)+'E16.7)',wavetmp,spectmp
         if wavetmp ge wavelo AND wavetmp le wavehi then begin
            fluxall[iflux,iz:iz+nages-1] = spectmp
            iflux++
         endif
      endwhile
      free_lun,readlun
;      for i=0,nages-1 do begin
;         fluxall[*,iz:iz+i] /= median(fluxall[*,iz:iz+i])
;      endfor
   endforeach
   
   template = {lambda: waveall, flux: fluxall, ages: tall, zs: zall}
   save,template,filename=outfile
      
end
