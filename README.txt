# History
#  09aug15  DSNR  created
#  13jun14  DSNR  updated libraries and added version comments

VERSIONS FOR DIFFERENT INSTRUMENTS

/gmos: Gemini Multi-Object Spectrograph (GMOS)
       Optimized for ~5500-7500 A IFU data
/lris: Keck Low Resolution Imaging Spectrograph (LRIS)
       Optimized for red+blue arms, old CCD, higher-res gratings
/nirsifs: Gemini Near-Infrared Integral Field Spectrograph (NIFS) or Keck OSIRIS
/slitspec: Optimized for long-slit spectra
/sp1: Optimized for single spectra (e.g., SDSS)

SETTING UP UHSPECFIT FOR A GENERIC INSTRUMENT

[In the following instructions, I quote for example the files for a
specific instrument, LRIS.]

1. Download required libraries.

   IDLUTILS
   http://www.sdss3.org/dr8/software/idlutils.php

   MPFIT
   http://www.physics.wisc.edu/~craigm/idl/idl.html

   Coyote
   http://www.idlcoyote.com/documents/programs.php#COYOTE_LIBRARY_DOWNLOAD

2. Start by customizing the wrapper program that actually calls the
fit for each spectrum [LRIS_FIT_SPECTRA].

3. Next, customize the routine that initializes the parameter guesses
and constraints [LRIS_INITPARINFO].  For LRIS, this consists of a
series of global parameters, including redshifts and line widths,
followed by the 3 parameters of a Gaussian for each emission line and
velocity component.  NOTE: the "ppoff0" parameter must be set to
reflect the number of global parameters.

4. Set up the linelist for the emission lines you want to fit
[LRIS_INITLINELIST].

5. Set up the continuum fitting routine [LRIS_FIT_CONTINUUM].

6. If desired, customize the routine that sets variables to initialize
the fit [LRIS_INITFIT_SPECTRA].

7. Finally, customize the routine to analyze the resulting fits
[LRIS_ANL_SPECTRA].

-------

NOTES

1. MANYGAUSS, the routine that evaluates the emission line Gaussians,
assumes constant dispersion (in A/pix) and that there are no "holes"
in the spectrum.  If you have a variable dispersion or gaps in
coverage, either modify MANYGAUSS to fit your needs or use
MANYGAUSS_SLOW (which, as its name implies, is slower).
