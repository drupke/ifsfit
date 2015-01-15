-------------------------------------------------------------------------
OVERVIEW
-------------------------------------------------------------------------

IFSFIT is a general-purpose library for fitting the continuum,
emission lines, and absorption lines in integral field spectra
(IFS). It uses PPXF (the Penalized Pixel-Fitting method developed by
Cappellari & Emsellem, 2004, PASP, 116, 138) to find the best fit
stellar continuum (using a user-defined library of stellar templates
and including additive polynomials), or optionally a user-defined
method to find the best fit continuum. It uses MPFIT to simultaneously
fit Gaussians to any number of emission lines and emission line
velocity components. It will also fit the NaI D feature using analytic
absorption and/or emission-line profiles.

-------------------------------------------------------------------------
REQUIREMENTS
-------------------------------------------------------------------------

IDL v8.0 or higher (tested with v8.3)

IDL libraries:
- IDL Astronomy User's Library, for various routines
  http://idlastro.gsfc.nasa.gov
- MPFIT, for non-linear least-squares fitting
  http://www.physics.wisc.edu/~craigm/idl/idl.html
- Coyote, for graphics
  http://www.idlcoyote.com/documents/programs.php#COYOTE_LIBRARY_DOWNLOAD
  [or from the subversion repository: https://code.google.com/p/idl-coyote/]
- PPXF
  http://www-astro.physics.ox.ac.uk/~mxc/software/#ppxf
- IDLUTILS, primarily for structure manipulation tasks (e.g.,
  STRUCT_ADDTAGS).
  http://www.sdss3.org/dr8/software/idlutils.php

Note that the IDL Astronomy User's Library ships with some Coyote
routines, and IDLUTILS ships with the IDL Astronomy User's Library and
MPFIT. However, it's not clear how well these libraries keep track of
each other, so it may be preferable to download each package
separately and delete the redundant routines that ship within other
packages.

To fit stellar continua, templates are required. E.g., the population
synthesis models from Gonzalez-Delgado et al. (2005, MNRAS, 357, 945)
are available at

http://www.iaa.csic.es/~rosa/research/synthesis/HRES/ESPS-HRES.html

The enclosed routine IFSF_GDTEMP can be used to convert these tables
into a form readable by IFSF.

-------------------------------------------------------------------------
USAGE
-------------------------------------------------------------------------

Usage is in principle straigtforward. The continuum and emission-line
fitting is run from the command line as

> IFSF,'initialization_file' [,cols=[low,high],rows=[low,high], etc.]

and then the results are processed using

> IFSFA,'initialization_file' [,cols=[low,high],rows=[low,high], etc.]

The latter produces plots of the lines fit, a table of emission line
parameters and a table of fit results, and an IDL save file containing
a "data cube" of emission line parameters. If desired, it will also
normalize the continuum around NaI D 5890, 5896 and estimate various
parameters of this feature.

The initialization file defines an IDL function that returns an
initialization structure, INITDAT. The tags of this structure, and a
description of each one, are found in INITTAGS.txt. Several tags are
required, but most are optional.

The initialization function should also define two other structures,
INITMAPS and INITNAD, and also call these as keyword arguments. The
first of these, INITMAPS, controls how IFSF_MAKEMAPS processes the
output of IFSFA into parameter maps. This routine produces various
emission-line, continuum, and absorption-line maps. The possible tags
for INITMAPS are described in INITTAGS.txt. If IFSF_MAKEMAPS is not
used, INITMAPS can be set to an empty structure.

The second of these other initialization structures, INITNAD, controls
how the region around the NaI D feature is fit using IFSF_FITNAD. This
routine fits emission and absorption-line models to HeI 5876 and NaI D
5890, 5896 and produces plots of the fits and a "data cube" with fit
parameters. The routine also estimates errors using Monte Carlo
methods. The error estimation can be sped up using multi-core parallel
processing. The possible tags for INITNAD are described in
INITTAGS.txt. If IFSF_MAKEMAPS is not used, INITNAD can be set to an
empty structure.

There are a considerable number of knobs that can be turned to
optimize/customize the fits and customize the outputs. An example
initialization file is included (init/IFSF_F05189.pro).

The other user-modifiable initialization procedure that is required is
one that sets the initial emission-line parameter structure for input
into MPFIT. Included in this version is a procedure that is optimized
for the GMOS instrument (init/IFSF_GMOS.pro) and contains line
parameters pertaining to strong emission lines in the wavelength range
4000-7000 A. However, this code can be relatively easily adapated to
other instruments and to include parameters (like fixed line ratios)
for other emission lines. The code has in the past been used to
successfully fit data from many instruments: GMOS, LRIS, NIFS, OSIRIS,
and WiFeS.

-------------------------------------------------------------------------
IMPORTANT NOTES
-------------------------------------------------------------------------

IFSF_MANYGAUSS, the routine that evaluates the emission line
Gaussians, assumes constant dispersion (in A/pix) and that there are
no "holes" in the spectrum (i.e., that each wavelength follows the
next by simply adding the dispersion).  If you have a variable
dispersion or gaps in your wavelength array, either modify
IFSF_MANYGAUSS to fit your needs or use IFSF_MANYGAUSS_SLOW (which, as
its name implies, is slower).

-------------------------------------------------------------------------
QUESTIONS? BUGS? WANT TO MODIFY THE CODE?
-------------------------------------------------------------------------

Feel free to contact David Rupke at drupke@gmail.com with questions,
bug reports, etc.

Modifications are encouraged, but subject to the license.

-------------------------------------------------------------------------
LICENSE AND COPYRIGHT
-------------------------------------------------------------------------

Copyright (C) 2013-2014 David S. N. Rupke

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License or any
later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see http://www.gnu.org/licenses/.
