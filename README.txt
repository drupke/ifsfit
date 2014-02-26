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
velocity components.

-------------------------------------------------------------------------
REQUIREMENTS
-------------------------------------------------------------------------

IDL v8.0 or higher (tested with v8.3)

IDL libraries:
- IDL Astronomy User's Library
  http://idlastro.gsfc.nasa.gov
- MPFIT
  http://www.physics.wisc.edu/~craigm/idl/idl.html
- Coyote
  http://www.idlcoyote.com/documents/programs.php#COYOTE_LIBRARY_DOWNLOAD
  [or from the subversion repository: https://code.google.com/p/idl-coyote/]
- PPXF
  http://www-astro.physics.ox.ac.uk/~mxc/software/#ppxf

To fit stellar continua, templates are required. E.g., the population
synthesis models from Gonzalez-Delgado et al. (2005, MNRAS, 357, 945)
are available at

http://www.iaa.csic.es/~rosa/research/synthesis/HRES/ESPS-HRES.html

The enclosed routine IFSF_GDTEMP can be used to convert these tables
into a form readable by IFSF.

-------------------------------------------------------------------------
USAGE
-------------------------------------------------------------------------

Usage is in principle straigtforward. The routine is run from the
command line as

> IFSF,'initialization_file' [,cols=[low,high],rows=[low,high], etc.]

and then the results are processed using

> IFSFA,'initialization_file' [,cols=[low,high],rows=[low,high], etc.]

The latter produces plots of the lines fit, a table of emission line
parameters and a table of fit results, and an IDL save file containing
a "data cube" of emission line parameters.

The initialization file defines an IDL function that returns an
initialization structure. The tags of this structure, and a
description of each one, are found in INITTAGS.txt. Several tags are
required, but most are optional.

In practice, there are a considerable number of knobs that can be
turned to optimize/customize the fit and customize the output. An
example initialization file is included (init/IFSF_F05189.pro).

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

Copyright (C) 2014 David S. N. Rupke

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
