getslit1.pro generates a slit from 10257 data cube
this command line program uses these new routines:
     slitcens.pro
     cornersbox.pro
     overlap.pro

slitcens.pro:
	Given cube, slit PA & center & element length generates
	list of slit center x & y coordinates in cube spaxel grid

cornersbox.pro
	fed input from slitcens as well as angle & cube
	gets coordinates for corners of rectangle of given length,
	width, angle and center-get corners of slit elements
	used to define polygons that are the slit elements

overlap.pro
	uses input from cornersbox (corners of each slit element,
	the polygons defining the slit element) and data cube
	spaxel grid. Generates list of overlapping spaxels
	and percent overlap of each spaxel for a given slit elemnt
	
	This list can be used to generate the simple weighted
	spectrum by coadding weighted spectra from all spaxels
	overlapping a spaxel element.


testslits?.pro:
	used for testing development of this code. fun for plotting.