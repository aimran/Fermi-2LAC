Fermi 2nd-Year AGN Catalog Parser
==================================

(1) This module defines functions for parsing the 2nd-year Fermi LAT AGN Catalog. 
The web page for the catalog is 

  http://heasarc.gsfc.nasa.gov/W3Browse/fermi/fermilac.html

and a full description of the point source analysis and the contents of the 
catalog is provided in http://adsabs.harvard.edu/abs/2011ApJ...743..171A 

The accompanying fits file: fermi_2lac_v01.fits generated using Table 3, 
Table 4, Table 7, Table 8 of the above paper. 

Please contact Asif Imran (asif.imran at icecube.wisc.edu) questions/issues. 

A class is provided to read in the catalog in FITS format and give access to 
the data columns as NumPy arrays. 


(2) I'm also including two ipython notebooks to showcase the use ipython to do
quick analysis. I can not gush enough about ipython notebooks and how awesome it
is. I hope you will agree...


* Author: Asif Imran (aimran@icecube.wisc.edu)
* License: BSD 3-Clause


