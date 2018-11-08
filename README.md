# IR-LIME
Line radiative transfer package for infrared transitions, built upon the existing code LIME.


## Prerequisites
First, install LIME 1.43 from Christian Brinch's page: https://bitbucket.org/brinch/limepackage
Note that this requires a C compiler as well as three library packages: qhull, gsl, and cfitsio. 

Second, if interested in generating new grid input files for IR-LIME, install RADLite from Klaus Pontoppidan's page here: https://bitbucket.org/pontoppi/radlite


## What's Included
This package is in early stages of development, but it is intended to be a branch/modification of LIME 1.43 focused on the treatment of 3D, non-LTE radiative transfer in the infrared. While many future upgrades are planned, at the moment this holds only a basic working example of functionality.

Included in this repository are python scripts to generate parameter grids for input to LIME, a set of input files for the grid creation, and modified LIME source files to implement these changes. Additionally, an example model is included to demonstrate functionality.

