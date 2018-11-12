# IR-LIME
Line radiative transfer package for infrared transitions, built upon the existing code LIME.

This package is in early stages of development, but it is intended to be a branch/modification of LIME 1.43 focused on the treatment of 3D non-LTE radiative transfer in infrared wavelengths. While many future upgrades are planned, at the moment this holds only a basic working example of functionality.


## What's Included
Included in this repository are python scripts to generate parameter grids for input to LIME, a set of input files for the grid creation, and modified LIME source files to implement these changes. Additionally, an example model is included to demonstrate functionality.

## Prerequisites
First, install LIME 1.43 from Christian Brinch's page: https://bitbucket.org/brinch/limepackage
Note that this requires a C compiler as well as three library packages: qhull, gsl, and cfitsio. Installation instructions for these packages are given in the LIME documentation.

Second, if interested in generating new grid input files for IR-LIME, you may want to install RADLite from Klaus Pontoppidan's page here: https://bitbucket.org/pontoppi/radlite


## Installation of IR-LIME
IR-LIME exists now as a slight modification of base LIME 1.43 - this means that installation involves replacing base LIME source files with those provided here. I would suggest you test the LIME installation with its provided example model, then archive the replaced LIME files when replacing them with IR-LIME. Then, the example provided here can be tested.

### Notes
I ran into an issue during LIME installation - the Makefile assumes your C compiler is gcc and that it is compatible with the -fopenmp flag. 
