# IR-LIME
Line radiative transfer package for infrared transitions, built upon the existing code LIME.

This package is in early stages of development, but it is intended to be a branch/modification of LIME 1.43 focused on the treatment of 3D non-LTE radiative transfer in infrared wavelengths. While many future upgrades are planned, at the moment this repo contains a basic working example of functionality.


## What's Included
Included in this repository are python scripts to generate parameter grids for input to LIME, a set of input files for the grid creation, and modified LIME source files to implement these changes. Additionally, an example model is included to demonstrate functionality.

## Prerequisites
First, follow the package installation instructions for qhull, gsl, and cfitsio from LIME 1.43 on Christian Brinch's page: https://bitbucket.org/brinch/limepackage
Note that this requires a C compiler, ideally gcc 6.xx or newer. Then download the files here, and follow the installation instructions for LIME 1.43 with these files. Not much has changed in the installation steps between branches.

Second, if interested in generating new grid input files for IR-LIME, you may want to install RADLite from Klaus Pontoppidan's page here: https://bitbucket.org/pontoppi/radlite


## Installation of IR-LIME
IR-LIME exists now as a slight modification of base LIME 1.43.

### Notes
The Makefile assumes your C compiler is gcc and that it is compatible with the -fopenmp flag. Additionally, ensure you change the file locations in model.c to the location of the files in your local installaion.
