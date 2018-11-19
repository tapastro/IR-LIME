/*
 *  model.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
 */

#include "lime.h"
#include "radlite.h"
/******************************************************************************/

/* Added by TAP to select closest RADMC grid point */
/* Moved nearestGridPt to radlite.h */


void
input(inputPars *par, image *img){
  int i;

  /*
   * Basic parameters. See cheat sheet for details.
   */
  par->radius                   = 0.9*AU;
  par->minScale                 = 0.13*AU;
  par->pIntensity               = 36000;
  par->sinkPoints               = 3000;


  par->dust                     = "radlite_dust.tab";
  //par->moldatfile[0]            = "co_zz_numax2.dat";//"12CO_H_H2@lamda.dat";
  par->moldatfile[0]            = "co_zz_numax2_jmax21.dat";//"12CO_H_H2@lamda.dat";
  par->antialias                = 1;
  par->sampling                 = 0; // log distr. for radius, directions distr. uniformly on a sphere.
  par->outputfile               = "populations_nu2j21_gtd100_abun1_iter12_tgas_allcofrac.pop";
  //par->binoutputfile            = "restart.pop";
  //par->gridfile                 = "grid_interp_nu2j21_gtd100_abun1_iter12_tgas.vtk";
  par->lte_only                 = 0;

  /*
   * Definitions for image #0. Add blocks with successive values of i for additional images.
  */ 
  i=0;
  img[i].nchan                  = 500;             // Number of channels
  img[i].velres                 = 3000.;           // Channel resolution in m/s
  //img[i].trans                  = 213;              // zero-indexed J quantum number
  img[i].freq                   = 6.33E+13;       // Center frequency of bandwidth
  //img[i].bandwidth               = 4.02734E+11;              // zero-indexed J quantum number
  img[i].pxls                   = 200;            // Pixels per dimension
  img[i].imgres                 = 0.008;            // Resolution in arc seconds
  img[i].distance               = 1*PC;         // source distance in m
  img[i].source_vel             = 0;              // source velocity in m/s
  img[i].theta 			        = PI*37./180.;
   /*
    For each set of image parameters above, numerous images with different units can be outputted. This is done by
   * setting img[].units to a delimited (space, comma, colon, underscore) string of the required outputs, where:
   *        0:Kelvin
   *        1:Jansky/pixel
   *        2:SI
   *        3:Lsun/pixel
   *        4:tau
   * If multiple units are specified for a single set of image parameters (e.g. "0 1 2") then the unit name is added
   * automatically at the end of the given filename, but before the filename extension if it exists. Otherwise if a
   * single unit is specified then the filename is unchanged.

   * A single image unit can also be specified for each image using img[].unit as in previous LIME versions. Note that
   * only img[].units or img[].unit should be set for each image.
    */  
  img[i].unit                   = 1;
  img[i].filename               = "co_line_nu2j21_gtd100_abun1_iter12_200px_1p5au_0p75tgas_allcofrac.fits";   // Output filename

}


/******************************************************************************/

void
density(double x, double y, double z, double *density){
    double tmp, den, h, h2;

    den = 1.e3/AMU*bilinInterpVal(x,y,z,DENCOL);

    h = 2.*bilinInterpVal(x,y,z,HCOL);
    h2 = 2.*bilinInterpVal(x,y,z,H2COL);

    /* Density profiles for H2 and electrons - very rough distribution */
    density[0] = 0.25*den*h2;
    density[1] = 0.75*den*h2;
    density[2] = den*h;
    //density[3] = 0.1*den;
}

/******************************************************************************/

void
temperature(double x, double y, double z, double *temperature){
    double temratio, tem,dtem;
    temratio = 0.75;
    tem = bilinInterpVal(x,y,z,TEMCOL);
    dtem = bilinInterpVal(x,y,z,DTEMCOL);
    temperature[0]=((temratio*tem)+((1.-temratio)*dtem));
    temperature[1]=bilinInterpVal(x,y,z,DTEMCOL);
    //temperature[0]=ingridarray[5][nearestGridPt(x,y,z)];
    //temperature[1]=ingridarray[6][nearestGridPt(x,y,z)];
}

/******************************************************************************/

void
abundance(double x, double y, double z, double *abundance){
  /*
   * Here we use a constant abundance. Could be a
   * function of (x,y,z).
   */
    abundance[0] = 1.*bilinInterpVal(x,y,z,ABNCOL);
}

/******************************************************************************/

void
doppler(double x, double y, double z, double *doppler){
  /*
   * 200 m/s as the doppler b-parameter. This
   * can be a function of (x,y,z) as well.
   * Note that *doppler is a pointer, not an array.
   * Remember the * in front of doppler.
   */
  *doppler = bilinInterpVal(x,y,z,DOPCOL);
}

/******************************************************************************/

void
velocity(double x, double y, double z, double *vel){
  double gridv,phi,r;
  double mstar=2.9829e30;

  r = sqrt(pow(x,2)+pow(y,2));
  /*Shortcut taken - assumed that grid velocity non-zero only in y direction */
  gridv = sqrt(GRAV*mstar/r);
  phi = acos(x/sqrt(pow(x,2)+pow(y,2)));

  vel[0] = -sin(phi)*gridv;
  vel[1] = cos(phi)*gridv;
  vel[2] = 0.;
}

/******************************************************************************/

