/*
 *  grid.c
 *  LIME, The versatile 3D line modeling environment
 *
 *  Created by Christian Brinch on 11/16/06.
 *  Copyright 2006-2017, Christian Brinch,
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */

#ifdef QHULL_INC_QHULL
#include <qhull/qhull_a.h>
#else
#include <libqhull/qhull_a.h>
#endif
#include "lime.h"
#include "radlite.h"

void
write_defgrid(inputPars *par, struct grid *g){
  FILE *fp;
  double hfrac,h2frac,densum;
  int i,j,l=0;
  char flags[255];

  if((fp=fopen(par->defgridfile, "w"))==NULL){
    if(!silent) bail_out("Error writing grid file!");
    exit(1);
  }

  for(i=0; i<par->pIntensity; i++) {
    densum=0.;
    for(j=0; j<par->collPart; j++) {
      densum+=g[i].dens[j];
    }
    /*hfrac = g[i].dens[2]/densum;
    h2frac = 1.00E+0-hfrac;
    */
    //fprintf(fp,"%.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e\n", g[i].x[0], g[i].x[1], g[i].x[2], densum, g[i].abun[0], g[i].t[0], g[i].t[1], g[i].dopb, g[i].vel[0], g[i].vel[1], g[i].vel[2], hfrac, h2frac);
    //fprintf(fp,"%e %e %e %e %e %e %e %e %e %e %e %e %e\n", g[i].x[0], g[i].x[1], g[i].x[2], densum, g[i].abun[0], g[i].t[0], g[i].t[1], g[i].dopb, g[i].vel[0], g[i].vel[1], g[i].vel[2], hfrac, h2frac);
    fprintf(fp,"%e %e %e %e %e %e %e %e %e %e %e\n", g[i].x[0], g[i].x[1], g[i].x[2], densum, g[i].abun[0], g[i].t[0], g[i].t[1], g[i].dopb, g[i].vel[0], g[i].vel[1], g[i].vel[2]);
  }
  fprintf(fp, "\n");

  fclose(fp);
}

