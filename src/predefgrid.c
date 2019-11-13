/*
 *  predefgrid.c
 *  LIME, The versatile 3D line modeling environment
 *
 *  Created by Christian Brinch on 08/26/10.
 *  Copyright 2006-2017, Christian Brinch,
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */

#include "lime.h"

void
predefinedGrid(inputPars *par, struct grid *g){
  FILE *fp;
  int i,j,k;
  double x,y,z,scale,dummy,hfrac,h2frac,densum;
  const gsl_rng *ran = gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(ran,time(0));

  fp=fopen(par->pregrid,"r");
  par->ncell=par->pIntensity+par->sinkPoints;

  /* Reading in x, y, z, density, molabun, gastemp, dusttemp, dopb, vx, vy, vz */
  for(i=0;i<par->pIntensity;i++){
    fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &g[i].x[0], &g[i].x[1], &g[i].x[2], &densum, &g[i].abun[0], &g[i].t[0], &g[i].t[1], &g[i].dopb, &g[i].vel[0], &g[i].vel[1], &g[i].vel[2]);

    g[i].sink=0;
    g[i].id = i;
    //g[i].abun[0] = 1.e-8;
     
    /* Must specify how the total density is distributed amongst collision partners.
       Here, we have 3:1 ortho/para for H2
     */
    par->collPart = 2;
    g[i].dens[0] = densum*0.25;
    g[i].dens[1] = densum*0.75;

    //Next line redundant with molinit?
	g[i].nmol[0]=g[i].abun[0]*(g[i].dens[0]+g[i].dens[1]);
	

	/* Below comment and code from Brinch */
    /* This next step needs to be done, even though it looks stupid */
	g[i].dir=malloc(sizeof(point)*1);
	g[i].ds =malloc(sizeof(double)*1);
	g[i].neigh =malloc(sizeof(struct grid *)*1);
	if(!silent) progressbar((double) i/((double)par->pIntensity-1), 4);
  }

  for(i=par->pIntensity;i<par->ncell;i++){
	x=2*gsl_rng_uniform(ran)-1.;
	y=2*gsl_rng_uniform(ran)-1.;
	z=2*gsl_rng_uniform(ran)-1.;
	if(x*x+y*y+z*z<1){
  	  scale=par->radius*sqrt(1/(x*x+y*y+z*z));
	  g[i].id=i;
	  g[i].x[0]=scale*x;
	  g[i].x[1]=scale*y;
	  g[i].x[2]=scale*z;
	  g[i].sink=1;
	  g[i].dens[0]=1e-30;
	  g[i].t[0]=par->tcmb;
	  g[i].t[1]=par->tcmb;
	  g[i].dopb=0.;
    } else i--;
  }
  fclose(fp);


  qhull(par,g);
  distCalc(par,g);
//  getArea(par,g, ran);
//  getMass(par,g, ran);
//  getVelosplines_lin(par,g);
  getVelosplines(par,g);
  dumpGrid(par, g);
}
