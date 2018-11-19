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
  double x,y,z,scale,dummy;
  const gsl_rng *ran = gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(ran,time(0));

  fp=fopen(par->pregrid,"r");
  par->ncell=par->pIntensity+par->sinkPoints;

  //for(j=0;j<100000; j++) fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &k, &dummy, &dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy);
  for(i=0;i<par->pIntensity;i++){
//    fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &g[i].id, &g[i].x[0], &g[i].x[1], &g[i].x[2],  &g[i].dens[0], &g[i].t[0], &abun, &g[i].dopb, &g[i].vel[0], &g[i].vel[1], &g[i].vel[2]);
//    fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf\n", &g[i].id, &g[i].x[0], &g[i].x[1], &g[i].x[2],  &g[i].dens[0], &g[i].t[0], &abun, &g[i].dopb);
    fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &g[i].id, &g[i].x[0], &g[i].x[1], &g[i].x[2],  &g[i].dens[0], &g[i].t[0], &g[i].abun[0], &g[i].vel[0], &g[i].vel[1], &g[i].vel[2]);



    g[i].dopb=200;
    g[i].abun[0]=1e-9;

    g[i].x[0]=(g[i].x[0]+ (g[i].x[0]*0.1)*(2*gsl_rng_uniform(ran)-1.))/100.;
    g[i].x[1]=(g[i].x[1]+ (g[i].x[1]*0.1)*(2*gsl_rng_uniform(ran)-1.))/100.;
    g[i].x[2]=(g[i].x[2]+ (g[i].x[2]*0.1)*(2*gsl_rng_uniform(ran)-1.))/100.;


    g[i].sink=0;
    g[i].t[0]=50;
    g[i].dens[0]=1e12;
    g[i].t[1]=g[i].t[0];
	  g[i].nmol[0]=g[i].abun[0]*g[i].dens[0];
	  velocity(g[i].x[0],g[i].x[1],g[i].x[2],g[i].vel);

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
  if(par->gridfile) write_VTK_unstructured_Points(par, g);
}
