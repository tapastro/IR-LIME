/*
 *  blowpops.c
 *  LIME, The versatile 3D line modeling environment
 *
 *  Created by Christian Brinch on 14/11/07.
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

void
popsout(inputPars *par, struct grid *g, molData *m){
  FILE *fp;
  int i,j,k,l,best=-1;
  double dens, mass=0.,dist,*vol;
  double dp,dpbest,*farea,suma;
  char flags[255];
  boolT ismalloc = False;
  facetT *facet, *neighbor, **neighborp;
  vertexT *vertex;
  coordT *pt_array;
  typedef struct {coordT *pt_array;int vps;int *flag;} S;
  S *pts;
  int curlong, totlong;

  gsl_rng *ran = gsl_rng_alloc(gsl_rng_ranlxs2);    /* Random number generator */
  gsl_rng_set(ran,time(0));


  /* int i,mi,c,q=0,best; */
  /* double vel[3],ra[100],rb[100],za[100],zb[100],min; */
  vol=malloc(sizeof(double)*par->pIntensity);
  pts=malloc(sizeof(S)*par->pIntensity);
  for(i=0;i<par->pIntensity;i++){
    pts[i].vps=0;
  }
  pt_array=malloc(DIM*sizeof(coordT)*par->ncell);
  for(i=0;i<par->ncell;i++) {
    for(j=0;j<DIM;j++) {
      pt_array[i*DIM+j]=g[i].x[j];
    }
  }

  sprintf(flags,"qhull v Qbb");
  if (!qh_new_qhull(DIM, par->ncell, pt_array, ismalloc, flags, NULL, NULL)) {
    qh_setvoronoi_all();
    FORALLvertices {
      i=qh_pointid(vertex->point);
      if(i<par->pIntensity){
        pts[i].vps=0;
        FOREACHneighbor_(vertex) {
          if (!neighbor->upperdelaunay) pts[i].vps++;
        }
        if(pts[i].vps > 0){
          pts[i].pt_array=malloc(DIM*sizeof(coordT)*pts[i].vps);
          pts[i].flag=malloc(DIM*sizeof(int)*pts[i].vps);
        } else {
          if(!silent) bail_out("Qhull error");
          exit(0);
        }
        k=0;
        FOREACHneighbor_(vertex) {
          if (!neighbor->upperdelaunay) {
            for(j=0;j<DIM;j++) pts[i].pt_array[k*DIM+j]=neighbor->center[j]+(gsl_rng_uniform(ran)*2-1)*par->radius/1e6;
            k+=1;
          }
        }
      }
    }
  } else {
    if(!silent) bail_out("Qhull error");
    exit(0);
  }
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort (&curlong, &totlong);

  sprintf(flags,"qhull ");
  for(i=0;i<par->pIntensity;i++){
    g[i].w=malloc(sizeof(double)*g[i].numNeigh);
    vol[i]=0.;
    suma=0.;
    if(pts[i].vps>0){
      farea=malloc(sizeof(double)*pts[i].vps);
    } else {
      if(!silent) bail_out("Qhull error");
      exit(0);
    }
    if (!qh_new_qhull(DIM, pts[i].vps, pts[i].pt_array, ismalloc, flags, NULL, NULL)) {
      FORALLfacets {
        dpbest=0.;
        for(j=0;j<g[i].numNeigh;j++){
          dp=facet->normal[0]*g[i].dir[j].xn[0]+facet->normal[1]*g[i].dir[j].xn[1]+
          facet->normal[2]*g[i].dir[j].xn[2];
          if(fabs(dp)>dpbest){
            dpbest=fabs(dp);
            best=j;
          }
        }
        if (!facet->normal)
          continue;
        if (facet->upperdelaunay && qh ATinfinity)
          continue;
        farea[best]=qh_facetarea (facet);
        suma+=farea[best];
        if (!qh DELAUNAY) {
          qh_distplane (qh interior_point, facet, &dist);
          vol[i] += -dist * farea[best]/ qh hull_dim;
        }
      }
      for(j=0;j<g[i].numNeigh;j++) g[i].w[j]=farea[j]/suma; //if(g[i].w[j]<1e-2) g[i].w[j]=1e-2;
      free(pts[i].flag);
      free(pts[i].pt_array);
      mass+=vol[i]*g[i].dens[0];
    }
    free(farea);
  }

  if((fp=fopen(par->outputfile, "w"))==NULL){
    if(!silent) bail_out("Error writing output populations file!");
    exit(1);
  }
  fprintf(fp,"# Column definition: x, y, z, H2 density, kinetic gas temperature, molecular abundance, convergence flag, pops_0...pops_n\n");
  for(j=0;j<par->ncell-par->sinkPoints;j++){
    dens=0.;
    for(l=0;l<par->collPart;l++) dens+=g[j].dens[l];
    //fprintf(fp,"%e %e %e %e %e %e %d ", g[j].x[0], g[j].x[1], g[j].x[2], dens, g[j].t[0], g[j].nmol[0]/dens, g[j].conv);
    fprintf(fp,"%e %e %e %e %e %e %e %e %e %e %e %d ", g[j].x[0], g[j].x[1], g[j].x[2], g[j].dens[0], g[j].dens[1], g[j].dens[2], g[j].t[0], g[j].t[1], g[j].nmol[0], g[j].mol[0].dust[217], vol[j], g[j].conv);
    for(k=0;k<m[0].nlev;k++) fprintf(fp,"%e ",g[j].mol[0].pops[k]);
    fprintf(fp,"\n");
  }
  fclose(fp);
}


void
binpopsout(inputPars *par, struct grid *g, molData *m){
  FILE *fp;
  int i,j;

  if((fp=fopen(par->binoutputfile, "wb"))==NULL){
    if(!silent) bail_out("Error writing binary output populations file!");
    exit(1);
  }

  fwrite(&par->radius,   sizeof(double), 1, fp);
  fwrite(&par->ncell,    sizeof(int), 1, fp);
  fwrite(&par->nSpecies, sizeof(int), 1, fp);

  for(i=0;i<par->nSpecies;i++){
    fwrite(&m[i].nlev,  sizeof(int),               1,fp);
    fwrite(&m[i].nline, sizeof(int),               1,fp);
    fwrite(&m[i].npart, sizeof(int),               1,fp);
    fwrite(m[i].ntrans, sizeof(int)*m[i].npart,    1,fp);
    fwrite(m[i].lal,    sizeof(int)*m[i].nline,    1,fp);
    fwrite(m[i].lau,    sizeof(int)*m[i].nline,    1,fp);
    fwrite(m[i].aeinst, sizeof(double)*m[i].nline, 1,fp);
    fwrite(m[i].freq,   sizeof(double)*m[i].nline, 1,fp);
    fwrite(m[i].beinstl,sizeof(double)*m[i].nline, 1,fp);
    fwrite(m[i].beinstu,sizeof(double)*m[i].nline, 1,fp);
    fwrite(m[i].local_cmb, sizeof(double)*m[i].nline,1,fp);
    fwrite(&m[i].norm,  sizeof(double),1,fp);
    fwrite(&m[i].norminv,sizeof(double),1,fp);
  }

  for(i=0;i<par->ncell;i++){
    fwrite(&g[i].id,   sizeof(int),      1, fp);
    fwrite(&g[i].x,  3*sizeof(double),   1, fp);
    fwrite(&g[i].vel,3*sizeof(double),   1, fp);
    fwrite(&g[i].sink, sizeof(int),      1, fp);
    fwrite(g[i].nmol,  sizeof(double)*par->nSpecies,1, fp);
    fwrite(&g[i].dopb, sizeof g[i].dopb, 1, fp);
    for(j=0;j<par->nSpecies;j++){
      fwrite(g[i].mol[j].pops,  sizeof(double)*m[j].nlev, 1, fp);
      fwrite(g[i].mol[j].knu,   sizeof(double)*m[j].nline,1, fp);
      fwrite(g[i].mol[j].dust,  sizeof(double)*m[j].nline,1, fp);
      fwrite(&g[i].mol[j].dopb, sizeof(double),           1, fp);
      fwrite(&g[i].mol[j].binv, sizeof(double),           1, fp);
    }
    fwrite(&g[i].dens[0], sizeof(double), 1, fp);
    fwrite(&g[i].t[0],    sizeof(double), 1, fp);
    fwrite(&g[i].abun[0], sizeof(double), 1, fp);
 }


  fclose(fp);

}
