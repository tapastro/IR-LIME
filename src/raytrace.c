/*
 *  raytrace.c
 *  LIME, The versatile 3D line modeling environment
 *
 *  Created by Christian Brinch on 16/12/06.
 *  Copyright 2006-2017, Christian Brinch,
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */


typedef struct {double x,y, *intensity, *tau;} data;

#ifdef QHULL_INC_QHULL
#include <qhull/qhull_a.h>
#else
#include <libqhull/qhull_a.h>
#endif
#include "lime.h"

void
velocityspline2(double x[3], double dx[3], double ds, double binv, double deltav, double *vfac){
  int i,steps=10;
  double v,d,val,vel[3];
  *vfac=0.;
  for(i=0;i<steps;i++){
    d=i*ds/steps;
    velocity(x[0]+(dx[0]*d),x[1]+(dx[1]*d),x[2]+(dx[2]*d),vel);
    v=deltav-veloproject(dx,vel);
    val=fabs(v)*binv;
    if(val <=  2500.){
      *vfac+= exp(-(val*val));
    }
  }
  *vfac=*vfac/steps;
  return;
}

void
velocityspline3(double x[3], double dx[3], double ds, double binv, double vel[3], double deltav, double *vfac){
  int i,steps=10;
  double v,d,val;
  *vfac=0.;
  for(i=0;i<steps;i++){
    d=i*ds/steps;
    //velocity(x[0]+(dx[0]*d),x[1]+(dx[1]*d),x[2]+(dx[2]*d),vel);
    v=deltav-veloproject(dx,vel);
    val=fabs(v)*binv;
    if(val <=  2500.){
      *vfac+= exp(-(val*val));
    }
  }
  *vfac=*vfac/steps;
  return;
}


void
line_plane_intersect(struct grid *g, double *ds, int posn, int *nposn, double *dx, double *x){
  double newdist, numerator, denominator ;
  int i;

  for(i=0;i<g[posn].numNeigh;i++) {
    /* Find the shortest distance between (x,y,z) and any of the posn Voronoi faces */
    /* ds=(p0-l0) dot n / l dot n */

    numerator=((g[posn].x[0]+g[posn].dir[i].x[0]/2. - x[0]) * g[posn].dir[i].x[0]+
               (g[posn].x[1]+g[posn].dir[i].x[1]/2. - x[1]) * g[posn].dir[i].x[1]+
               (g[posn].x[2]+g[posn].dir[i].x[2]/2. - x[2]) * g[posn].dir[i].x[2]);

    denominator=(dx[0]*g[posn].dir[i].x[0]+dx[1]*g[posn].dir[i].x[1]+dx[2]*g[posn].dir[i].x[2]);

    if(fabs(denominator) > 0){
      newdist=numerator/denominator;
      if(newdist<*ds && newdist > 1e4){
        *ds=newdist;
        *nposn=g[posn].neigh[i]->id;
      }
	}
  }
}


void
tracerays(int px, data *ray, int tmptrans, int im, inputPars *par, struct grid *g, molData *m, image *img, int nlinetot, int *counta, int *countb){
  double *tau, *subintens, shift, deltav;
  double xp,yp,zp, x[3], dx[3],dist,ndist,col,ds,dtau,jnu,alpha,snu,snu_pol[3],vfac=0.;
  int ichan, i, posn, nposn, iline;

  tau = malloc(sizeof(double)*img[im].nchan);
  subintens = malloc(sizeof(double)*img[im].nchan);
  for(ichan=0;ichan<img[im].nchan;ichan++){
    tau[ichan]=0.0;
    subintens[ichan]=0.0;
  }


  xp=ray[px].x;
  yp=ray[px].y;


  /* Rotation matrix

          |1	  0		    0   |
   R_x(a)=|0	cos(a)	sin(a)	|
          |0 -sin(a)	cos(a)	|

          |cos(b)	0	-sin(b) |
   R_y(b)=|  0		1	   0	|
          |sin(b)	0	 cos(b) |

          |cos(b)  		    0 	   sin(b)   |
   Rot =  |sin(a)sin(b)	cos(a)	sin(a)cos(b)|
          |cos(a)sin(b)   -sin(a)  cos(a)cos(b)|

   */
  if(sqrt(xp*xp+yp*yp)/par->radius <= 1 ) {
    zp=par->radius*cos(asin(sqrt(xp*xp+yp*yp)/par->radius));

    x[0]=xp*cos(img[im].phi)                   +yp*0.                -zp*sin(img[im].phi);
    x[1]=xp*sin(img[im].theta)*sin(img[im].phi)+yp*cos(img[im].theta)+zp*sin(img[im].theta)*cos(img[im].phi);
    x[2]=xp*cos(img[im].theta)*sin(img[im].phi)-yp*sin(img[im].theta)+zp*cos(img[im].theta)*cos(img[im].phi);

    dx[0]= sin(img[im].phi);
    dx[1]=-sin(img[im].theta)*cos(img[im].phi);
    dx[2]=-cos(img[im].theta)*cos(img[im].phi);

    dist=1e60;
    posn=-1;
    for(i=0;i<par->ncell;i++){
      ndist=sqrt(pow(x[0]-g[i].x[0],2)+pow(x[1]-g[i].x[1],2)+pow(x[2]-g[i].x[2],2));
      if(ndist<dist){
        posn=i;
        dist=ndist;
      }
    }

    col=0;
    do{
      ds=2.*zp-col;
      nposn=-1;
      line_plane_intersect(g,&ds,posn,&nposn,dx,x);
//      if(nposn==-1){
//        if(!silent) bail_out("Error: Ray-tracer failed to detect Voronoi edge");
 //       exit(1);
 //     }
      if(par->polarization){
        for(ichan=0;ichan<img[im].nchan;ichan++){
          sourceFunc_pol(snu_pol,&dtau,ds,m,vfac,g,posn,0,0,img[im].theta);
          subintens[ichan]+=exp(-tau[ichan])*(1.-exp(-dtau))*snu_pol[ichan];
          tau[ichan]+=dtau;
        }
      } else {
        for(ichan=0;ichan<img[im].nchan;ichan++){
          jnu=.0;
          alpha=0.;
          snu=0.;
          dtau=0.;

          for(iline=0;iline<nlinetot;iline++){
            if(img[im].doline && m[counta[iline]].freq[countb[iline]] > img[im].freq-img[im].bandwidth/2. && m[counta[iline]].freq[countb[iline]] < img[im].freq+img[im].bandwidth/2.){
              if(img[im].trans > -1){
                shift=(m[counta[iline]].freq[countb[iline]]-m[counta[iline]].freq[img[im].trans])/m[counta[iline]].freq[img[im].trans]*CLIGHT;
              } else {
                shift=(m[counta[iline]].freq[countb[iline]]-img[im].freq)/img[im].freq*CLIGHT;
              }
              deltav=(ichan-(int)(img[im].nchan/2.))*img[im].velres-img[im].source_vel + shift;

              if(!par->pregrid && !par->discreteVelocities) velocityspline2(x,dx,ds,g[posn].mol[counta[iline]].binv,deltav,&vfac);
              else vfac=gaussline(deltav-veloproject(dx,g[posn].vel),g[posn].mol[counta[iline]].binv);
//              else velocityspline3(x,dx,ds,g[posn].mol[counta[iline]].binv,g[posn].vel,deltav,&vfac);

              sourceFunc_line(&jnu,&alpha,m,vfac,g,posn,counta[iline],countb[iline]);
            }
          }

          if(img[im].doline && img[im].trans > -1) sourceFunc_cont(&jnu,&alpha,g,posn,0,img[im].trans);
          else if(img[im].doline && img[im].trans == -1) sourceFunc_cont(&jnu,&alpha,g,posn,0,tmptrans);
          else sourceFunc_cont(&jnu,&alpha,g,posn,0,0);
          if(fabs(alpha)>0.){
            snu=(jnu/alpha)*m[0].norminv;
            dtau=alpha*ds;
          }
          subintens[ichan]+=exp(-tau[ichan])*(1.-exp(-dtau))*snu;
          tau[ichan]+=dtau;

        }
      }

      /* new coordinates */
      for(i=0;i<3;i++) x[i]+=ds*dx[i];
      col+=ds;
      posn=nposn;
    } while(col < 2*zp);

    /* add or subtract cmb */
    if(!par->polarization){
      for(ichan=0;ichan<img[im].nchan;ichan++){
        subintens[ichan]+=(exp(-tau[ichan])-1.)*m[0].local_cmb[tmptrans];
      }
    } else {
        subintens[0]+=(exp(-tau[0])-1.)*m[0].local_cmb[tmptrans];
    }

  }
#pragma omp critical
  {
  for(ichan=0;ichan<img[im].nchan;ichan++){
    ray[px].intensity[ichan]=subintens[ichan];
    ray[px].tau[ichan]=tau[ichan];
  }
  }
  free(tau);
  free(subintens);
}



void
raytrace(int im, inputPars *par, struct grid *g, molData *m, image *img){
  int *counta, *countb,nlinetot;
  int ichan,i,px,iline,tmptrans,count;
  double size,xp,yp,minfreq;
  const gsl_rng *ran = gsl_rng_alloc(gsl_rng_ranlxs2);	/* Random number generator */
  gsl_rng_set(ran,time(0));
  data *ray;

  int sg,n;
  double cx,cy;

  double x1,x2,x3,y1,y2,y3,z1,z2,z3,xt[3],yt[3],di,p,d1,d2,d3;
  int zt[3];
  int c;

  char flags[255];
  boolT ismalloc = False;
  facetT *facet, *neighbor, **neighborp;;
  vertexT *vertex,**vertexp;
  coordT *pt_array;

  int id;
  coordT point[3];
  boolT isoutside;
  realT bestdist;

  /* Determine whether there are blended lines or not */
  if(img[im].doline==0) nlinetot=1;
  else lineCount(par->nSpecies, m, &counta, &countb, &nlinetot);

  /* Fix the image parameters */
  if(img[im].freq < 0) img[im].freq=m[0].freq[img[im].trans];
  if(img[im].nchan == 0 && img[im].bandwidth>0){
    img[im].nchan=(int) (img[im].bandwidth/(img[im].velres/CLIGHT*img[im].freq));
  } else if (img[im].velres<0 && img[im].bandwidth>0){
    img[im].velres = img[im].bandwidth*CLIGHT/img[im].freq/img[im].nchan;
  } else img[im].bandwidth = img[im].nchan*img[im].velres/CLIGHT * img[im].freq;

  if(img[im].trans<0){
    minfreq=1e30;
    tmptrans=-1;
    for(iline=0;iline<m[0].nline;iline++){
      if(fabs(img[im].freq-m[0].freq[iline])<minfreq){
        minfreq=fabs(img[im].freq-m[0].freq[iline]);
        tmptrans=iline;
      }
    }
  } else tmptrans=img[im].trans;

  /* Allocate dynamical arrays */
  ray= malloc(sizeof(data) * (par->pIntensity));

  for(i=0;i<par->pIntensity;i++){
    ray[i].x=g[i].x[0];
    ray[i].y=g[i].x[1];
    ray[i].intensity=malloc(sizeof(double) * img[im].nchan);
    ray[i].tau=malloc(sizeof(double) * img[im].nchan);
    for(ichan=0;ichan<img[im].nchan;ichan++) {
      ray[i].intensity[ichan]=0.0;
      ray[i].tau[ichan]=0.0;
    }
  }


  /* Smooth out the distribution of rays */
  for(sg=0;sg<20;sg++){
    pt_array=malloc(2*sizeof(coordT)*par->pIntensity);

    for(i=0;i<par->pIntensity;i++) {
      pt_array[i*2+0]=ray[i].x;
      pt_array[i*2+1]=ray[i].y;
    }

    sprintf(flags,"qhull v s Qbb T0");
    if (!qh_new_qhull(2, par->pIntensity, pt_array, ismalloc, flags, NULL, NULL)) {

      qh_setvoronoi_all();

      FORALLvertices {
        i=qh_pointid(vertex->point);

        cx=0.;
        cy=0.;
        n=0;
        FOREACHneighbor_(vertex) {
          if (!neighbor->upperdelaunay) n++;
        }
        if(n>0){


        } else {
          if(!silent) bail_out("Qhull error");
          exit(0);
        }

        FOREACHneighbor_(vertex) {
          if (!neighbor->upperdelaunay) {
            cx+=neighbor->center[0];
            cy+=neighbor->center[1];
          }
        }

        ray[i].x = ray[i].x - (ray[i].x-cx/ (double) n)*0.1;
        ray[i].y = ray[i].y - (ray[i].y-cy/ (double) n)*0.1;
      }
    } else {
      printf("qhull error\n");
    }

    qh_freeqhull(!qh_ALL);
    free(pt_array);
  }



  /* Ray-trace in parallel */
  count=0;
  omp_set_dynamic(0);
  #pragma omp parallel private(px) num_threads(par->threads)
  {
    #pragma omp for
    /* Main loop through rays */
    for(px=0;px<par->pIntensity;px++){

      tracerays(px, ray, tmptrans, im, par, g, m, img, nlinetot, counta, countb);

      #pragma omp critical
      {
        ++count;
      }

      if (omp_get_thread_num() == 0){
        if(!silent) progressbar((double)(count)/(double)(par->pIntensity-1), 13);
      }
    }
  }


  /* Remap rays onto pixel grid */
  pt_array=malloc(2*sizeof(coordT)*par->pIntensity);

  for(i=0;i<par->pIntensity;i++) {
    pt_array[i*2+0]=ray[i].x;
    pt_array[i*2+1]=ray[i].y;
  }

/* This allocation belongs to "Shepard's method" below
  d=malloc(sizeof(double)*par->pIntensity);
*/
  size=img[im].distance*img[im].imgres;

  sprintf(flags,"qhull d Qbb");
  if (!qh_new_qhull(2, par->pIntensity, pt_array, ismalloc, flags, NULL, NULL)) {


    for(px=0;px<img[im].pxls*img[im].pxls;px++){
      for(ichan=0;ichan<img[im].nchan;ichan++){
        img[im].pixel[px].intense[ichan]=0.0;
        img[im].pixel[px].tau[ichan]=0.0;
      }
      xp=size*(0.5+px%img[im].pxls)-size*img[im].pxls/2.;
      yp=size*(0.5+px/img[im].pxls)-size*img[im].pxls/2.;

/*
This part works great! This is "Shepard's method" with a weight of 8. Slow, unfortunately.
      for(ichan=0;ichan<img[im].nchan;ichan++){
        img[im].pixel[px].intense[ichan] = 0.;
        di=0;
        for(i=0;i<par->pIntensity;i++){
          d[i]=1./pow(sqrt(pow(xp-ray[i].x,2)+ pow(yp-ray[i].y,2)),8.);
          img[im].pixel[px].intense[ichan] += ray[i].intensity[ichan]*d[i];
          di+=d[i];
        }
        img[im].pixel[px].intense[ichan] /= di;
      }
*/


      point[0]=xp;
      point[1]=yp;
      point[2]=0.;

      qh_setdelaunay (3, 1, point);
      facet= qh_findbestfacet (point, qh_ALL, &bestdist, &isoutside);

      c=0;
      FOREACHvertex_( facet->vertices ) {
        id=qh_pointid(vertex->point);
        xt[c]=ray[id].x; yt[c]=ray[id].y; zt[c]=id;
        c++;
      }


      x1=xt[0];x2=xt[1];x3=xt[2];
      y1=yt[0];y2=yt[1];y3=yt[2];

      for(ichan=0;ichan<img[im].nchan;ichan++){
        z1=ray[zt[2]].intensity[ichan];z2=ray[zt[1]].intensity[ichan];z3=ray[zt[2]].intensity[ichan];


        p=1.;
        d1=1./pow(sqrt(pow(xp-x1,2)+ pow(yp-y1,2)),p);
        d2=1./pow(sqrt(pow(xp-x2,2)+ pow(yp-y2,2)),p);
        d3=1./pow(sqrt(pow(xp-x3,2)+ pow(yp-y3,2)),p);
        di=d1+d2+d3;
        img[im].pixel[px].intense[ichan] = 1./di * (z1*d1 + z2*d2 + z3*d3);
        
	z1=ray[zt[2]].tau[ichan];z2=ray[zt[1]].tau[ichan];z3=ray[zt[2]].tau[ichan];
	img[im].pixel[px].tau[ichan] = 1./di * (z1*d1 + z2*d2 + z3*d3);
      }


    }
  } else {
	if(!silent) bail_out("Qhull failed to triangulate");
	exit(1);
  }



  img[im].trans=tmptrans;
//  free(counta);
//  free(countb);
}
