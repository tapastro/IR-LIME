/*
 *  weights.c
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

#include "lime.h"
#include "radlite.h"
int
pointEvaluation(inputPars *par,double ran, double x, double y, double z){
  FILE *fp;
  double weight1, weight2, weight3, r, val1[9],n0,densgrad0,tem0,tgrad0,val2[9],val3[9];
  double sphr,az,cut,powr,varpowr;
  az = fabs(z);
  r = sqrt(pow(x,2)+pow(y,2));
  sphr = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
  cut = 2.7; //Reasonable r/z cut for young RNo90 disk
  //cut = 6.0; //Reasonable r/z cut for SR 21
  //cut = 8.7; //Reasonable r/z cut for first flared_test grid
  powr = 2.0;//was 2.5 for RNo90 disk grid, 1.8 for flared
  varpowr = 1.0+1.5*(fabs(r-6.5*AU)/93.5*AU);
  if (r/az < cut){
  weight1 = pow((r/az)/cut,varpowr);
  } else {
  weight1 = pow((r/az)/cut,-1.*powr); 
  }

  n0=4e20; // Appropriate max dens for RNo90
  //n0=8.0e17; // Appropriate max dens for SR 21
  density(x,y,z,val1);
  weight2=pow(val1[0]/n0,0.5);

  //weight3 = pow(1./((r-4.7*AU)/AU),5.0)*(r/az+az/r);
  //weight3 = 0.0;
  //if(sphr<6.4*AU){
  //  weight3 = 0.3;
  //}

  
  //if(ran < weight1 || ran < weight2 || ran < weight3){
  if(ran < weight1 || ran < weight2){
  //fp = fopen("weights.txt","a");
  //fprintf(fp,"%e %e   %e   %e   %e\n",r,az,weight1,weight2,weight3);
  //fprintf(fp,"%e  %e\n",weight1,weight2);
  //fclose(fp);
  return 1;
  }
  else return 0;
}
