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
  double weight1, weight2, weight3, r, val1[9],n0,val2[9],val3[9];
  double polr,sphr,az,cut,powr,varpowr,minr,mingridr;
  
  az = fabs(z);
  r = sqrt(pow(x,2)+pow(y,2));
  polr = sqrt(pow(x,2)+pow(y,2));
  sphr = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
  
  cut = (1./11.)*pow((polr/(AU)),1.5/7.);
  powr = 6.0;//was 2.5 for RNo90 disk grid, 1.8 for flared
  
  varpowr = 1.0+1.5*(fabs(r-6.5*AU)/93.5*AU);
  
  if (az/r < cut){
  weight1 = pow((az/r)/cut,powr);
  } else {
  weight1 = pow((az/r)/cut,-1.*powr); 
  }

  maxdens(&n0); // Appropriate max dens for SR 21
  density(x,y,z,val1);
  weight2=pow(val1[0]/n0,1.8);


  minscale(&minr);
  //minr=0.01*AU;
  mingridr=log10(minr/0.95/AU);
  minr=log10(minr/AU);
  weight3 = 0.2*pow((mingridr-minr)/(log10(sphr/AU)-minr),4.0);//*(r/az+az/r);
  //weight3 = 0.0;
  //if(sphr<6.4*AU){
  //  weight3 = 0.3;
  //}

  
  if(ran < weight1 || ran < weight2 || ran < weight3){
  //if(ran < weight1 || ran < weight3){
  //fp = fopen("weights.txt","a");
  //fprintf(fp,"%e %e   %e   %e   %e\n",r,az,weight1,weight2,weight3);
  //fclose(fp);
  return 1;
  }
  else return 0;
}
