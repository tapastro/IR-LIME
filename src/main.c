/*
 *  main.c
 *  LIME, The versatile 3D line modeling environment
 *
 *  Created by Christian Brinch on 16/11/06.
 *  Copyright 2006-2017, Christian Brinch,
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 *  LIME is derived from RATRAN by Michiel Hogerheijde and Floris van der Tak,
 *  Copyright 2000, Hogerheijde and van der Tak, A&A, 362, 697, 2000.
 *
 *  DISCLAIMER: LIME is provided as is and the author does not resume any
 *  responsibility for erroneous results due to bugs in the code.
 *
 *  Any publication with results obtain using LIME should refer to
 *  Brinch & Hogerheijde, A&A, 523, A25, 2010
 *
 */

#include "lime.h"
/*Added header for radlite input */
#include "radlite.h"

int main () {
  int i;
  int initime=time(0);
  int popsdone=0;
  molData	  *m;
  inputPars	  par;
  struct grid *g;
  image		  *img;
  if(!silent) greetings();
  if(!silent) screenInfo();
  
  /*Added by TAP to read radlite grid */
  gridread(&ingridarray,&gridlength,&thtvec,&radvec);
  parseInput(&par,&img,&m);
  gridAlloc(&par,&g);

  if(par.pregrid) predefinedGrid(&par,g);
  else if(par.restart) popsin(&par,&g,&m,&popsdone);
  else buildGrid(&par,g);
  for(i = 0;i < par.nImages; i++){
    if(img[i].doline == 1 && popsdone == 0) {
      levelPops(m,&par,g,&popsdone);
    }
    if(img[i].doline == 0) {
      continuumSetup(i,img,m,&par,g);
    }
    raytrace(i,&par,g,m,img);
    writefits(i,&par,m,img);
  }
  if(!silent) goodnight(initime,img[0].filename);
  return 0;
}


double bilinInterpVal(double x1, double y1, double z1, int col) {
    double r1, t1, pz; 
    double wr1, wr2, wt1, wt2; 
    double d1, d2, val, tmp1, tmp2;
    int rind=-1,tind=-1,rmax=0,tmax=0;

    //Assume +/- z symmetry    
    pz = fabs(z1);
    //Find r, theta for interpolation
    r1 = sqrt(pow(x1,2)+pow(y1,2));
    t1 = atan(r1/pz);
    if(t1>PI/2.) t1 = PI - t1;
 
    for (int i = RLEN;i>-1;i--){
      if(radvec[i]>r1) rind=i; /*rind set to index of first rvec bigger than r1*/
    }
    if(rind<0) {
      //r1 is larger than all values in rvec
      rmax = 1;
      rind = RLEN-1;
      wr1 = 1.;
      wr2 = 0.;
    }
    else if(rind==0){
      //r1 is smaller than all values in rvec
      wr1 = 1.;
      wr2 = 0.;
    }
    else {
      //Find distances to nearest rvecs, then find weights based on distances
      d1 = r1-radvec[rind-1];
      d2 = radvec[rind]-r1;
      wr1 = d2/(d1+d2);
      wr2 = d1/(d1+d2);
    }

    for (int j = TLEN;j>-1;j--){
      if(thtvec[j]>t1) tind=j; /*tind set to index of first tvec bigger than t1*/
    }
    if(tind<0) {
      //t1 is larger than all values in tvec
      tmax = 1;
      tind = TLEN-1;
      wt1 = 1.;
      wt2 = 0.;
    }
    else if(tind==0) {
      //t1 is smaller than all values in tvec
      wt1 = 1.;
      wt2 = 0.;
    }
    else {
      //Find distances to nearest tvecs, then find weights based on distances
      d1 = t1 - thtvec[tind-1];
      d2 = thtvec[tind]-t1;
      wt1 = d2/(d1+d2);
      wt2 = d1/(d1+d2);
    }

    if((rind==0)||(rmax==1)) {
      if((tind==0)||(tmax==1)) {
        val = ingridarray[col][tind+(rind*TLEN)];
      }
      else {
        tmp1 = wt1 * ingridarray[col][tind-1+(rind*TLEN)];
        tmp2 = wt2 * ingridarray[col][tind+(rind*TLEN)];
        val = tmp1 + tmp2;
      }
    }
    else if((tind==0)||(tmax==1)) {
      tmp1 = wr1 * ingridarray[col][tind+((rind-1)*TLEN)];
      tmp2 = wr2 * ingridarray[col][tind+(rind*TLEN)];
      val = tmp1 + tmp2;
    }
    else {
      tmp1 = wr1 * ingridarray[col][tind-1+((rind-1)*TLEN)] + wr2 * ingridarray[col][tind-1+(rind*TLEN)];
      tmp2 = wr1 * ingridarray[col][tind+((rind-1)*TLEN)] + wr2 * ingridarray[col][tind+(rind*TLEN)];
      val = wt1 * tmp1 + wt2 * tmp2;
    } 
    if(val<0) {
      //fpb = fopen("error.log","a");
      //fprintf(fpb,"%lf %lf %lf %lf %lf %d %d \n", x1, y1, z1, r1, t1, rind, tind);
      //fprintf(fpb,"%lf %lf %lf %lf \n", wr1, wr2, wt1, wt2);
      //fprintf(fpb,"%lf %lf %lf \n \n", tmp1, tmp2, val);
      //fprintf(fpb,"%lf %lf \n", radvec[0], radvec[129]);
      //fclose(fpb);
      printf("Interp LESS THAN ZERO, col = %d",col); 
    }
    return val;
}
/* Added by TAP - functions for reading in RADMC grid */

void gridread(double*** array, int *glen, double** tvec, double** rvec) {

    FILE *fp,*fpr,*fpt;
    int i,j,k,gridlen=0;
    char string[280];
    double** arr;
    double* tv;
    double *rv;

    tv = malloc(TLEN*sizeof(double));
    rv = malloc(RLEN*sizeof(double));

    //fpr=fopen("radiuscopy_sr21.dat","r");
    fpr=fopen(RADFILE,"r");
    for(i=0;i<RLEN;i++){
      fscanf(fpr,"%lf\n", &rv[i]);
    }
    fclose(fpr);

    //fpt=fopen("thetacopy_sr21.dat","r");
    fpt=fopen(THTFILE,"r");
    for(j=0;j<TLEN;j++){
      fscanf(fpt,"%lf\n", &tv[j]);
    }
    fclose(fpt);

/*  MAKE THIS PART MORE SMOOTH - RIGHT NOW, ENTER GRID FILE HERE */
    fp=fopen(GRIDFILE, "r");
    
    while(fgets(string,280,fp) != NULL){
      gridlen++;
    }
    rewind(fp);

    arr = make2DDoubleArray(NCOL, gridlen);

    for(k=0;k<gridlen;k++){
      fscanf(fp,"%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n", &arr[0][k], &arr[1][k], &arr[2][k], &arr[3][k], &arr[4][k], &arr[5][k], &arr[6][k], &arr[7][k], &arr[8][k], &arr[9][k], &arr[10][k], &arr[11][k], &arr[12][k]);
    }
    fclose(fp);

    *rvec = rv;
    *tvec = tv;
    *glen = gridlen;
    *array = arr;
}

double** make2DDoubleArray(int arraySizeX, int arraySizeY) {
    double** theArray;
    theArray = malloc(arraySizeX*sizeof(double*));
    for (int i = 0; i < arraySizeX; i++)
    theArray[i] = malloc(arraySizeY*sizeof(double));
    return theArray;
}
