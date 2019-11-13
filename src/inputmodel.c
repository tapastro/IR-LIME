#include "lime.h"
#include "radlite.h"

void __attribute__((weak))
gridsize(int *rlen,int *tlen){
  if(!silent) bail_out("Inputmodel is on, but size of input grid not specified in gridsize!");
  exit(1);
}

void __attribute__((weak))
filenames(char *rname,char *tname,char *gname){
  if(!silent) bail_out("Inputmodel is on, but filenames for input rad,theta,grid files not specified!");
  exit(1);
}

double bilinInterpVal(double x1, double y1, double z1, int col) {
    double rpol,r1, t1, pz; 
    double wr1, wr2, wt1, wt2; 
    double d1, d2, val, tmp1, tmp2;
    int rind=-1,tind=-1,rmax=0,tmax=0;
    int RLEN=0,TLEN=0;

    gridsize(&RLEN,&TLEN);

    //Assume +/- z symmetry    
    pz = fabs(z1);
    //Find r, theta for interpolation
    rpol = sqrt(pow(x1,2)+pow(y1,2));
    t1 = atan(rpol/pz);
    r1 = sqrt(pow(x1,2)+pow(y1,2)+pow(z1,2));
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
    char rname[200],tname[200],gname[200];
    int RLEN=0,TLEN=0;

    gridsize(&RLEN,&TLEN);

    tv = malloc(TLEN*sizeof(double));
    rv = malloc(RLEN*sizeof(double));

    //Give 1D vector of radii in grid, same for theta, finally the grid itself
    //Grid: x, y, z, dens, abn, Tg, Td, Vdisp, Vx, Vy, Vz, Hfrac, H2frac
    filenames(rname,tname,gname);

    //fpr=fopen("radiuscopy_sr21.dat","r");
    fpr=fopen(rname,"r");
    for(i=0;i<RLEN;i++){
      fscanf(fpr,"%lf\n", &rv[i]);
    }
    fclose(fpr);

    //fpt=fopen("thetacopy_sr21.dat","r");
    fpt=fopen(tname,"r");
    for(j=0;j<TLEN;j++){
      fscanf(fpt,"%lf\n", &tv[j]);
    }
    fclose(fpt);

    fp=fopen(gname, "r");
    
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
