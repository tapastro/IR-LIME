#define RLEN 131
#define RADFILE "radius.dat"

#define TLEN 80
#define THTFILE "theta.dat"

#define GRIDFILE "radlite_griddata_posval.dat"
#define NCOL 13

#define DENCOL 3
#define ABNCOL 4
#define TEMCOL 5
#define DTEMCOL 6
#define DOPCOL 7
#define HCOL 11
#define H2COL 12

void gridread(double*** array, int* glen, double** thtvec, double** radvec);
double bilinInterpVal(double x, double y, double z, int col);
double** make2DDoubleArray(int arraySizeX, int arraySizeY);
int gridlength;
double** ingridarray;
double* thtvec;
double* radvec;

