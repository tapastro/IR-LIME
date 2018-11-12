#define rlen 131
#define tlen 80
#define dencol 3

void gridread(double*** array, int* glen, double** thtvec, double** radvec);
double bilinInterpVal(double x, double y, double z, int col);
double** make2DDoubleArray(int arraySizeX, int arraySizeY);
int gridlength;
double** ingridarray;
double* thtvec;
double* radvec;

