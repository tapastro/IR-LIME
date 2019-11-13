#define NCOL 13

void gridread(double*** array, int* glen, double** thtvec, double** radvec);
double bilinInterpVal(double x, double y, double z, int col);
double** make2DDoubleArray(int arraySizeX, int arraySizeY);
int gridlength;
double** ingridarray;
double* thtvec;
double* radvec;
void filenames(char* rname, char* tname, char* gname);
void gridsize(int* RLEN,int* TLEN);
void maxdens(double* mdens);
void minscale(double* minr);
