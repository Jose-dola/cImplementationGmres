/* headers of arnoldi's functions */
void arnoldi( double* b, int dim, int n, double** H, double** Q, double* (*A_product)(double*,int) );
void arnoldiRefined( double* b, int dim, int n, double** H, double** Q, double rc, double* (*A_product)(double*,int) );
void givensTest( double* b, int dim, int n, double** H, double** Q, double** R,
                 double* sin, double* cos, double rc, double* (*A_product)(double*,int) );
double** Qalloc(int n);
double** Halloc(int n);
void free_arnoldi(int n, double** Hn, double** Q);
