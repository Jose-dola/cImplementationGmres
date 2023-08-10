/* gmres headers */
double* gmresUnlimited( int dim, double* b, void (*A_product)(double*,int,double*), double rc, double tol );
double* gmresM( int m, int dim, double* b, void (*A_product)(double*,int,double*), double rc, double tol, double** error );
