/*
 * Example of gmres function usage
 *
 * by Jose Luis Dorado
 * University of Barcelona
 * March 2018
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gmres.h"

/* A_product header */
void A_product(double* v, int dim, double* Av);
/* global variables */
double** A; /* matrix A */

int main(void)
{
  int      i;       /* auxiliary variable */
  int      j;       /* auxiliary variable */
/*  double   aux;*/     /* auxiliary variable */
  FILE*    f;       /* pointer to file */
  double*  b;       /* starting vector */
  int      dim;     /* dimension of vector b */
  double*  x;       /* solution */
/*  double*  check;*/   /* auxiliary vector to check the solution */
  double*  err[1];  /* error vector */


  /* open the file with the vector b */
  f = fopen("b.txt","r");
  if (f == NULL) { puts("problem to open 'b.txt' file"); exit(1); }
  /* get the dimension of the vector b */
  fscanf(f,"%d",&dim); 
  /* get the vector b */
  b = (double*)malloc(dim*sizeof(double));
  if (b == NULL) { puts("memory problem"); exit(1); }
  for (i=0; i<dim; i++) fscanf(f,"%lf",b+i);
  fclose(f);
 
  /* open the file with the matrix A */
  f = fopen("A.txt","r");
  if (f == NULL) { puts("problem to open 'A.txt' file"); exit(1); }
  /* get the matrix A */
  A = (double**)malloc(dim*sizeof(double*));
  if (A == NULL) { puts("memory problem"); exit(1); }
  for (i=0; i<dim; i++)
  {
    A[i] = (double*)malloc(dim*sizeof(double));
    if (A[i] == NULL) { puts("memory problem"); exit(1); }
    for (j=0; j<dim; j++) fscanf(f,"%lf",&A[i][j]); 
  }
  fclose(f);

/* allocate x */
x = (double*)malloc(dim*sizeof(double));
  
  /* show A */
/*  puts("A:");
  for (i=0; i<dim; i++)
  {
    for (j=0; j<dim; j++) printf("% 10lf ",A[i][j]);
    puts(""); 
  }
  puts("");
*/  /* show b */
/*  puts("b:");
  for (i=0; i<dim; i++) printf("% 10lf ", b[i]);
  puts("");
*/
  /* GMRES Unlimited to solve the system */
/*  puts("doing Unlimited GMRES ...");
  x = gmresUnlimited(dim,b,A_product,0.25,1e-15);
*/
  /* show the solution */
/*  if (x == NULL) puts("no solution found for the system Ax=b");
  else
  {
    puts("solution found of Ax=b; x= "); 
    for (i=0; i<dim; i++) printf("% 10le ", x[i]);
    puts("");
  }
*/
  /* check */
/*  puts("checking the solution of Unlimited GMRES...");
  check = (double*)malloc(dim*sizeof(double));
  if (check == NULL) { puts("memory problem"); exit(1); } 
  A_product(x,dim,check);
  printf("error_check= ");
  aux=0;
  for (i=0; i<dim; i++) aux+=(b[i]-check[i])*(b[i]-check[i]);
  printf("% 10le\n",sqrt(aux));

  puts("");
*/
  /* GMRESm to solve the system */
/*  puts("doing Restarted GMRES ...");
*/  if ( gmresM(1,dim,b,A_product,0.25,1e-15,err,x) < 0 ) printf("GMRESm stagnates\n");

  /* show the solution */
/*  if (x == NULL) puts("no solution found for the system Ax=b");
  else
  {
    puts("solution found of Ax=b; x= "); 
    for (i=0; i<dim; i++) printf("% 10le ", x[i]);
    puts("");
  }
*/  
  /* show the error */
/*  printf("error_found= ");
  aux=0;
  for (i=0; i<dim; i++) aux+=err[0][i]*err[0][i];
  printf("% 10le\n",sqrt(aux));
*/  
  free(x); /*free(check);*/ free(b); free(err[0]);

  return 0;
 
}

/* matrix-vector product subroutine */
void A_product(double* v, int dim, double* Av)
{
  int     i;
  int     j;
  
  /* make the product */
  for (i=0; i<dim; i++)
  {
    Av[i] = 0;
    for (j=0; j<dim; j++) Av[i]+=A[i][j]*v[j];
  }
}
