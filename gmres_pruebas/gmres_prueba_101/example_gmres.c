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
  int     i;      /* auxiliary variable */
  int     j;      /* auxiliary variable */
  FILE*   f;      /* pointer to file */
  double* b;      /* starting vector */
  int     dim;    /* dimension of vector b */
  double* x;      /* solution */
 
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
/*gmresUnlimited(dim,b,A_product,0.25,1e-13);*/

  /* GMRESm 'm' error text */
  for (i=2; i<=dim; i++)
  {
/*    printf("errors Restarted GMRES(%d) ...\n",i);
*/    x = gmresM(i,dim,b,A_product,0.25,1e-13,NULL);
  }

  /* show the solution */
/*  if (x == NULL) puts("no solution found for the system Ax=b");
  else
  {
    puts("solution found of Ax=b; x= "); 
    for (i=0; i<dim; i++) printf("% 10le ", x[i]);
    puts("");
  }
*/ 
  /* free memory */
  for (i=0; i<dim; i++) free(A[i]);
  free(A); free(x); free(b);

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
