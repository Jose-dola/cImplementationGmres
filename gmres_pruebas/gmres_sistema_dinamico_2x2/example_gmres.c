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
 
  /* allocate vector b */
  b = (double*)malloc(dim*sizeof(double));
  if (b == NULL) { puts("memory problem"); exit(1); }
  /* open the file with the matrix A */
  f = fopen("A.txt","r");
  if (f == NULL) { puts("problem to open 'A.txt' file"); exit(1); }
  /* get the dimension */
  fscanf(f,"%d",&dim); 
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

  /* GMRESm dynamic system test */
  for (b[0]=-1; b[0]<=1; b[0]+=0.01)
    for (b[1]=-1; b[1]<=1; b[1]+=0.01)
      printf("%.2lf  %.2lf  %13le\n",b[0],b[1],gmresM(1,dim,b,A_product,0.25,1e-13,NULL));

  /* free memory */
  for (i=0; i<dim; i++) free(A[i]);
  free(A); free(b);

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
