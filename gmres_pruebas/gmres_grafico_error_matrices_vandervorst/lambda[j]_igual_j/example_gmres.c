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

/* beta of A=SAS^(-1) definition */
double beta;
/* vector of lambdas of A=SAS^(-1) definition */
double* l;

int main(void)
{
  int     i;      /* auxiliary variable */
  FILE*   f;      /* pointer to file */
  double* b;      /* starting vector */
  int     dim;    /* dimension of vector b */
 
  /* open the file with the parameters in this order: dim, beta, l[0], ..., l[dim-1] */
  f = fopen("parameters.txt","r");
  if (f == NULL) { puts("problem to open the file"); exit(1); }
  /* get the dimension of the system */
  fscanf(f,"%d",&dim); 
  /* get beta parameter */
  fscanf(f,"%lf",&beta); 
  /* get the vector lambda */
  l = (double*)malloc(dim*sizeof(double));
  if (l == NULL) { puts("memory problem"); exit(1); }
  for (i=0; i<dim; i++) fscanf(f,"%lf",l+i);
  /* close the file */
  fclose(f);

  /* make the vector b */
  b = (double*)malloc(dim*sizeof(double));
  if (b == NULL) { puts("memory problem"); exit(1); }
  for (i=0; i<dim-1; i++) b[i]=l[i]+beta*l[i+1];
  b[dim]=l[dim];
 
  /* GMRES error test */
  gmresUnlimited(dim,b,A_product,0.25,1e-14);
 
  /* GMRESm 'm' error text */
/*
  for (i=1; i<=dim; i++)
  {
    printf("errors Restarted GMRES(%d) ...\n",i);
    x = gmresM(i,dim,b,A_product,0.25,1e-14,NULL);
  }
*/

  free(l); free(b);

  return 0;
 
}

/* matrix-vector product subroutine */
void A_product(double* v, int dim, double* Av)
{
  int     i;
  int     j;
  double  aux;
  
  /* make the product */
  /* S^-1 product */
  for (i=0; i<dim; i++)
  {
    Av[i] = 0;
    aux   = 1;
    for(j=i; j<dim; j++) 
    {
      Av[i] += aux*v[j];
      aux   *= -beta;
    }
  }
  /* B product */
  for (i=0; i<dim; i++) Av[i] *= l[i];
  /* S product */
  for (i=0; i<dim-1; i++) Av[i] += beta*Av[i+1];
}
