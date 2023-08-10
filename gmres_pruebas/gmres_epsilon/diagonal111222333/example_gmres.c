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

/* X of XDX^-1=A */
/*double* X;*/
/* D of XDX^-1=A */
double* D;
/* epsilon test */  
double  epsilon;

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
  /* get the vector D */
  D = (double*)malloc(dim*sizeof(double));
  if (D == NULL) { puts("memory problem"); exit(1); }
  for (i=0; i<dim; i++) fscanf(f,"%lf",D+i);
  /* close the file */
  fclose(f);

  /* make the vector b */
  b = (double*)malloc(dim*sizeof(double));
  if (b == NULL) { puts("memory problem"); exit(1); }

  for (i=0; i<dim-1; i++) b[i]=D[i]+D[i+1];
  b[dim-1]=D[dim-1]+D[0];

/*
  b[0]=1;
  b[1]=2;
  b[2]=3;
  b[3]=4;
  b[4]=5;
  b[5]=6;
  b[6]=7;
  b[7]=8;
  b[8]=9;
*/

  /* GMRESm 'm' epsilon text */
  epsilon=0;
/* 
  int j;
  double* aux;
  aux = (double*)malloc(dim*sizeof(double));
  if (aux == NULL) { puts("memory problem"); exit(1); }
  aux[0]=1;
  for(i=1;i<dim;i++) aux[i]=0;
  A_product(aux,dim,b);
  for(j=0; j<dim; j++) printf("% 10le ",b[j]);
  puts("");
  for(i=1;i<dim;i++)
  {
    aux[i-1]=0;
    aux[i]=1;
    A_product(aux,dim,b);
    for(j=0; j<dim; j++) printf("% 10le ",b[j]);
    puts("");
  }

  return 0;
*/
  printf("%le ",epsilon);
  for (i=1; i<=dim; i++)
  {
    printf("%d ",i);
    gmresM(i,dim,b,A_product,0.25,1e-14,NULL);
  }
  puts("");
  for (epsilon=1e-8; epsilon < 1e-1+1e-2; epsilon*=10)
  {
    printf("%le ",epsilon);
    for (i=1; i<=dim; i++)
    {
      printf("%d ",i);
      gmresM(i,dim,b,A_product,0.25,1e-14,NULL);
    }
    puts("");
  }

  free(D); free(b);

  return 0;
 
}

/* matrix-vector product subroutine */
void A_product(double* v, int dim, double* Av)
{
  int     i;
  
  /* make the product */
  /* D product */
  for (i=0; i<dim; i+=3) Av[i] = v[i]*(D[i]-epsilon);
  for (i=1; i<dim; i+=3) Av[i] = v[i]*D[i];
  for (i=2; i<dim; i+=3) Av[i] = v[i]*(D[i]+epsilon);
}
