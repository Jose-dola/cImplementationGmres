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

/* beta of X of XDX^-1=A */
double beta;
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
  beta=0.9;
  /* GMRESm 'm' epsilon text */
/*  for(beta=0; beta<=1; beta+=0.01)
  {
    epsilon=0;
    printf("#beta=%lf\n",beta);
    printf("%le %le ",beta,epsilon);
    for (i=1; i<=dim; i++)
    {
      printf("%d ",i);
      gmresM(i,dim,b,A_product,0.25,1e-12,NULL);
    }
    puts("");
*/
      epsilon=0;    
      for (i=1; i<=dim; i++)
      {
        printf("%le %d ",epsilon,i);
        gmresM(i,dim,b,A_product,0.25,1e-12,NULL);
        puts("");
      }
      puts("");

    for (epsilon=1e-9; epsilon < 1+1e-2; epsilon*=10)
    {
      for (i=1; i<=dim; i++)
      {
        printf("%le %d ",epsilon,i);
        gmresM(i,dim,b,A_product,0.25,1e-12,NULL);
        puts("");
      }
      puts("");
    }
/*    puts("");
  }
*/
  free(D); free(b);

  return 0;
}

/* matrix-vector product subroutine */
void A_product(double* v, int dim, double* Av)
{
  int     i;
  int     j;
  double  aux;
 
  /* make the product */
  /* X^-1 product */
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
  /* D product */
  for (i=0; i<dim; i+=3) Av[i] = Av[i]*(D[i]-epsilon);
  for (i=1; i<dim; i+=3) Av[i] = Av[i]*D[i];
  for (i=2; i<dim; i+=3) Av[i] = Av[i]*(D[i]+epsilon);
  /* X product */
  for (i=0; i<dim-1; i++) Av[i] += beta*Av[i+1];
}
