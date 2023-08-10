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

void showAperturbed(int dim);

/* A_product header */
void A_product(double* v, int dim, double* Av);

/* matrix */
double** A;
double*  w;

/* epsilon test */  
double  epsilon;

int main(void)
{
  int     i;      /* auxiliary variable */
  int     j;
  FILE*   f;      /* pointer to file */
  double* b;      /* starting vector */
  int     dim;    /* dimension of vector b */

  f = fopen("parameters.txt","r");
  if (f == NULL) { puts("problem to open the file"); exit(1); }
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
  /* make vector w to perturbate A */
  w = (double*)malloc(dim*sizeof(double));
  if (w == NULL) { puts("memory problem"); exit(1); }
  for (i=0; i<dim; i++) fscanf(f,"%lf",&w[i]);
  /* close the file */
  fclose(f);

  /* make the vector b */
  b = (double*)malloc(dim*sizeof(double));
  if (b == NULL) { puts("memory problem"); exit(1); }
  for (i=0; i<dim; i++) b[i]=i*0.1;

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
*/  /* show w */
/*  puts("w:");
  for (i=0; i<dim; i++) printf("% 10lf ", w[i]);
  puts("");
*/
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
    /* GMRESm 'm' epsilon text */
    for (epsilon=-1e-4; epsilon < -1e-6+1e-12; epsilon/=sqrt(10))
    {
      A[4][5]=2*(1+epsilon); A[5][4]=-2*(1+epsilon);
      A[6][7]=2*(1-epsilon); A[7][6]=-2*(1-epsilon);
/*      showAperturbed(dim);*/
      for (i=2; i<=dim; i++)
      {
        printf("%d ",i);
        gmresM(i,dim,b,A_product,0.25,1e-13,NULL);
        puts("");
      }
      puts("");
    }
    epsilon=0;
    A[4][5]=2*(1+epsilon); A[5][4]=-2*(1+epsilon);
    A[6][7]=2*(1-epsilon); A[7][6]=-2*(1-epsilon);
/*    showAperturbed(dim);*/
    for (i=2; i<=dim; i++)
    {
      printf("%d ",i);
      gmresM(i,dim,b,A_product,0.25,1e-13,NULL);
      puts("");
    }
    puts("");
    for (epsilon=1e-6; epsilon < 1e-4+1e-12; epsilon*=sqrt(10))
    {
      A[4][5]=2*(1+epsilon); A[5][4]=-2*(1+epsilon);
      A[6][7]=2*(1-epsilon); A[7][6]=-2*(1-epsilon);
/*      showAperturbed(dim);*/
      for (i=2; i<=dim; i++)
      {
        printf("%d ",i);
        gmresM(i,dim,b,A_product,0.25,1e-13,NULL);
        puts("");
      }
      puts("");
    }


  free(A); free(b);

  return 0;
}


void showAperturbed(int dim)
{
  int i;
  int j;

  puts("A:");
  for (i=0; i<dim; i++)
  {
    for(j=0; j<1; j++) printf("% 10lf ", A[i][j]);
    printf("% 10lf ", A[i][1]+epsilon*w[i]);
    for (j=2; j<dim; j++) printf("% 10lf ",A[i][j]);
    puts("");
  }
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
    for (j=0; j<1; j++) Av[i]+=A[i][j]*v[j];
    Av[i]+=(A[i][1]+epsilon*w[i])*v[1];
    for (j=2; j<dim; j++) Av[i]+=A[i][j]*v[j];
  }
}
