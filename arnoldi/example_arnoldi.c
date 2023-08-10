/*
 * Example of Arnoldi iteration
 *
 * by Jose Luis Dorado
 * University of Barcelona
 * February 2018
*/

#include <stdio.h>
#include <stdlib.h>
#include "arnoldi.h"

/* A_product header */
double* A_product(double* v, int dim);
/* global variables */
double** A; /* matrix A */

int main(void)
{
  int       i;   /* auxiliary variable */
  int       j;   /* auxiliary variable */
  FILE*     f;   /* pointer to file */
  double*   b;   /* starting vector */
  int       dim; /* dimension of vector b */
  int       n;   /* number of iterations */
  double**  Hn;  /* submatrix of the associated Hessemberg matrix
                  * OBSERVATION: the first index represent the column and the second the row */  
  double**  Q;   /* the Krylov orthonormal subspace which is resulted by the Arnoldi iteration 
                  * OBSERVATION: the first index represent the column and the second the row 
                  * the columns are the vectors of the Krylov orthonormal subspace */
  double**  R;   /* R of QR=H computed by Givens rotation */
  double*   s;   /* It is an array with the sine of the Givens rotations */
  double*   c;   /* It is an array with the cosine of the Givens rotations */
 
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

  /* define as many iterations as the dimension */
  n=dim;

  /* show A */
  puts("A:");
  for (i=0; i<dim; i++)
  {
    for (j=0; j<dim; j++) printf("% 10lf ",A[i][j]);
    puts(""); 
  }
  puts("");

  /************ do the arnoldi iteration *************/
  Q  = Qalloc(n);
  Hn = Halloc(n);
  arnoldi(b, dim, n, Hn, Q, A_product);
 
  puts("");
  puts("ARNOLDI ITERATION:");

  /* show Q */
  puts("Q:");
  for (i=0; i<dim; i++)
  {
    for (j=0; j<dim+1; j++) printf("% 10lf ",Q[j][i]);
    puts(""); 
  }
  puts("");
  /* show H */
  puts("H:");
  for (i=0; i<dim+1; i++)
  { 
    j = 0;
    if  (0<=i-2) { j++; printf(" 000000000 "); }
    for (; j<=i-2; j++) printf(" 000000000 ");
    for (; j<dim; j++) printf("% 10lf ",Hn[j][i]);
    puts(""); 
  }

  /* free memory */
  free_arnoldi(n,Hn,Q);
  /****************************************************/

  /************ do the arnoldi iteration with refined modified Gram–Schmidt *************/ 
  Q  = Qalloc(n);
  Hn = Halloc(n);
  arnoldiRefined(b, dim, n, Hn, Q, 0.25, A_product);

  puts("");
  puts("ARNOLDI ITERATION with refined modified Gram–Schmidt:");
  
  /* show Q */
  puts("Q:");
  for (i=0; i<dim; i++)
  {
    for (j=0; j<dim+1; j++) printf("% 10lf ",Q[j][i]);
    puts("");
  }
  puts("");
  /* show H */
  puts("H:");
  for (i=0; i<dim+1; i++)
  {
    j = 0;
    if  (0<=i-2) { j++; printf(" 000000000 "); }
    for (; j<=i-2; j++) printf(" 000000000 ");
    for (; j<dim; j++) printf("% 10lf ",Hn[j][i]);
    puts("");
  }

  /* free memory */
  free_arnoldi(n,Hn,Q);
  /***************************************************************************************/

  /************ do the arnoldi iteration with refined modified Gram–Schmidt and test givens rotations *********/
  R = (double**)malloc(n*sizeof(double*));
  if (R == NULL) { puts("memory problem"); exit(1); }
  s = (double*)malloc(n*sizeof(double));
  if (s == NULL) { puts("memory problem"); exit(1); }
  c = (double*)malloc(n*sizeof(double));
  if (c == NULL) { puts("memory problem"); exit(1); }
  Q  = Qalloc(n);
  Hn = Halloc(n);
  givensTest(b, dim, n, Hn, Q, R, s, c, 0.25, A_product);

  puts("");
  puts("ARNOLDI ITERATION with refined modified Gram–Schmidt and test givens rotations:");

  /* show Q */
  puts("Q:");
  for (i=0; i<dim; i++)
  {
    for (j=0; j<dim+1; j++) printf("% 10lf ",Q[j][i]);
    puts("");
  }
  puts("");
  /* show H */
  puts("H:");
  for (i=0; i<dim+1; i++)
  {
    j = 0;
    if  (0<=i-2) { j++; printf(" 000000000 "); }
    for (; j<=i-2; j++) printf(" 000000000 ");
    for (; j<dim; j++)  printf("% 10lf ",Hn[j][i]);
    puts("");
  }
  /* show R */
  puts("R:");
  for (i=0; i<dim; i++)
  {
    j = 0;
    if  (0<=i-1) { j++; printf(" 000000000 "); }
    for (; j<=i-1; j++) printf(" 000000000 ");
    for (; j<dim; j++)  printf("% 10lf ",R[j][i]);
    puts("");
  }
  /* show QR */
  puts("QR: (is equal than H?)");
  double a;
  double d;
  double aux;
  for (i=n-1; i>=0; i--)
  {
    aux = 0;
    for (j=i; j>=0; j--)
    {
      a   = R[i][j]*c[j] - aux*s[j];
      d   = R[i][j]*s[j] + aux*c[j];
      aux = a;
      printf("% 10lf ",d);
    }
    printf("% 10lf \n",a);
  }

  /* free memory */
  free_arnoldi(n,Hn,Q);
  for (i=0; i<n; i++) free(R[i]); 
  free(R); free(s); free(c);
  /***************************************************************************************/


  /* free matrix A */
  for (i=0; i<dim; i++) free(A[i]);
  free(A);
  
  return 0;
}

/* matrix-vector product subroutine */
double* A_product(double* v, int dim)
{
  int     i;
  int     j;
  double* res;
  
  /* allocate result vector */
  res = (double*)malloc(dim*sizeof(double));
  if (res == NULL) { puts("memory problem"); exit(1); }
  /* make the product */
  for (i=0; i<dim; i++)
  {
    res[i] = 0;
    for (j=0; j<dim; j++) res[i]+=A[i][j]*v[j];
  }
  
  return res; 
}
