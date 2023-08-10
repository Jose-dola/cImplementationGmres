/*
 * Functions to do Arnoldi iteration
 *
 * by Jose Luis Dorado
 * University of Barcelona
 * February 2018
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * Arnoldi’s method with modified Gram–Schmidt orthogonalization.
 *
 * PARAM:
 *   b:         starting vector
 *   dim:       dimension of vector b
 *   n:         number of iterations
 *   H:         we put here the submatrix of the associated Hessemberg matrix
 *              OBSERVATION: the first index represent the column and the second the row 
 *   Q:         we put here the Krylov orthonormal subspace which is resulted by the Arnoldi iteration 
 *              OBSERVATION: the first index represent the column and the second the row 
 *              the columns are the vectors of the Krylov orthonormal subspace 
 *   A_product: the user defines the matrix A like a function that represents the product of A with a vector (Av). 
 *              the parameters of this function are:
                  first:  pointer to the vector (double precision)
                  second: dimension of the vector (integer)
*/
void arnoldi( double* b, double dim, int n, double** H, double** Q, double* (*A_product)(double*,int) )
{
  int      i;
  int      j;
  int      k;
  double   aux;

  /* Q[0] = ||b|| */
  Q[0] = (double*)malloc((dim)*sizeof(double));
  if (Q[0] == NULL) { puts("memory problem"); exit(1); }
  aux = 0;
  for (i=0; i<dim; i++) aux+=b[i]*b[i];
  aux = sqrt(aux);
  for (i=0; i<dim; i++) Q[0][i]=b[i]/aux;

  /* ARNOLDI ITERATION */
  for (i=0; i<n; i++)
  {
    Q[i+1] = A_product(Q[i],dim);
    H[i]   = (double*)malloc((2+i)*sizeof(double));
    if (H[i] == NULL) { puts("memory problem"); exit(1); }
    
    for (j=0; j<i; j++)
    {
      /* product of Q[j]*Q[i+1] */
      H[i][j] = 0;
      for (k=0; k<dim; k++) H[i][j]+=Q[j][k]*Q[i+1][k];
      /* Q[i+1] = Q[i+1] - H[i][j]*Q[j] */
      for (k=0; k<dim; k++) Q[i+1][k]=Q[i+1][k]-H[i][j]*Q[j][k]; 
    }
    /* last iteration of the previous "for" with j=i */
    /* product of Q[j]*Q[i+1] */
    H[i][i] = 0;
    for (k=0; k<dim; k++) H[i][i]+=Q[i][k]*Q[i+1][k];
    /* Q[i+1] = Q[i+1] - H[i][j]*Q[j] and get ||Q[i+1]||^2 */
    aux = 0;
    for (k=0; k<dim; k++) { Q[i+1][k]=Q[i+1][k]-H[i][i]*Q[i][k]; aux+=Q[i+1][k]*Q[i+1][k]; }

    /* H[i][i+1] = ||Q[i+1]|| */
    H[i][i+1] = sqrt(aux);
    /* normalize Q[i+1] */
    for (k=0; k<dim; k++) Q[i+1][k] = Q[i+1][k]/H[i][i+1];
  }
}

/*
 * The Arnoldi Method with refined modified Gram–Schmidt
 *
 * PARAM:
 *   b:         starting vector
 *   dim:       dimension of vector b
 *   n:         number of iterations
 *   H:         we put here the submatrix of the associated Hessemberg matrix
 *              OBSERVATION: the first index represent the column and the second the row 
 *   Q:         we put here the Krylov orthonormal subspace which is resulted by the Arnoldi iteration 
 *              OBSERVATION: the first index represent the column and the second the row 
 *              the columns are the vectors of the Krylov orthonormal subspace 
 *   rc:        refinement condition. For example 0.25
 *   A_product: the user defines the matrix A like a function that represents the product of A with a vector (Av). 
 *              the parameters of this function are:
                  first:  pointer to the vector (double precision)
                  second: dimension of the vector (integer)
*/
void arnoldiRefined( double* b, double dim, int n, double** H, double** Q, double rc, double* (*A_product)(double*,int) )
{ 
  int    i;
  int    j;
  int    k;
  double aux;
  
  /* Q[0] = ||b|| */
  Q[0] = (double*)malloc((dim)*sizeof(double));
  if (Q[0] == NULL) { puts("memory problem"); exit(1); }
  aux = 0;
  for (i=0; i<dim; i++) aux+=b[i]*b[i];
  aux = sqrt(aux); 
  for (i=0; i<dim; i++) Q[0][i]=b[i]/aux;
  
  /*** ARNOLDI ITERATION ***/
  for (i=0; i<n; i++)
  {
    Q[i+1] = A_product(Q[i],dim);
    /* save ||Q[i+1]|| for the refinement check */
    aux = 0;
    for (k=0; k<dim; k++) aux+=Q[i+1][k]*Q[i+1][k];
    aux = sqrt(aux);
 
    /* allocate space for next column of H */
    H[i]   = (double*)malloc((2+i)*sizeof(double));
    if (H[i] == NULL) { puts("memory problem"); exit(1); }

    /** orthogonalization **/
    for (j=0; j<i; j++)
    {
      /* product of Q[j]*Q[i+1] */
      H[i][j] = 0;
      for (k=0; k<dim; k++) H[i][j]+=Q[j][k]*Q[i+1][k];
      /* Q[i+1] = Q[i+1] - H[i][j]*Q[j] */
      for (k=0; k<dim; k++) Q[i+1][k]=Q[i+1][k]-H[i][j]*Q[j][k];
    }
    /* last iteration of the previous "for" with j=i (and code changes) to calculate H[i][i] and H[i][i+1]=||Q[i+1]|| */
    /* product of Q[j]*Q[i+1] */
    H[i][i] = 0;
    for (k=0; k<dim; k++) H[i][i]+=Q[i][k]*Q[i+1][k];
    /* Q[i+1] = Q[i+1] - H[i][j]*Q[j]  and  H[i][i+1] = ||Q[i+1]|| */
    H[i][i+1] = 0;
    for (k=0; k<dim; k++) { Q[i+1][k]=Q[i+1][k]-H[i][i]*Q[i][k]; H[i][i+1]+=Q[i+1][k]*Q[i+1][k]; }
    H[i][i+1] = sqrt(H[i][i+1]);

    /*** refinement ***/
    if ( H[i][i+1]/aux <= rc)
    {
      for (j=0; j<i; j++)
      {
        /* aux = Q[j]*Q[i+1] */
        aux = 0;
        for (k=0; k<dim; k++) aux+=Q[j][k]*Q[i+1][k];
        /* Q[i+1] = Q[i+1] - aux*Q[j] */
        for (k=0; k<dim; k++) Q[i+1][k]=Q[i+1][k]-aux*Q[j][k];
        H[i][j] = H[i][j] + aux;
      }     
      /* last iteration of the previous "for" with j=i (and code changes) to calculate H[i][i] and H[i][i+1]=||Q[i+1]|| */
      /* aux = Q[j]*Q[i+1] */
      aux = 0;
      for (k=0; k<dim; k++) aux+=Q[i][k]*Q[i+1][k];
      /* Q[i+1] = Q[i+1] - aux*Q[j] and get ||Q[i+1]|| */
      H[i][i+1] = 0;
      for (k=0; k<dim; k++) { Q[i+1][k]=Q[i+1][k]-aux*Q[i][k]; H[i][i+1]+=Q[i+1][k]*Q[i+1][k]; }
      H[i][i] = H[i][i] + aux;
      /* H[i][i+1] = ||Q[i+1]|| */
      H[i][i+1] = sqrt(H[i][i+1]);
    }
    
    /* normalize Q[i+1] */
    for (k=0; k<dim; k++) Q[i+1][k] = Q[i+1][k]/H[i][i+1];
  }
}

/* function to allocate Q matrix. This function should be called before arnoldi iteration 
 * PARAM: number of iterations
*/
double** Qalloc(int n)
{
  double** Q;
  Q = (double**)malloc((n+1)*sizeof(double));
  if (Q == NULL) { puts("memory problem"); exit(1); }
  return Q;
}

/* function to allocate H matrix. This function should be called before arnoldi iteration 
 * PARAM: number of iterations
*/
double** Halloc(int n)
{
  double** H;
  H = (double**)malloc(n*sizeof(double));
  if (H == NULL) { puts("memory problem"); exit(1); }
  return H;
}

/* function to free memory
 * PARAM:
 *   n:  number of iterations
 *   Q:  pointer to Q matrix
 *   H:  pointer to Hn matrix
*/
void free_arnoldi(int n, double** H, double** Q)
{
  int i;
  for (i=0; i<n; i++) { free(Q[i]); free(H[i]); } 
  free(Q[n]); 
  free(Q); free(H);
}

