/*
 * Functions to solve Ax=b with gmres method
 * 
 * by Jose Luis Dorado
 * University of Barcelona
 * February 2018
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define TOL 1e-15

/* 
 * GMRES algorithm to solve Ax=b
 *
 * PARAM:
 *   dim:       dimension of the system
 *   b:         b vector
 *   A_product: the user define the matrix A like a subrutine that computes the product of the matrix A and a vector v
 *              the parameters of this function are:
 *                first:  pointer to double that represents the vector v
 *                second: integer that is the dimension of the system (wich is the same as dim)
 *                third:  pointer to double where the function save the result (the space is already reserved 
 *                        before the fucntion called)
 *   rc:        refinement condition of arnoldi (for example 0.25)
 *   tol:       error tolerance
 *
 * RETURN:
 *   the solution of the system
*/
double* gmresUnlimited( int dim, double* b, void (*A_product)(double*,int,double*), double rc, double tol )
{
  int      i;
  int      j;
  int      k;
  double   aux;
  double*  x;      /* global solution */
  double*  bLS;    /* after QR givens (H=QR), Hy=||b||*e1 --> QRy=||b||*e1 --> Ry=bLS  */
  double*  y;      /* solution of Hy=||b||*e1 */
  double*  sine;   /* sine and cosine represent Q of QR givens */
  double*  cosine; /* sine and cosine represent Q of QR givens */
  double** R;      /* R of QR givens */
  double** Q;      /* arnoldi's Q */
  double** H;      /* arnoldi's H */

  /* allocate Q and H for the arnoldi iteration */
  /* Q has dim+1 columns (maximum) */
  Q = (double**)malloc((dim+1)*sizeof(double*));
  if (Q == NULL)      { puts("memory problem"); exit(1); }
  H = (double**)malloc((dim)*sizeof(double*));
  if (H == NULL)      { puts("memory problem"); exit(1); }

  /* allocate Q[0] */
  Q[0] = (double*)malloc((dim)*sizeof(double));
  if (Q[0] == NULL)   { puts("memory problem"); exit(1); }

  /* allocate R, sine and cosine for the QR givens */
  R = (double**)malloc((dim)*sizeof(double*));
  if (R == NULL)      { puts("memory problem"); exit(1); }
  sine = (double*)malloc((dim)*sizeof(double));
  if (sine == NULL)   { puts("memory problem"); exit(1); }
  cosine = (double*)malloc((dim)*sizeof(double));
  if (cosine == NULL) { puts("memory problem"); exit(1); }

  /* we use the vector bLS to solve least square problem  
   * we put the solution in the vector y 
  */
  bLS = (double*)malloc((dim)*sizeof(double));
  if (bLS == NULL)    { puts("memory problem"); exit(1); }
  y   = (double*)malloc((dim)*sizeof(double));
  if (y == NULL)      { puts("memory problem"); exit(1); }
 
  /* Q[0] = b/||b|| and bLS[0] = ||b|| */
  bLS[0] = 0;
  for (i=0; i<dim; i++) bLS[0]+=b[i]*b[i];
  bLS[0] = sqrt(bLS[0]); 
  for (i=0; i<dim; i++) Q[0][i]=b[i]/bLS[0];
  
  /* set iterator 'j' */
  j=0;
  /*************GMRES UNLIMITED ITERATIONS*******************/
  while (1)
  {
    /***************** ARNOLDI ITERATION ********************/
    /********************************************************/

    Q[j+1] = (double*)malloc((dim)*sizeof(double));
    if (Q[j+1] == NULL)   { puts("memory problem"); exit(1); }
    A_product(Q[j],dim,Q[j+1]);
    /* save ||Q[j+1]|| for the refinement check */
    aux = 0;
    for (k=0; k<dim; k++) aux+=Q[j+1][k]*Q[j+1][k];
    aux = sqrt(aux);

    /* allocate space for next column of H */
    H[j] = (double*)malloc((2+j)*sizeof(double));
    if (H[j] == NULL) { puts("memory problem"); exit(1); }

    /** orthogonalization **/
    for (k=0; k<j; k++)
    {
      /* product of Q[k]*Q[j+1] */
      H[j][k] = 0;
      for (i=0; i<dim; i++) H[j][k]+=Q[k][i]*Q[j+1][i];
      /* Q[j+1] = Q[j+1] - H[j][k]*Q[k] */
      for (i=0; i<dim; i++) Q[j+1][i]=Q[j+1][i]-H[j][k]*Q[k][i];
    }
    /* last iteration of the previous "for" with k=j (and code changes) to calculate H[j][j] and H[j][j+1]=||Q[j+1]|| */
    /* product of Q[j]*Q[j+1] */
    H[j][j] = 0;
    for (i=0; i<dim; i++) H[j][j]+=Q[j][i]*Q[j+1][i];
    /* Q[j+1] = Q[j+1] - H[j][i]*Q[i]  and  H[j][j+1] = ||Q[j+1]|| */
    H[j][j+1] = 0;
    for (i=0; i<dim; i++) { Q[j+1][i]=Q[j+1][i]-H[j][j]*Q[j][i]; H[j][j+1]+=Q[j+1][i]*Q[j+1][i]; }
    H[j][j+1] = sqrt(H[j][j+1]);

    /*** refinement ***/
    if ( H[j][j+1]/aux <= rc)
    {
      for (k=0; k<j; k++)
      {
        /* aux = Q[k]*Q[j+1] */
        aux = 0;
        for (i=0; i<dim; i++) aux+=Q[k][i]*Q[j+1][i];
        /* Q[j+1] = Q[j+1] - aux*Q[k] */
        for (i=0; i<dim; i++) Q[j+1][i]=Q[j+1][i]-aux*Q[k][i];
        H[j][k] = H[j][k] + aux;
      }
      /* last iteration of the previous "for" with k=j (and code changes) to calculate H[j][j] and H[j][j+1]=||Q[j+1]|| */
      /* aux = Q[j]*Q[j+1] */
      aux = 0;
      for (i=0; i<dim; i++) aux+=Q[j][i]*Q[j+1][i];
      /* Q[j+1] = Q[j+1] - aux*Q[j] and get ||Q[j+1]|| */
      H[j][j+1] = 0;
      for (i=0; i<dim; i++) { Q[j+1][i]=Q[j+1][i]-aux*Q[j][i]; H[j][j+1]+=Q[j+1][i]*Q[j+1][i]; }
      H[j][j] = H[j][j] + aux;
      /* H[j][j+1] = ||Q[j+1]|| */
      H[j][j+1] = sqrt(H[j][j+1]);
    }
    
    /************************** BREAKDOWN ***************************/
    /**** if H[j][j+1]=0 or we have a 'dim' iterations --> break ****/
    if ( fabs(H[j][j+1]) < tol || j == dim-1) { break; }

    /* normalize Q[j+1] */
    for (i=0; i<dim; i++) Q[j+1][i] = Q[j+1][i]/H[j][j+1];

    /** GIVENS ROTATIONS TO FIND R AND bLS TO SOLVE LEAST SQUARE PROBLEM  **/
    /***********************************************************************/
    
    /* allocate space for next column of R */
    R[j] = (double*)malloc((1+j)*sizeof(double));
    if (R[j] == NULL) { puts("memory problem"); exit(1); }

    /* compute the product of the previous Givens rotations with H[j] */
    R[j][0]=H[j][0];
    for (k=0; k<j; k++)
    {
      aux       = R[j][k];
      R[j][k]   = aux*cosine[k] + H[j][k+1]*sine[k];
      R[j][k+1] = H[j][k+1]*cosine[k] - aux*sine[k];
    }

    /* compute and apply Givens rotation to put 0 in the subdiagonal */
    if ( fabs(H[j][j+1]) > fabs(R[j][j]) ) { aux=R[j][j]/H[j][j+1]; sine[j]=1./sqrt(1+aux*aux); cosine[j]=sine[j]*aux;  }
    else                                   { aux=H[j][j+1]/R[j][j]; cosine[j]=1./sqrt(1+aux*aux); sine[j]=cosine[j]*aux; }
    /* apply the rotation found (remember that bLS[i+1]=0) */
    R[j][j]  = R[j][j]*cosine[j] + H[j][j+1]*sine[j];
    bLS[j+1] = -bLS[j]*sine[j];
    bLS[j]   = bLS[j]*cosine[j];
    
    /* increase iterator for the new iteration */
    j++;
  }
  /************ END GMRES UNLIMITED ITERATIONS *************/

  /* compute last R column */
  /*************************/
  /* allocate space */
  R[j] = (double*)malloc((1+j)*sizeof(double));
  if (R[j] == NULL) { puts("memory problem"); exit(1); }
  /* compute */
  R[j][0]=H[j][0];
  for (k=0; k<j; k++)
  {
    aux       = R[j][k];
    R[j][k]   = aux*cosine[k] + H[j][k+1]*sine[k];
    R[j][k+1] = H[j][k+1]*cosine[k] - aux*sine[k];
  }

  /*********** SOLVE LEAST SQUARE PROBLEM ****************/
  /*******************************************************/
  /* solve LS problem solving the system Ry=bLS */
  for (k=j; k>=0; k--)
  {
    y[k] = bLS[k];
    for (i=k+1; i<=j; i++) y[k]-=y[i]*R[i][k];
    y[k] = y[k]/R[k][k];
  }

  /* we use the space allocated in Q[j+1] for the solution */
  x = Q[j+1];

  /* store the solution in x // x=Qy */
  for (k=0; k<dim; k++)
  {
    x[k] = 0;
    for (i=0; i<=j; i++) x[k]+=Q[i][k]*y[i];
  }

  /*** FREE MEMORY ***/
  for(i=0; i<=j; i++) { free(R[i]); free(H[i]); free(Q[i]); }
  free(Q);
  free(R);
  free(H);
  free(sine);
  free(cosine);
  free(bLS);
  free(y);

  return x;
}


/* 
 * RESTARTED GMRES dynamic system text
 *
 * PARAM:
 *   m:         restart parameter
 *   dim:       dimension of the system
 *   b:         b vector
 *   A_product: the user define the matrix A like a subrutine that computes the product of the matrix A and a vector v
 *              the parameters of this function are:
 *                first:  pointer to double that represents the vector v
 *                second: integer that is the dimension of the system (wich is the same as dim)
 *                third:  pointer to double where the function save the result (the space is already reserved 
 *                        before the fucntion called)
 *   rc:        refinement condition of arnoldi (for example 0.25)
 *   tol:       error tolerance
 *   error:     we put in this pointer the error vector of the solution returned (it is not necessary allocate)
 *              we ignore if the input is NULL
 *
 * RETURN:
 *   error
*/
double gmresM( int m, int dim, double* b, void (*A_product)(double*,int,double*), double rc, double tol, double** error )
{
  int      i;
  int      j;
  int      k;
  double   aux;
  double   previusBnorm=-1;
  int      mLS;     /* least square system dimension */
  double*  err;     /* iterative error vector */
  double*  x;       /* global solution */
  double*  bLS;     /* vector of Least Square system Ry=bLS  */
  double*  y;       /* Least Square problem solution */
  double*  sine;    /* sine and cosine represent Q of QR givens */
  double*  cosine;  /* sine and cosine represent Q of QR givens */
  double** R;       /* R of QR givens to solve the Least Square system */
  double** Q;       /* arnoldi's Q */
  double** H;       /* arnoldi's H */
 
  /* set least square system dimension as m-1 */
  mLS     = m-1;
 
  /* allocate error vector and copy b vector inside */
  err = (double*)malloc(dim*sizeof(double));
  if (err == NULL)      { puts("memory problem"); exit(1); }
  for(i=0; i<dim; i++)  { err[i]=b[i]; }

  /* allocate solution vector and set it equal to 0 */
  x = (double*)malloc(dim*sizeof(double));
  if (x == NULL)        { puts("memory problem"); exit(1); }
  for(i=0; i<dim; i++)  { x[i]=0; }

  /* allocate Q and H for the arnoldi iteration */
  Q = (double**)malloc((m+1)*sizeof(double*));
  if (Q == NULL)        { puts("memory problem"); exit(1); }
  H = (double**)malloc((m)*sizeof(double*));
  if (H == NULL)        { puts("memory problem"); exit(1); }
  /* allocate columns of Q */
  for (i=0; i<=m; i++)
  {
    Q[i] = (double*)malloc((dim)*sizeof(double));
    if (Q[i] == NULL)   { puts("memory problem"); exit(1); }
  }
  /* allocate columns of H */
  for (i=0; i<m; i++)
  {
    H[i] = (double*)malloc((i+2)*sizeof(double));
    if (H[i] == NULL)   { puts("memory problem"); exit(1); }
  }

  /* allocate R, sine and cosine for the QR givens and least square problem */
  R = (double**)malloc((m)*sizeof(double*));
  if (R == NULL)        { puts("memory problem"); exit(1); }
  for (i=0; i<m; i++)
  {
    R[i] = (double*)malloc((i+1)*sizeof(double));
    if (R[i] == NULL)   { puts("memory problem"); exit(1); }
  }
  sine = (double*)malloc((m)*sizeof(double));
  if (sine == NULL)     { puts("memory problem"); exit(1); }
  cosine = (double*)malloc((m)*sizeof(double));
  if (cosine == NULL)   { puts("memory problem"); exit(1); }

  /* we use the vector bLS to solve least square problem */
  /* we put the solution in the vector y */
  bLS = (double*)malloc((m+1)*sizeof(double));
  if (bLS == NULL)      { puts("memory problem"); exit(1); }
  y   = (double*)malloc((m)*sizeof(double));
  if (y == NULL)        { puts("memory problem"); exit(1); }

  /*** PREPARING THE FIRST GMRES ITERATION ***/
  /* bLS[0] = ||b|| */
  bLS[0] = 0;
  for (i=0; i<dim; i++) bLS[0]+=b[i]*b[i];
  bLS[0] = sqrt(bLS[0]);

  /****** RESTARTED GMRES LOOP ******/
  /* BREAKDOWN if overall error is less than the tolerance */
  while( bLS[0] > tol ) 
  {

    /* Q[0] = err/||err|| */
    for (i=0; i<dim; i++) Q[0][i]=err[i]/bLS[0];

    /******************* GMRES m ITERATIONS *******************/
    for(j=0; j<m; j++)
    {
      /***************** ARNOLDI ITERATION ********************/
      /********************************************************/
      A_product(Q[j],dim,Q[j+1]);
      /* save ||Q[j+1]|| for the refinement check */
      aux = 0;
      for (k=0; k<dim; k++) aux+=Q[j+1][k]*Q[j+1][k];
      aux = sqrt(aux);

      /** orthogonalization **/
      for (k=0; k<j; k++)
      {
        /* product of Q[k]*Q[j+1] */
        H[j][k] = 0;
        for (i=0; i<dim; i++) H[j][k]+=Q[k][i]*Q[j+1][i];
        /* Q[j+1] = Q[j+1] - H[j][k]*Q[k] */
        for (i=0; i<dim; i++) Q[j+1][i]=Q[j+1][i]-H[j][k]*Q[k][i];
      }
      /* last iteration of the previous "for" with k=j (and code changes) to calculate H[j][j] and H[j][j+1]=||Q[j+1]|| */
      /* product of Q[j]*Q[j+1] */
      H[j][j] = 0;
      for (i=0; i<dim; i++) H[j][j]+=Q[j][i]*Q[j+1][i];
      /* Q[j+1] = Q[j+1] - H[j][i]*Q[i]  and  H[j][j+1] = ||Q[j+1]|| */
      H[j][j+1] = 0;
      for (i=0; i<dim; i++) { Q[j+1][i]=Q[j+1][i]-H[j][j]*Q[j][i]; H[j][j+1]+=Q[j+1][i]*Q[j+1][i]; }
      H[j][j+1] = sqrt(H[j][j+1]);

      /*** refinement ***/
      if ( H[j][j+1]/aux <= rc)
      {
        for (k=0; k<j; k++)
        {
          /* aux = Q[k]*Q[j+1] */
          aux = 0;
          for (i=0; i<dim; i++) aux+=Q[k][i]*Q[j+1][i];
          /* Q[j+1] = Q[j+1] - aux*Q[k] */
          for (i=0; i<dim; i++) Q[j+1][i]=Q[j+1][i]-aux*Q[k][i];
          H[j][k] = H[j][k] + aux;
        }
        /* last iteration of the previous "for" with k=j (and code changes) to calculate H[j][j] and H[j][j+1]=||Q[j+1]|| */
        /* aux = Q[j]*Q[j+1] */
        aux = 0;
        for (i=0; i<dim; i++) aux+=Q[j][i]*Q[j+1][i];
        /* Q[j+1] = Q[j+1] - aux*Q[j] and get ||Q[j+1]|| */
        H[j][j+1] = 0;
        for (i=0; i<dim; i++) { Q[j+1][i]=Q[j+1][i]-aux*Q[j][i]; H[j][j+1]+=Q[j+1][i]*Q[j+1][i]; }
        H[j][j] = H[j][j] + aux;
        /* H[j][j+1] = ||Q[j+1]|| */
        H[j][j+1] = sqrt(H[j][j+1]);
      }

      /***** if H[j][j+1] = 0 --> go to solve Least Square problem *****/
      if ( fabs(H[j][j+1]) < tol ) 
      {
        /* compute last R column */
        /*************************/
        /* compute */
        R[j][0]=H[j][0];
        for (k=0; k<j; k++)
        {
          aux       = R[j][k];
          R[j][k]   = aux*cosine[k] + H[j][k+1]*sine[k];
          R[j][k+1] = H[j][k+1]*cosine[k] - aux*sine[k];
        }
        mLS     = j;
        break; 
      }

      /* normalize Q[j+1] */
      for (i=0; i<dim; i++) Q[j+1][i] = Q[j+1][i]/H[j][j+1];

      /** GIVENS ROTATIONS TO FIND R AND bLS TO SOLVE LEAST SQUARE PROBLEM  **/
      /***********************************************************************/

      /* compute the product of the previous Givens rotations with H[j] */
      R[j][0]=H[j][0];
      for (k=0; k<j; k++)
      {
        aux       = R[j][k];
        R[j][k]   = aux*cosine[k] + H[j][k+1]*sine[k];
        R[j][k+1] = H[j][k+1]*cosine[k] - aux*sine[k];
      }
 
      /* compute and apply Givens rotation to put 0 in the subdiagonal */
      if ( fabs(H[j][j+1]) > fabs(R[j][j]) ) { aux=R[j][j]/H[j][j+1]; sine[j]=1./sqrt(1+aux*aux); cosine[j]=sine[j]*aux;  }
      else                                   { aux=H[j][j+1]/R[j][j]; cosine[j]=1./sqrt(1+aux*aux); sine[j]=cosine[j]*aux; }
      /* apply the rotation found (remember that bLS[i+1]=0) */
      R[j][j]  = R[j][j]*cosine[j] + H[j][j+1]*sine[j];
      bLS[j+1] = -bLS[j]*sine[j];
      bLS[j]   = bLS[j]*cosine[j];

    }/***************** END m ITERATIONS OF GMRES ****************/

    /*********** SOLVE LEAST SQUARE PROBLEM ****************/
    /*******************************************************/
    /* solve LS problem solving the system Ry=bLS */
    for (k=mLS; k>=0; k--)
    {
      y[k] = bLS[k];
      for (i=k+1; i<=mLS; i++) y[k]-=y[i]*R[i][k];
      y[k] = y[k]/R[k][k];
    } 

    /* compute next aproximation of the solution */
    for (k=0; k<dim; k++)
      for (i=0; i<=mLS; i++) x[k]+=Q[i][k]*y[i];

    /* compute the error and prepare the next iteration */
    A_product(x,dim,err);
    for(i=0; i<dim; i++) err[i]=b[i]-err[i];

    /****** PREPARING THE NEXT ITERATION *******/
    /* Q[0] = err/||err|| and bLS[0] = ||err|| */
    bLS[0] = 0;
    for (i=0; i<dim; i++) bLS[0]+=err[i]*err[i];
    bLS[0] = sqrt(bLS[0]);

    for(i=0;i<dim;i++) printf("%le ",err[i]);
    printf(" ");    

    if(fabs(previusBnorm-bLS[0]) < tol) return bLS[0];
    previusBnorm = bLS[0];

  }/** END RESTARTED GMRES LOOP **/

  /*** save error vector or free it ***/
  if(error) *error=err;
  else      free(err);

  /*** FREE MEMORY ***/
  for(i=0; i<m; i++) { free(R[i]); free(H[i]); free(Q[i]); }
  free(Q[m]);
  free(Q);
  free(R);
  free(H);
  free(sine);
  free(cosine);
  free(bLS);
  free(y);

  return 0;
}

