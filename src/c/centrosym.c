#include "symtrx.h"

#define PRECISION 1e-12
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) > (b) ? (b) : (a))


/*!
 * Compute the actual size of a centrosym square matrix.
 *
 * \param[in]  dim The dimensions.
 * \retval Its size in memory.
 */
long centrosym_size(int dim)
{

  return dim * (dim+1) / 2;

}


/*!
 * Allocate space for centrosymmetric square matrix.
 *
 * \param[out]  mat The matrix.
 * \param[in]  dim Its dimension.
 * \retval none
 */
void centrosym_alloc(double **mat, int dim)
{

  *mat = (double*)calloc(dim * (dim+1) / 2, sizeof(double));

}

/*!
 * Return index for the (i,j)th elements of a centrosymmetric square matrix (fast; must have j <= i).
 * Naive indexing : row by row
 *
 * \param[in]  i Row index.
 * \param[in]  j Column index.
 * \param[in]  dim Matrix dimension.
 * \retval The index of the j-th element of the i-th row.
 */
int centrosym_ind(int i, int j, int dim){  // NAIVE INDEXING: ROW by ROW

  return i * (i + 1) / 2 + j ;

}

/*!
 * Return index for the (i,j)th elements of a centrosymmetric square matrix (fast; must have j <= i).
 * Second indexing : diagonal by diagonal
 *
 * \param[in]  i Row index.
 * \param[in]  j Column index.
 * \param[in]  dim Matrix dimension.
 * \retval The index of the j-th element of the i-th row.
 */
int centrosym_ind2(int i, int j, int dim){  // SECOND INDEXING: DIAG by DIAG

  int N = dim - 2;
  int k = i - j;

  return 2 * i - j + N * k - k * (k + 1) / 2 + k;

}


/*!
 * Create random centrosymmetric square matrix (full form, full square matrix).
 *
 * \param[out]  mat The matrix.
 * \param[in]  dim Its dimension.
 * \retval none
 */
void centrosym_full_random(double *mat, int dim)
{

  const int seed = (int)(10000.0*(double)clock()/(double)CLOCKS_PER_SEC);
  int i, j;
  double val;
  for(i = 0; i < dim; i++){
    for(j = 0; j <= i; j++){
      val = ran2_dp(seed);
      mat[ square_ind(i,j,dim) ] = val;
      mat[ square_ind(dim-i-1,dim-j-1,dim) ] = val ;
    }
  }

}


/*!
 * Extract the compressed form of a centrosymmetric matrix from its square form.
 *
 * \param[out]  matcomp The matrix in compressed form (triangle).
 * \param[in]  matfull The matrix in full form (square matrix).
 * \param[in]  dim Its dimension.
 * \retval none
 */
void centrosym_full_extractcomp(double *matcomp, double *matfull, int dim)
{

  int i, j;
  for(i = 0; i < dim; i++)
    for(j = 0; j <= i; j++)
      matcomp[ centrosym_ind(i,j,dim) ] = matfull[ square_ind(i,j,dim) ];

}


/*!
 * Assert equality of two centrosymmetric matrices in compressed form.
 *
 * \param[in]  matcomp1 The first matrix.
 * \param[in]  matcomp2 The secont matrix.
 * \param[in]  dim Their dimensions.
 * \retval 1 if they are equal.
 */
int centrosym_assertequal(double *matcomp1, double *matcomp2, int dim)
{

  int i, j;
  for(i = 0; i < dim; i++)
    for(j = 0; j <= i; j++)
      if( abs(matcomp1[ centrosym_ind(i,j,dim) ] - matcomp2[ centrosym_ind(i,j,dim)]) > PRECISION  ){
        //printf("Problem : %2.1f - %2.1f = %2.1e\n",matcomp1[ centrosym_ind(i,j,dim) ],matcomp2[ centrosym_ind(i,j,dim) ], matcomp1[ centrosym_ind(i,j,dim) ] - matcomp2[ centrosym_ind(i,j,dim)]);
        return 0;
      }
  return 1;   

}


/*!
 * Print a centrosymmetric square matrix in compressed form.
 *
 * \param[in]  mat The matrix.
 * \param[in]  dim Its dimension.
 * \retval none
 */
void centrosym_print(double *mat, int dim)
{
  int i, j;
  for(i = 0; i < dim; i++){
    for(j = 0; j <= i; j++){
      printf(" %2.1f ", mat[ centrosym_ind(i,j,dim) ]);
    }
    printf("\n");
  }
}


/*!
 * Compute the product of two centrosymmetric square matrices in compressed form (v1, fast).
 *
 * \param[out]  outmat The resulting matrix.
 * \param[in]  mat1 The first matrix.
 * \param[in]  mat2 The second matrix.
 * \param[in]  dim Their dimensions.
 * \retval none
 */
void centrosym_product(double *outmat, double *mat1, double *mat2, int dim)
{

  int i, j, k;
  for(i = 0; i < dim; i++){

    for(j = 0; j <= i; j++)
      outmat[ centrosym_ind(i,j,dim) ] = 0.0;

    for(k = 0; k <= i; k++){

      for(j = 0; j <= k; j++)
        outmat[ centrosym_ind(i,j,dim) ] += 
          mat1[ centrosym_ind(i,k,dim) ] * mat2[ centrosym_ind(k,j,dim) ];

      for(j = k+1; j <= i; j++)
        outmat[ centrosym_ind(i,j,dim) ] += 
          mat1[ centrosym_ind(i,k,dim) ] * mat2[ centrosym_ind(dim-k-1,dim-j-1,dim) ];

    }

    for(k = i+1 ; k < dim; k++){

      for(j = 0; j <= k; j++)
        outmat[ centrosym_ind(i,j,dim) ] += 
          mat1[ centrosym_ind(dim-i-1,dim-k-1,dim) ] * mat2[ centrosym_ind(k,j,dim) ];

      for(j = k+1; j <= i; j++)
        outmat[ centrosym_ind(i,j,dim) ] += 
          mat1[ centrosym_ind(dim-i-1,dim-k-1,dim) ] * mat2[ centrosym_ind(dim-k-1,dim-j-1,dim) ];

    }  
  }

}


/*!
 * Compute the product of two centrosymmetric square matrices in compressed form (v2, slower).
 *
 * \param[out]  outmat The resulting matrix.
 * \param[in]  mat1 The first matrix.
 * \param[in]  mat2 The second matrix.
 * \param[in]  dim Their dimensions.
 * \retval none
 */
void centrosym_product2(double *outmat, double *mat1, double *mat2, int dim)
{
  int i, j, k;
  for(i = 0; i < dim; i++){
    for(j = 0; j <= i; j++){
      outmat[ centrosym_ind2(i,j,dim) ] = 0.0;
      for(k = 0; k < dim; k++){
        outmat[ centrosym_ind2(i,j,dim) ] += 
          mat1[ centrosym_ind2(i,k,dim) ] * mat2[ centrosym_ind2(k,j,dim) ];
      }
    }
  }
}


/*!
 * Check if a square matrix is a valid centrosymmetric matrix.
 *
 * \param[in]  mat The square matrix.
 * \param[in]  dim Its dimensions.
 * \retval 1 if valid centrosymmetric matrix.
 */
int centrosym_isvalid(double *mat, int dim)
{

  int i, j;
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      if( abs(mat[ square_ind(i,j,dim) ] - mat[ square_ind(dim-i-1,dim-j-1,dim) ]) > PRECISION ){
        //printf("Problem : %2.3e != %2.3e\n",mat[ square_ind(i,j,dim) ],mat[ square_ind(dim-i-1,dim-j-1,dim) ]);
        return 0;
      }
    }
  }
  return 1;

}


/*!
 * Compute the trace of a centrosymmetric square matrix in compressed form.
 *
 * \param[in]  mat The input matrix.
 * \param[in]  dim Its dimensions.
 * \retval The trace of the matrix.
 */
double centrosym_trace(double *mat, int dim)
{

  int i;
  double res = 0.0;
  for(i = 0; i < dim; i++){
     res += mat[ centrosym_ind(i,i,dim) ];
  }
  return res;

}


/*!
 * Compute the trace of the product of two centrosymmetric square matrices in compressed form (indirect, slowo).
 *
 * \param[in]  mat1 The first matrix.
 * \param[in]  mat2 The second matrix.
 * \param[in]  dim Their dimensions.
 * \retval The trace of mat1  *mat2.
 */
double centrosym_traceprod2(double *mat1, double *mat2, int dim)
{

  double *matprod;
  centrosym_alloc(&matprod, dim);
  centrosym_product(matprod, mat1, mat2, dim);
  double res = centrosym_trace(matprod, dim);
  free(matprod);
  return res;

}

/*!
 * Compute the trace of the product of two centrosymmetric square matrices in compressed form (direct, fast).
 *
 * \param[in]  mat1 The first matrix.
 * \param[in]  mat2 The second matrix.
 * \param[in]  dim Their dimensions.
 * \retval The trace of mat1  *mat2.
 */
double centrosym_traceprod(double *mat1, double *mat2, int dim)
{

  int i, j;
  double res = 0.0;
  for(i = 0; i < dim; i++){

    for(j = 0; j <= i; j++){
      res += mat1[ centrosym_ind(i,j,dim) ] * mat2[ centrosym_ind(dim-j-1,dim-i-1,dim) ];
    }
    for(j = i+1; j < dim; j++){
      res += mat1[ centrosym_ind(dim-i-1,dim-j-1,dim) ] * mat2[ centrosym_ind(j,i,dim) ];
    }

  }
  return res;

}


/*!
 * Compute the quadratic form of a centrosymmetric square matrix in compressed form and two vectors (x^t * A * y).
 *
 * \param[in]  x The first vector.
 * \param[in]  mat2 The square matrix.
 * \param[in]  y The second vector.
 * \param[in]  dim Their dimensions.
 * \retval The quadratic form (x^t * A * y).
 */
double centrosym_quadform(double *x, double *mat, double *y, int dim)
{

  int i, j;
  double res = 0.0;
  for(i = 0; i < dim; i++){

    for(j = 0; j <= i; j++){
      res += x[i] * y[j] * mat[ centrosym_ind(i,j,dim) ] ;
    }
    for(j = i+1; j < dim; j++){
      res += x[i] * y[j] * mat[ centrosym_ind(dim-i-1,dim-j-1,dim) ] ;
    }

  }
  return res;

}
