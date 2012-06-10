#include "symtrx.h"


/*!
 * Allocate space for square matrix.
 *
 * \param[out]  mat The matrix.
 * \param[in]  dim Its dimension.
 * \retval none
 */
 void square_alloc(double **mat, int dim)
{

  *mat = (double*)calloc(dim * dim, sizeof(double));

}


/*!
 * Create random square matrix.
 *
 * \param[out]  mat The matrix.
 * \param[in]  dim Its dimension.
 * \retval none
 */
void square_random(double *mat, int dim)
{

  const int seed = (int) (10000.0 * (double) clock() / (double) CLOCKS_PER_SEC);
  int i;
  for (i = 0; i < dim * dim; i++)
    mat[i] = ran2_dp(seed);

}


/*!
 * Symmetrise a square matrix.
 *
 * \param[inout]  mat The matrix.
 * \param[in]  dim Its dimension.
 * \retval none
 */
void square_symmetrise(double *mat, int dim)
{

  int i,j;
  for (i = 0; i < dim; i++)
    for (j = i + 1; j < dim; j++)
      mat[ square_ind(j,i,dim) ] = mat[ square_ind(i,j,dim) ];
}


/*!
 * Return index for the (i,j)th elements of a square matrix.
 *
 * \param[in]  i Row index.
 * \param[in]  j Column index.
 * \param[in]  dim Matrix dimension.
 * \retval The index of the j-th element of the i-th row.
 */
int square_ind(int i, int j, int dim){
  return i * dim + j;
}


/*!
 * Print a square matrix.
 *
 * \param[in]  mat The matrix.
 * \param[in]  dim Its dimension.
 * \retval none
 */
void square_print(double *mat, int dim)
{

  int i, j;
  for (i = 0; i < dim; i++){
    for (j = 0; j < dim; j++){
      printf(" %2.3e ", mat[ square_ind(i,j,dim) ]);
    }
    printf("\n");
  }

}


/*!
 * Compute the product of two square matrices.
 *
 * \param[out]  outmat The resulting matrix.
 * \param[in]  mat1 The first matrix.
 * \param[in]  mat2 The second matrix.
 * \param[in]  dim Their dimensions.
 * \retval none
 */
void square_product(double *outmat, double *mat1, double *mat2, int dim)
{

  int i, j, k;
  for(i = 0; i < dim; i++){
    for (j = 0; j < dim; j++){
      outmat[ square_ind(i,j,dim) ] = 0.0;
      for (k = 0; k < dim; k++){
        outmat[ square_ind(i,j,dim) ] += 
          mat1[ square_ind(i,k,dim) ] * mat2[ square_ind(k,j,dim) ];
      }
    }
  }

}


/*!
 * Compute the trace of a square matrix.
 *
 * \param[in]  mat The input matrix.
 * \param[in]  dim Its dimensions.
 * \retval The trace of the matrix.
 */
double square_trace(double *mat, int dim)
{

  int i;
  double res = 0.0;
  for(i = 0; i < dim; i++){
    res += mat[ square_ind(i,i,dim) ];
  }
  return res;

}


/*!
 * Compute the trace of the product of two square matrices.
 *
 * \param[in]  mat1 The first matrix.
 * \param[in]  mat2 The second matrix.
 * \param[in]  dim Their dimensions.
 * \retval The trace of mat1  *mat2.
 */
double square_traceprod(double *mat1, double *mat2, int dim)
{

  int i, j;
  double res = 0;
  for (i = 0; i < dim; i++){
    for (j = 0; j < dim; j++){
      res += mat1[ square_ind(i,j,dim) ] * mat2[ square_ind(j,i,dim) ];
    }
  }
  return res;

}


/*!
 * Compute the quadratic form of a square matrix and two vectors (x^t * A * y).
 *
 * \param[in]  x The first vector.
 * \param[in]  mat2 The square matrix.
 * \param[in]  y The second vector.
 * \param[in]  dim Their dimensions.
 * \retval The quadratic form (x^t * A * y).
 */
double square_quadform(double *x, double *mat, double *y, int dim)
{

  int i, j;
  double res = 0.0;
  for(i = 0; i < dim; i++){
    for (j = 0; j < dim; j++){
      res += x[i] * y[j] * mat[ square_ind(i,j,dim) ] ;
    }
  }
  return res;

}
