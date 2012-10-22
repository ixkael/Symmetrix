#include "symtrx.h"

#define PRECISION 1e-12
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) > (b) ? (b) : (a))


/*!
 * Compute the actual size of a bisymmetric square matrix.
 *
 * \param[in]  dim The dimensions.
 * \retval Its size in memory.
 */
int bisym_size(int dim){

  if( dim % 2 == 0 ) // even 
    return dim*(dim+2)/4;
  else{ //odd
    int n = floor(dim/2);
    return n*(n+1) + n + 1;
  }

}


/*!
 * Allocate space for bisymmetric square matrix.
 *
 * \param[out]  mat The matrix.
 * \param[in]  dim Its dimension.
 * \retval none
 */
void bisym_alloc(double **mat, int dim)
{

  *mat = (double*)calloc( bisym_size(dim), sizeof(double));

}


/*!
 * Return index for the (i,j)th elements of a bisymmetric square matrix (fast; must have j <= i).
 *
 * \param[in]  i Row index.
 * \param[in]  j Column index.
 * \param[in]  dim Matrix dimension.
 * \retval The index of the j-th element of the i-th row.
 */
int bisym_ind(int i, int j, int dim){

  /*
  int normalcount = i*(i+1)/2 + j;
  int activator = i / ( dim / 2 + dim % 2 + 1 );
  int numberofterms = i - ( dim / 2 + dim % 2 ) ;
  int firstterm = 1 + dim % 2 ;
  int lastterm = 2 * (i - ( dim / 2 + dim % 2 )) - 1 + dim % 2;
  int count =  (numberofterms * (firstterm + lastterm)) / 2 ;
  printf(" %i ", normalcount - activator*count);
  */
  if( i / ( dim / 2 + dim % 2 + 1 ) )
    return i*(i+1)/2 + j - (i - ( dim / 2 + dim % 2 )) * ( 2 * dim % 2 + 2 * (i - ( dim / 2 + dim % 2 )) );
  else
    return i*(i+1)/2 + j;

}



/*!
 * Create random bisymmetric square matrix (full form, full square matrix).
 *
 * \param[out]  mat The matrix.
 * \param[in]  dim Its dimension.
 * \retval none
 */
void bisym_full_random(double *mat, int dim)
{

  const int seed = (int)(10000.0*(double)clock()/(double)CLOCKS_PER_SEC);
  int i, j;
  double val;
  for(i = 0; i < dim; i++){
    for(j = 0; j <= MIN(i, dim-i-1) ; j++){
      val = ran2_dp(seed);
      mat[ square_ind(i,j,dim) ] = val;
      mat[ square_ind(dim-i-1,dim-j-1,dim) ] = val ;
      mat[ square_ind(j,i,dim) ] = val ;
      mat[ square_ind(dim-j-1,dim-i-1,dim) ] = val ;
    }
  }

}


/*!
 * Extract the compressed form of a bisymmetric matrix from its square form.
 *
 * \param[out]  matcomp The matrix in compressed form (triangle).
 * \param[in]  matfull The matrix in full form (square matrix).
 * \param[in]  dim Its dimension.
 * \retval none
 */
void bisym_full_extractcomp(double *matcomp, double *matfull, int dim)
{

  int i, j;
  for(i = 0; i < dim; i++)
    for(j = 0; j <= MIN(i, dim-i-1); j++)
      matcomp[ bisym_ind(i,j,dim) ] = matfull[ square_ind(i,j,dim) ];

}


/*!
 * Assert equality of two bisymmetric matrices in compressed form.
 *
 * \param[in]  matcomp1 The first matrix.
 * \param[in]  matcomp2 The secont matrix.
 * \param[in]  dim Their dimensions.
 * \retval 1 if they are equal.
 */
int bisym_assertequal(double *matcomp1, double *matcomp2, int dim)
{

  int i, j;
  for(i = 0; i < dim; i++)
    for(j = 0; j <= MIN(i, dim-i-1); j++)
      if( abs(matcomp1[ bisym_ind(i,j,dim) ] - matcomp2[ bisym_ind(i,j,dim)]) > PRECISION  ){
        //WARNING: printf("Problem : %2.3e - %2.3e = %2.3e\n",matcomp1[ centrosym_ind(i,j,dim) ],matcomp2[ centrosym_ind(i,j,dim) ], matcomp1[ centrosym_ind(i,j,dim) ] - matcomp2[ centrosym_ind(i,j,dim)]);
        return 0;
      }
  return 1;  

}


/*!
 * Print a bisymmetric square matrix in compressed form.
 *
 * \param[in]  mat The matrix.
 * \param[in]  dim Its dimension.
 * \retval none
 */
void bisym_print(double *mat, int dim)
{

  int i, j;
  for(i = 0; i < dim; i++){
    for(j = 0; j <= MIN(i, dim-i-1); j++){
      printf(" %2.3e ", mat[ bisym_ind(i,j,dim) ]);
    }
    printf("\n");
  }

}


/*!
 * Compute the product of two bisymmetric square matrices in compressed form (v1, fast).
 *
 * \param[out]  outmat The resulting matrix.
 * \param[in]  mat1 The first matrix.
 * \param[in]  mat2 The second matrix.
 * \param[in]  dim Their dimensions.
 * \retval none
 */
void bisym_product(double *outmat, double *mat1, double *mat2, int dim)
{

  int i, j, k;
  for(i = 0; i < dim; i++){

    for(j = 0; j <= MIN(i, dim-i-1); j++)
      outmat[ bisym_ind(i,j, dim) ] = 0.0;

    for(k = 0; k <= MIN(i, dim-i-1); k++){

      for(j = 0; j <= MIN(k, dim-k-1); j++)
        outmat[ bisym_ind(i,j,dim) ] += 
          mat1[ bisym_ind(i,k,dim) ] * mat2[ bisym_ind(k,j,dim) ];

      //...

    }

    //...
    
  }

}


/*!
 * Check if a square matrix is a valid bisymmetric matrix.
 *
 * \param[in]  mat The square matrix.
 * \param[in]  dim Its dimensions.
 * \retval 1 if valid bisymmetric matrix.
 */
int bisym_isvalid(double *mat, int dim)
{

  int i, j;
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      if( 
          abs(mat[ square_ind(i,j,dim) ] 
            - mat[ square_ind(dim-i-1,dim-j-1,dim) ]) > PRECISION 
          ||
          abs(mat[ square_ind(i,j,dim) ] 
            - mat[ square_ind(j,i,dim) ]) > PRECISION 
          ||
          abs(mat[ square_ind(i,j,dim) ] 
            - mat[ square_ind(dim-j-1,dim-i-1,dim) ]) > PRECISION 
          ){
        //WARNING: printf("Problem : %2.3e != %2.3e\n",mat[ square_ind(i,j,dim) ],mat[ square_ind(dim-i-1,dim-j-1,dim) ]);
        return 0;
      }
    }
  }
  return 1;

}