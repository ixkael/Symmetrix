#include "melissa.h"

#define PRECISION 1e-12
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) > (b) ? (b) : (a))

//
// General operations on matrices
//

void matrixops_symmetrise(double *mat, int dim)
{
  int i,j;
  for(i = 0; i < dim; i++)
    for(j = i+1; j < dim; j++)
      mat[ matrixops_indsq(j,i,dim) ] = mat[ matrixops_indsq(i,j,dim) ];
}

void matrixops_inverse(double* A, int N)
{
    int *IPIV = (int*)calloc(N+1, sizeof(int));
    int LWORK = N*N;
    double *WORK = (double*)calloc(LWORK, sizeof(double));
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    free(IPIV);
    free(WORK);
}


void matrixops_alloc(double **mat, int dim1, int dim2)
{
  *mat = (double*)calloc(dim1 * dim2, sizeof(double));
}

void matrixops_allocsq(double **mat, int dim)
{
  *mat = (double*)calloc(dim * dim, sizeof(double));
}

void matrixops_random(double *mat, int dim1, int dim2)
{
  const int seed = (int)(10000.0*(double)clock()/(double)CLOCKS_PER_SEC);
  int i;
  for(i = 0; i < dim1 * dim2; i++)
    mat[i] = ran2_dp(seed);
}

void matrixops_randomvector(double *x, int dim)
{
  const int seed = (int)(10000.0*(double)clock()/(double)CLOCKS_PER_SEC);
  int i;
  for(i = 0; i < dim; i++)
    x[i] = ran2_dp(seed);
}

int matrixops_ind(int i, int j, int dim1, int dim2){
  return i * dim2 + j;
}

int matrixops_indsq(int i, int j, int dim){
  return i * dim + j;
}

void matrixops_product(double *outmat, double *mat1, double *mat2, int dim1, int dim2, int dim3)
{
  int i, j, k;
  for(i = 0; i < dim1; i++){
    for(j = 0; j < dim3; j++){
      outmat[ matrixops_ind(i,j,dim1,dim3) ] = 0.0;
      for(k = 0; k < dim2; k++){
        outmat[ matrixops_ind(i,j,dim1,dim3) ] += 
          mat1[ matrixops_ind(i,k,dim1,dim2) ] * mat2[ matrixops_ind(k,j,dim2,dim3) ];
      }
    }
  }
}

void matrixops_productvec(double *outvec, double *mat, double *vec, int dim1, int dim2)
{
  int i, j, k;
  for(i = 0; i < dim1; i++){
    outvec[ i ] = 0.0;
    for(j = 0; j < dim2; j++){
      outvec[ i ] += mat[ matrixops_ind(i,j,dim1,dim2) ] * vec[ j ];
    }
  }
}

void matrixops_printsq(double *mat, int dim)
{
  int i, j;
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      printf(" %2.3e ", mat[ matrixops_indsq(i,j,dim) ]);
    }
    printf("\n");
  }
}
void matrixops_printvec(double *vec, int dim)
{
  int i, j;
  for(i = 0; i < dim; i++){
    printf(" %2.3e ", vec[ i ]);
  }
  printf("\n");
}

void matrixops_productsq(double *outmat, double *mat1, double *mat2, int dim)
{
  int i, j, k;
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      outmat[ matrixops_indsq(i,j,dim) ] = 0.0;
      for(k = 0; k < dim; k++){
        outmat[ matrixops_indsq(i,j,dim) ] += 
          mat1[ matrixops_indsq(i,k,dim) ] * mat2[ matrixops_indsq(k,j,dim) ];
      }
    }
  }
}

double matrixops_traceprodsq(double *mat1, double *mat2, int dim)
{
  int i, j;
  double res = 0;
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      res += mat1[ matrixops_indsq(i,j,dim) ] * mat2[ matrixops_indsq(j,i,dim) ];
    }
  }
  return res;
}

double matrixops_tracesq(double *mat, int dim)
{
  int i, j;
  double res = 0.0;
  for(i = 0; i < dim; i++){
    res += mat[ matrixops_indsq(i,i,dim) ];
  }
  return res;
}

double matrixops_quadform(double *x, double *mat, double *y, int dim)
{
  int i, j;
  double res = 0.0;
  for(i = 0; i < dim; i++){

    for(j = 0; j < dim; j++){
      res += x[i] * y[j] * mat[ matrixops_indsq(i,j,dim) ] ;
    }

  }
  return res;
}

//
// Operations on centrosymmetric matrices
//

long matrixops_centrosym_size(int dim)
{
  return dim * (dim+1) / 2;
}

void matrixops_centrosym_alloc(double **mat, int dim)
{
  *mat = (double*)calloc(dim * (dim+1) / 2, sizeof(double));
}

int matrixops_centrosym_ind2(int i, int j, int dim){
  if( j <= i )
    return i * (i + 1) / 2 + j ;
  else
    return matrixops_centrosym_ind(dim-i-1, dim-j-1);
}

int matrixops_centrosym_ind(int i, int j){
  return i * (i + 1) / 2 + j ;
}

void matrixops_centrosym_randomfull(double *mat, int dim)
{
  const int seed = (int)(10000.0*(double)clock()/(double)CLOCKS_PER_SEC);
  int i, j;
  double val;
  for(i = 0; i < dim; i++){
    for(j = 0; j <= i; j++){
      val = ran2_dp(seed);
      mat[ matrixops_indsq(i,j,dim) ] = val;
      mat[ matrixops_indsq(dim-i-1,dim-j-1,dim) ] = val ;
    }
  }
}

void matrixops_centrosym_extractcomp(double *matcomp, double *matfull, int dim)
{
  int i, j;
  for(i = 0; i < dim; i++)
    for(j = 0; j <= i; j++)
      matcomp[ matrixops_centrosym_ind(i,j) ] = matfull[ matrixops_indsq(i,j,dim) ];
}

int matrixops_centrosym_assertequal(double *matcomp1, double *matcomp2, int dim)
{
  int i, j;
  for(i = 0; i < dim; i++)
    for(j = 0; j <= i; j++)
      if( abs(matcomp1[ matrixops_centrosym_ind(i,j) ] - matcomp2[ matrixops_centrosym_ind(i,j)]) > PRECISION  ){
        //printf("Problem : %2.3e - %2.3e = %2.3e\n",matcomp1[ matrixops_centrosym_ind(i,j,dim) ],matcomp2[ matrixops_centrosym_ind(i,j,dim) ], matcomp1[ matrixops_centrosym_ind(i,j,dim) ] - matcomp2[ matrixops_centrosym_ind(i,j,dim)]);
        return 0;
      }
  return 1;    
}

void matrixops_centrosym_print(double *mat, int dim)
{
  int i, j;
  for(i = 0; i < dim; i++){
    for(j = 0; j <= i; j++){
      printf(" %2.3e ", mat[ matrixops_centrosym_ind(i,j) ]);
    }
    printf("\n");
  }
}

void matrixops_centrosym_product2(double *outmat, double *mat1, double *mat2, int dim)
{
  int i, j, k;
  for(i = 0; i < dim; i++){
    for(j = 0; j <= i; j++){
      outmat[ matrixops_centrosym_ind2(i,j,dim) ] = 0.0;
      for(k = 0; k < dim; k++){
        outmat[ matrixops_centrosym_ind2(i,j,dim) ] += 
          mat1[ matrixops_centrosym_ind2(i,k,dim) ] * mat2[ matrixops_centrosym_ind2(k,j,dim) ];
      }
    }
  }
}

void matrixops_centrosym_product(double *outmat, double *mat1, double *mat2, int dim)
{
  int i, j, k;
  for(i = 0; i < dim; i++){

    for(j = 0; j <= i; j++)
      outmat[ matrixops_centrosym_ind(i,j) ] = 0.0;

    for(k = 0; k <= i; k++){

      for(j = 0; j <= k; j++)
        outmat[ matrixops_centrosym_ind(i,j) ] += 
          mat1[ matrixops_centrosym_ind(i,k) ] * mat2[ matrixops_centrosym_ind(k,j) ];

      for(j = k+1; j <= i; j++)
        outmat[ matrixops_centrosym_ind(i,j) ] += 
          mat1[ matrixops_centrosym_ind(i,k) ] * mat2[ matrixops_centrosym_ind(dim-k-1,dim-j-1) ];

    }

    for(k = i+1 ; k < dim; k++){

      for(j = 0; j <= k; j++)
        outmat[ matrixops_centrosym_ind(i,j) ] += 
          mat1[ matrixops_centrosym_ind(dim-i-1,dim-k-1) ] * mat2[ matrixops_centrosym_ind(k,j) ];

      for(j = k+1; j <= i; j++)
        outmat[ matrixops_centrosym_ind(i,j) ] += 
          mat1[ matrixops_centrosym_ind(dim-i-1,dim-k-1) ] * mat2[ matrixops_centrosym_ind(dim-k-1,dim-j-1) ];

    }
    
  }
}

int matrixops_centrosym_isvalid(double *mat, int dim)
{
  int i, j;
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      if( abs(mat[ matrixops_indsq(i,j,dim) ] - mat[ matrixops_indsq(dim-i-1,dim-j-1,dim) ]) > PRECISION ){
        //printf("Problem : %2.3e != %2.3e\n",mat[ matrixops_indsq(i,j,dim) ],mat[ matrixops_indsq(dim-i-1,dim-j-1,dim) ]);
        return 0;
      }
    }
  }
  return 1;
}

double matrixops_centrosym_trace(double *mat, int dim)
{
  int i, j;
  double res = 0.0;
  for(i = 0; i < dim; i++){
     res += mat[ matrixops_centrosym_ind(i,i) ];
  }
  return res;
}

double matrixops_centrosym_traceprod(double *mat1, double *mat2, int dim)
{
  int i, j;
  double res = 0.0;
  for(i = 0; i < dim; i++){

    for(j = 0; j <= i; j++){
      res += mat1[ matrixops_centrosym_ind(i,j) ] * mat2[ matrixops_centrosym_ind(dim-j-1,dim-i-1) ];
    }
    for(j = i+1; j < dim; j++){
      res += mat1[ matrixops_centrosym_ind(dim-i-1,dim-j-1) ] * mat2[ matrixops_centrosym_ind(j,i) ];
    }

  }
  return res;
}

double matrixops_centrosym_quadform(double *x, double *mat, double *y, int dim)
{
  int i, j;
  double res = 0.0;
  for(i = 0; i < dim; i++){

    for(j = 0; j <= i; j++){
      res += x[i] * y[j] * mat[ matrixops_centrosym_ind(i,j) ] ;
    }
    for(j = i+1; j < dim; j++){
      res += x[i] * y[j] * mat[ matrixops_centrosym_ind(dim-i-1,dim-j-1) ] ;
    }

  }
  return res;
}

//
// Operations on bisymmetric matrices
//

int matrixops_bisym_getsize(dim){
  if( dim % 2 == 0 ) // even 
    return dim*(dim+2)/4;
  else{ //odd
    int n = floor(n/2);
    return n*(n+1) + n + 1;
  }
}

void matrixops_bisym_alloc(double **mat, int dim)
{
  *mat = (double*)calloc( matrixops_bisym_getsize(dim), sizeof(double));
}


int matrixops_bisym_ind(int i, int j, int dim){
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



void matrixops_bisym_randomfull(double *mat, int dim)
{
  const int seed = (int)(10000.0*(double)clock()/(double)CLOCKS_PER_SEC);
  int i, j;
  double val;
  for(i = 0; i < dim; i++){
    for(j = 0; j <= MIN(i, dim-i-1) ; j++){
      val = ran2_dp(seed);
      mat[ matrixops_indsq(i,j,dim) ] = val;
      mat[ matrixops_indsq(dim-i-1,dim-j-1,dim) ] = val ;
      mat[ matrixops_indsq(j,i,dim) ] = val ;
      mat[ matrixops_indsq(dim-j-1,dim-i-1,dim) ] = val ;
    }
  }
}

void matrixops_bisym_extractcomp(double *matcomp, double *matfull, int dim)
{
  int i, j;
  for(i = 0; i < dim; i++)
    for(j = 0; j <= MIN(i, dim-i-1); j++)
      matcomp[ matrixops_bisym_ind(i,j,dim) ] = matfull[ matrixops_indsq(i,j,dim) ];
}

int matrixops_bisym_assertequal(double *matcomp1, double *matcomp2, int dim)
{
  int i, j;
  for(i = 0; i < dim; i++)
    for(j = 0; j <= MIN(i, dim-i-1); j++)
      if( cabs(matcomp1[ matrixops_bisym_ind(i,j,dim) ] - matcomp2[ matrixops_bisym_ind(i,j,dim)]) > PRECISION  ){
        //WARNING: printf("Problem : %2.3e - %2.3e = %2.3e\n",matcomp1[ matrixops_centrosym_ind(i,j,dim) ],matcomp2[ matrixops_centrosym_ind(i,j,dim) ], matcomp1[ matrixops_centrosym_ind(i,j,dim) ] - matcomp2[ matrixops_centrosym_ind(i,j,dim)]);
        return 0;
      }
  return 1;    
}

void matrixops_bisym_print(double *mat, int dim)
{
  int i, j;
  for(i = 0; i < dim; i++){
    for(j = 0; j <= MIN(i, dim-i-1); j++){
      printf(" %2.3e ", mat[ matrixops_bisym_ind(i,j,dim) ]);
    }
    printf("\n");
  }
}

void matrixops_bisym_product(double *outmat, double *mat1, double *mat2, int dim)
{
  int i, j, k;
  for(i = 0; i < dim; i++){

    for(j = 0; j <= MIN(i, dim-i-1); j++)
      outmat[ matrixops_bisym_ind(i,j, dim) ] = 0.0;

    for(k = 0; k <= MIN(i, dim-i-1); k++){

      for(j = 0; j <= MIN(k, dim-k-1); j++)
        outmat[ matrixops_bisym_ind(i,j,dim) ] += 
          mat1[ matrixops_bisym_ind(i,k,dim) ] * mat2[ matrixops_bisym_ind(k,j,dim) ];

      //...

    }

    //...
    
  }
}

int matrixops_bisym_isvalid(double *mat, int dim)
{
  int i, j;
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      if( 
          cabs(mat[ matrixops_indsq(i,j,dim) ] 
            - mat[ matrixops_indsq(dim-i-1,dim-j-1,dim) ]) > PRECISION 
          ||
          cabs(mat[ matrixops_indsq(i,j,dim) ] 
            - mat[ matrixops_indsq(j,i,dim) ]) > PRECISION 
          ||
          cabs(mat[ matrixops_indsq(i,j,dim) ] 
            - mat[ matrixops_indsq(dim-j-1,dim-i-1,dim) ]) > PRECISION 
          ){
        //WARNING: printf("Problem : %2.3e != %2.3e\n",mat[ matrixops_indsq(i,j,dim) ],mat[ matrixops_indsq(dim-i-1,dim-j-1,dim) ]);
        return 0;
      }
    }
  }
  return 1;
}