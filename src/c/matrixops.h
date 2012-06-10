
#ifndef MATRIXOPS
#define MATRIXOPS

void matrixops_symmetrise(double *F, int dim);
void matrixops_inverse(double* A, int N);
long matrixops_centrosym_size(int dim);
void matrixops_alloc(double **mat, int dim1, int dim2);
void matrixops_allocsq(double **mat, int dim);
void matrixops_random(double *mat, int dim1, int dim2);
void matrixops_randomvector(double *x, int dim);
int matrixops_ind(int i, int j, int dim1, int dim2);
int matrixops_indsq(int i, int j, int dim);
void matrixops_product(double *outmat, double *mat1, double *mat2, int dim1, int dim2, int dim3);
void matrixops_printsq(double *mat, int dim);
void matrixops_productsq(double *outmat, double *mat1, double *mat2, int dim);
double matrixops_traceprodsq(double *mat1, double *mat2, int dim);
double matrixops_tracesq(double *mat, int dim);
double matrixops_quadform(double *x, double *mat, double *y, int dim);

void matrixops_centrosym_alloc(double **mat, int dim);
int matrixops_centrosym_ind(int i, int j);
int matrixops_centrosym_ind2(int i, int j, int dim);
void matrixops_centrosym_randomfull(double *mat, int dim);
void matrixops_centrosym_extractcomp(double *matcomp, double *matfull, int dim);
int matrixops_centrosym_assertequal(double *matcomp1, double *matcomp2, int dim);
void matrixops_centrosym_print(double *mat, int dim);
void matrixops_centrosym_product(double *outmat, double *mat1, double *mat2, int dim);
void matrixops_centrosym_product2(double *outmat, double *mat1, double *mat2, int dim);
int matrixops_centrosym_isvalid(double *mat, int dim);
double matrixops_centrosym_trace(double *mat, int dim);
double matrixops_centrosym_traceprod(double *mat1, double *mat2, int dim);
double matrixops_centrosym_quadform(double *x, double *mat, double *y, int dim);

int matrixops_bisym_getsize(dim);
void matrixops_bisym_alloc(double **mat, int dim);
int matrixops_bisym_ind(int i, int j, int dim);
void matrixops_bisym_randomfull(double *mat, int dim);
void matrixops_bisym_extractcomp(double *matcomp, double *matfull, int dim);
int matrixops_bisym_assertequal(double *matcomp1, double *matcomp2, int dim);
void matrixops_bisym_print(double *mat, int dim);
void matrixops_bisym_product(double *outmat, double *mat1, double *mat2, int dim);
int matrixops_bisym_isvalid(double *mat, int dim);

void packvector(double *x, double *map, int *mask, int dim);

#endif