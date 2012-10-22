// Symmetrix (SMTRX) package
// Copyright (C) 2012
// Boris Leistedt

#ifndef SQUARE
#define SQUARE

void square_alloc(double **mat, int dim);
void square_random(double *mat, int dim);
void square_symmetrise(double *mat, int dim);
int square_ind(int i, int j, int dim);
void square_print(double *mat, int dim);
void square_product(double *outmat, double *mat1, double *mat2, int dim);
double square_trace(double *mat, int dim);
double square_traceprod(double *mat1, double *mat2, int dim);
double square_quadform(double *x, double *mat, double *y, int dim);

#endif
