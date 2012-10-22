// Symmetrix (SMTRX) package
// Copyright (C) 2012
// Boris Leistedt

#ifndef CENTROSYM
#define CENTROSYM

void centrosym_alloc(double **mat, int dim);
int centrosym_ind(int i, int j, int dim);
int centrosym_ind2(int i, int j, int dim);
void centrosym_full_random(double *mat, int dim);
void centrosym_full_extractcomp(double *matcomp, double *matfull, int dim);
int centrosym_assertequal(double *matcomp1, double *matcomp2, int dim);
void centrosym_print(double *mat, int dim);
void centrosym_product(double *outmat, double *mat1, double *mat2, int dim);
int centrosym_isvalid(double *mat, int dim);
double centrosym_trace(double *mat, int dim);
double centrosym_traceprod(double *mat1, double *mat2, int dim);
double centrosym_traceprod2(double *mat1, double *mat2, int dim);
double centrosym_quadform(double *x, double *mat, double *y, int dim);

#endif
