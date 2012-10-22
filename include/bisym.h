// Symmetrix (SMTRX) package
// Copyright (C) 2012
// Boris Leistedt

#ifndef BISYM
#define BISYM

int bisym_size(int dim);
void bisym_alloc(double **mat, int dim);
int bisym_ind(int i, int j, int dim);
void bisym_full_random(double *mat, int dim);
void bisym_full_extractcomp(double *matcomp, double *matfull, int dim);
int bisym_assertequal(double *matcomp1, double *matcomp2, int dim);
void bisym_print(double *mat, int dim);
void bisym_product(double *outmat, double *mat1, double *mat2, int dim);
int bisym_isvalid(double *mat, int dim);

#endif
