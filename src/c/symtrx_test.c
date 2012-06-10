
#include "symtrx.h"
#include <assert.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) > (b) ? (b) : (a))

void test_centrosym(int NREPEAT, int dim)
{
  int res, irepeat;
  clock_t t1, t2;
  double tmean_product_full=0, tmean_product_comp=0;
  double tmean_traceprod_full=0, tmean_traceprod_comp=0;
  double tmean_traceprodnaive_full=0, tmean_traceprodnaive_comp=0;
  double tmean_quadform_full=0, tmean_quadform_comp=0;

  printf("\n==============================================\n");
  printf("Testing properties of centrosymmetric matrices\n");
  printf("----------------------------------------------\n");
  printf("Performing benchmark");

  for( irepeat = 0; irepeat < NREPEAT; irepeat++ ){
    fflush(NULL);
    printf(".");

    // Generate two random centrosymmetric matrices (full form)
    double *matfull1;
    square_alloc(&matfull1, dim);
    centrosym_full_random(matfull1, dim);
    res = centrosym_isvalid(matfull1, dim);
    if(res == 0) printf("matfull1 is not centrosymmetric\n");
    double *matfull2;
    square_alloc(&matfull2, dim);
    centrosym_full_random(matfull2, dim);
    res = centrosym_isvalid(matfull2, dim);
    if(res == 0) printf("matfull2 is not centrosymmetric\n");

    // Extract compressed forms
    double *matcomp1;
    centrosym_alloc(&matcomp1, dim);
    centrosym_full_extractcomp(matcomp1, matfull1, dim);
    double *matcomp2;
    centrosym_alloc(&matcomp2, dim);
    centrosym_full_extractcomp(matcomp2, matfull2, dim);

    // Perform the product in full form
    double *matfull3;
    square_alloc(&matfull3, dim);
    fflush(NULL); t1 = clock();
    square_product(matfull3, matfull1, matfull2, dim);
    fflush(NULL); t2 = clock();
    tmean_product_full += (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    tmean_traceprodnaive_full += (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    //printf("Matrix product in full form : %4.4e seconds\n",(t2 - t1) / (double)CLOCKS_PER_SEC);
    double *matcomp3;
    centrosym_alloc(&matcomp3, dim);
    centrosym_full_extractcomp(matcomp3, matfull3, dim);

    // Check that the full form is still centrosymmetric
    res = centrosym_isvalid(matfull3, dim);
    if(res == 0) printf("matfull3 is not centrosymmetric\n");

    // Perform the product in compressed fonm
    double *matcomp4;
    centrosym_alloc(&matcomp4, dim);
    fflush(NULL); t1 = clock();
    centrosym_product(matcomp4, matcomp1, matcomp2, dim);
    fflush(NULL); t2 = clock();
    tmean_product_comp += (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    tmean_traceprodnaive_comp += (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    //printf("Matrix product in comp form : %4.4e seconds\n",(t2 - t1) / (double)CLOCKS_PER_SEC);

    // Extract and compare the two results
    res = centrosym_assertequal(matcomp4, matcomp3, dim);
    if(res == 0) printf("matcomp4 is not equal to matcomp3\n");
    //centrosym_print(matcomp3, dim);
    //centrosym_print(matcomp4, dim);

    // Test the trace of a product
    fflush(NULL); t1 = clock();
    double traceprodfull = square_traceprod(matfull1, matfull2, dim);
    fflush(NULL); t2 = clock();
    tmean_traceprod_full += (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    fflush(NULL); t1 = clock();
    double traceprodcomp = centrosym_traceprod(matcomp1, matcomp2, dim);
    fflush(NULL); t2 = clock();
    tmean_traceprod_comp += (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    if( abs(traceprodfull - traceprodcomp) > 1e-10 ) printf("traceprodcomp is not equal to traceprodfull\n");

    fflush(NULL); t1 = clock();
    double traceprodfullnaive = square_trace(matfull3, dim);
    fflush(NULL); t2 = clock();
    tmean_traceprodnaive_comp += (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    fflush(NULL); t1 = clock();
    double traceprodcompnaive = centrosym_trace(matcomp4, dim);
    fflush(NULL); t2 = clock();
    tmean_traceprodnaive_full += (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    if( abs(traceprodfullnaive - traceprodcompnaive) > 1e-10 ) printf("traceprodfullnaive is not equal to traceprodcompnaive\n");

    // Generate random vectors
    double *x = (double*)calloc(dim, sizeof(double));
    double *y = (double*)calloc(dim, sizeof(double));
    vector_random(x, dim);
    vector_random(y, dim);

    // Test quadratic forms
    fflush(NULL);t1 = clock();
    double quadform_full = square_quadform(x, matfull3, y, dim);
    fflush(NULL);t2 = clock();
    tmean_quadform_full += (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    fflush(NULL);t1 = clock();
    double quadform_comp = centrosym_quadform(x, matcomp4, y, dim);
    fflush(NULL);t2 = clock();
    tmean_quadform_comp += (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    //printf("\n %f - %f = %2.2e\n",quadform_full,quadform_comp,quadform_full-quadform_comp);
    if( abs(quadform_full - quadform_comp) > 1e-6 ) printf("quadform_full is not equal to quadform_comp\n");


    free(x);
    free(y);
    free(matfull1);
    free(matcomp1);
    free(matfull2);
    free(matcomp2);
    free(matfull3);
    free(matcomp3);
    free(matcomp4);

  }

  printf("done\n");

  printf("> Memory gain due to compressed form : %2.2f \n", ((double)dim*dim)/(dim*(dim+1)/2) );

  //printf("Mean time for matrix product in full form : %4.4f seconds\n",tmean_product_full / (double)NREPEAT);
  //printf("Mean time for matrix product in comp form : %4.4f seconds\n",tmean_product_comp / (double)NREPEAT);
  printf("> Matrix product acceleration factor : %2.2f \n", (double)tmean_product_full / (double)tmean_product_comp );

  //printf("Mean time for trace product in full form : %4.4f seconds\n",tmean_traceprod_full / (double)NREPEAT);
  //printf("Mean time for trace product in comp form : %4.4f seconds\n",tmean_traceprod_comp / (double)NREPEAT);
  //printf("Mean time for naive trace product in full form : %4.4f seconds\n",tmean_traceprodnaive_full / (double)NREPEAT);
  //printf("Mean time for naive trace product in comp form : %4.4f seconds\n",tmean_traceprodnaive_comp / (double)NREPEAT);
  printf("> Trace-product acceleration factor  : %2.2f \n", (double)tmean_traceprod_full / (double)tmean_traceprod_comp );
  
  //printf("Mean time for full quadratic form : %4.4f seconds\n",tmean_quadform_full / (double)NREPEAT);
  //printf("Mean time for compressed quadratic form : %4.4f seconds\n",tmean_quadform_comp / (double)NREPEAT);
  printf("> Quadratic form acceleration factor : %2.2f \n", (double)tmean_quadform_full / (double)tmean_quadform_comp );

  printf("----------------------------------------------\n");

}



void test_bisym(int NREPEAT, int dim)
{
  int res, irepeat;
  clock_t t1, t2;
  double tmean_product_full=0, tmean_product_comp=0;

  
  
  printf("\n==============================================\n");
  printf("Testing properties of bisymmetric matrices\n");
  printf("----------------------------------------------\n");

  for( irepeat = 0; irepeat < NREPEAT; irepeat++ ){
    fflush(NULL);
    printf(".");

    // Generate two random bisymmetric matrices (full form)
    double *matfull1;
    square_alloc(&matfull1, dim);
    bisym_full_random(matfull1, dim);
    res = bisym_isvalid(matfull1, dim);
    if(res == 0) printf("matfull1 is not bisymmetric\n");
    double *matfull2;
    square_alloc(&matfull2, dim);
    bisym_full_random(matfull2, dim);
    res = bisym_isvalid(matfull2, dim);
    if(res == 0) printf("matfull2 is not bisymmetric\n");

    // Extract compressed forms
    double *matcomp1;
    bisym_alloc(&matcomp1, dim);
    bisym_full_extractcomp(matcomp1, matfull1, dim);
    double *matcomp2;
    bisym_alloc(&matcomp2, dim);
    bisym_full_extractcomp(matcomp2, matfull2, dim);


    // Perform the product in full form
    double *matfull3;
    square_alloc(&matfull3, dim);
    t1 = clock();
    square_product(matfull3, matfull1, matfull2, dim);
    t2 = clock();
    tmean_product_full += t2 - t1;
    //printf("Matrix product in full form : %4.4e seconds\n",(t2 - t1) / (double)CLOCKS_PER_SEC);
    double *matcomp3;
    bisym_alloc(&matcomp3, dim);
    bisym_full_extractcomp(matcomp3, matfull3, dim);

    // Check that the full form is still centrosymmetric
    res = bisym_isvalid(matfull3, dim);
    if(res == 0) printf("matfull3 is not bisymmetric\n");

    // Perform the product in compressed fonm
    double *matcomp4;
    bisym_alloc(&matcomp4, dim);
    t1 = clock();
    bisym_product(matcomp4, matcomp1, matcomp2, dim);
    t2 = clock();
    tmean_product_comp += t2 - t1;
    //printf("Matrix product in comp form : %4.4e seconds\n",(t2 - t1) / (double)CLOCKS_PER_SEC);

    // Extract and compare the two results
    res = bisym_assertequal(matcomp4, matcomp3, dim);
    if(res == 0) printf("matcomp4 is not equal to matcomp3\n");

    free(matfull1);
    free(matcomp1);
    free(matfull2);
    free(matcomp2);
    free(matfull3);
    free(matcomp3);
    free(matcomp4);

  }

  printf("done\n");
  printf("Mean time for matrix product in full form : %4.4f seconds\n",tmean_product_full / (double)CLOCKS_PER_SEC);
  printf("Mean time for matrix product in comp form : %4.4f seconds\n",tmean_product_comp / (double)CLOCKS_PER_SEC);
  printf("> Acceleration factor : %f \n", (double)tmean_product_full / tmean_product_comp );
  printf("> Storage size factor : %f \n", ((double)dim*dim)/(dim*(dim+1)/2) );

  printf("----------------------------------------------\n");

  

}

 
int main(int argc, char *argv[]) 
{

  const int dim = 6;
  const int NREPEAT = 5;

  /*
  int i, j;
  double *mat;
  centrosym_alloc(&mat, dim);
  for(i=0;i<dim;i++)
    for(j=0;j<dim;j++)
      mat[centrosym_ind(i,j,dim)] = centrosym_ind(i,j,dim);

  centrosym_print(mat, dim);
  */

  printf("\n==============================================\n");
  printf("SYMTRX TESTS\n");
  printf("----------------------------------------------\n");
  printf("\n");
  printf("Size of square matrices : %i x %i\n", dim, dim);
  printf("Number of realisations  : %i\n", NREPEAT);

  // Testing centrosymmetric matrices
  test_centrosym(NREPEAT, dim);

  // Testing bisymmetric matrices
  test_bisym(NREPEAT, dim);

  printf("\n");
  return 0;
}