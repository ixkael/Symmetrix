
#include "../modules/melissa.h"
#include <assert.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) > (b) ? (b) : (a))

void test_matrixops_centrosym(int NREPEAT, int dim)
{
  int res, irepeat;
  clock_t t1, t2;
  double tmean_product_full, tmean_product_comp;
  double tmean_traceprod_full, tmean_traceprod_comp, tmean_traceprodnaive_full, tmean_traceprodnaive_comp;
  double tmean_quadform_full, tmean_quadform_comp;

  printf("\n==============================================\n");
  printf("Testing properties of centrosymmetric matrices\n");
  printf("----------------------------------------------\n");
  printf("Testing matrix product and traceproduct");

  for( irepeat = 0; irepeat < NREPEAT; irepeat++ ){
    fflush(NULL);
    printf(".");

    // Generate two random centrosymmetric matrices (full form)
    double *matfull1;
    matrixops_allocsq(&matfull1, dim);
    matrixops_centrosym_randomfull(matfull1, dim);
    res = matrixops_centrosym_isvalid(matfull1, dim);
    if(res == 0) printf("matfull1 is not centrosymmetric\n");
    double *matfull2;
    matrixops_allocsq(&matfull2, dim);
    matrixops_centrosym_randomfull(matfull2, dim);
    res = matrixops_centrosym_isvalid(matfull2, dim);
    if(res == 0) printf("matfull2 is not centrosymmetric\n");

    // Extract compressed forms
    double *matcomp1;
    matrixops_centrosym_alloc(&matcomp1, dim);
    matrixops_centrosym_extractcomp(matcomp1, matfull1, dim);
    double *matcomp2;
    matrixops_centrosym_alloc(&matcomp2, dim);
    matrixops_centrosym_extractcomp(matcomp2, matfull2, dim);

    // Perform the product in full form
    double *matfull3;
    matrixops_allocsq(&matfull3, dim);
    fflush(NULL); t1 = clock();
    matrixops_productsq(matfull3, matfull1, matfull2, dim);
    fflush(NULL); t2 = clock();
    tmean_product_full += (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    tmean_traceprodnaive_full += (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    //printf("Matrix product in full form : %4.4e seconds\n",(t2 - t1) / (double)CLOCKS_PER_SEC);
    double *matcomp3;
    matrixops_centrosym_alloc(&matcomp3, dim);
    matrixops_centrosym_extractcomp(matcomp3, matfull3, dim);

    // Check that the full form is still centrosymmetric
    res = matrixops_centrosym_isvalid(matfull3, dim);
    if(res == 0) printf("matfull3 is not centrosymmetric\n");

    // Perform the product in compressed fonm
    double *matcomp4;
    matrixops_centrosym_alloc(&matcomp4, dim);
    fflush(NULL); t1 = clock();
    matrixops_centrosym_product(matcomp4, matcomp1, matcomp2, dim);
    fflush(NULL); t2 = clock();
    tmean_product_comp += (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    tmean_traceprodnaive_comp += (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    //printf("Matrix product in comp form : %4.4e seconds\n",(t2 - t1) / (double)CLOCKS_PER_SEC);

    // Extract and compare the two results
    res = matrixops_centrosym_assertequal(matcomp4, matcomp3, dim);
    if(res == 0) printf("matcomp4 is not equal to matcomp3\n");

    // Test the trace of a product
    fflush(NULL); t1 = clock();
    double traceprodfull = matrixops_traceprodsq(matfull1, matfull2, dim);
    fflush(NULL); t2 = clock();
    tmean_traceprod_full += (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    fflush(NULL); t1 = clock();
    double traceprodcomp = matrixops_centrosym_traceprod(matcomp1, matcomp2, dim);
    fflush(NULL); t2 = clock();
    tmean_traceprod_comp += (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    if( cabs(traceprodfull - traceprodcomp) > 1e-10 ) printf("traceprodcomp is not equal to traceprodfull\n");

    fflush(NULL); t1 = clock();
    double traceprodfullnaive = matrixops_tracesq(matfull3, dim);
    fflush(NULL); t2 = clock();
    tmean_traceprodnaive_comp += (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    fflush(NULL); t1 = clock();
    double traceprodcompnaive = matrixops_centrosym_trace(matcomp4, dim);
    fflush(NULL); t2 = clock();
    tmean_traceprodnaive_full += (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    if( cabs(traceprodfullnaive - traceprodcompnaive) > 1e-10 ) printf("traceprodfullnaive is not equal to traceprodcompnaive\n");

    // Generate random vectors
    double *x = (double*)calloc(dim, sizeof(double));
    double *y = (double*)calloc(dim, sizeof(double));
    matrixops_randomvector(x, dim);
    matrixops_randomvector(y, dim);

    // Test quadratic forms
    fflush(NULL);t1 = clock();
    double quadform_full = matrixops_quadform(x, matfull3, y, dim);
    fflush(NULL);t2 = clock();
    tmean_quadform_full += (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    fflush(NULL);t1 = clock();
    double quadform_comp = matrixops_centrosym_quadform(x, matcomp4, y, dim);
    fflush(NULL);t2 = clock();
    tmean_quadform_comp += (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    //printf("\n %f - %f = %2.2e\n",quadform_full,quadform_comp,quadform_full-quadform_comp);
    if( cabs(quadform_full - quadform_comp) > 1e-6 ) printf("quadform_full is not equal to quadform_comp\n");


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
  //printf("Mean time for matrix product in full form : %4.4e seconds\n",tmean_product_full / (double)NREPEAT);
  //printf("Mean time for matrix product in comp form : %4.4e seconds\n",tmean_product_comp / (double)NREPEAT);
  printf("> Product acceleration factor : %2.2e \n", (double)tmean_product_full / (double)tmean_product_comp );
  printf("> Product storage gain factor : %2.2e \n", ((double)dim*dim)/(dim*(dim+1)/2) );
  //printf("Mean time for trace product in full form : %4.4e seconds\n",tmean_traceprod_full / (double)NREPEAT);
  //printf("Mean time for trace product in comp form : %4.4e seconds\n",tmean_traceprod_comp / (double)NREPEAT);
  //printf("Mean time for naive trace product in full form : %4.4e seconds\n",tmean_traceprodnaive_full / (double)NREPEAT);
  //printf("Mean time for naive trace product in comp form : %4.4e seconds\n",tmean_traceprodnaive_comp / (double)NREPEAT);
  printf("> Traceproduct acceleration factor : %2.2e \n", (double)tmean_traceprod_full / (double)tmean_traceprod_comp );
  //printf("Mean time for full quadratic form : %4.4e seconds\n",tmean_quadform_full / (double)NREPEAT);
  //printf("Mean time for compressed quadratic form : %4.4e seconds\n",tmean_quadform_comp / (double)NREPEAT);
  printf("> Quadratic form acceleration factor : %2.2e \n", (double)tmean_quadform_full / (double)tmean_quadform_comp );

  printf("----------------------------------------------\n");

}



void test_matrixops_bisym(int NREPEAT, int dim)
{
  int res, irepeat;
  clock_t t1, t2, tmean_product_full, tmean_product_comp;

  /*
  
  printf("\n==============================================\n");
  printf("Testing properties of bisymmetric matrices\n");
  printf("----------------------------------------------\n");

  for( irepeat = 0; irepeat < NREPEAT; irepeat++ ){
    fflush(NULL);
    printf(".");

    // Generate two random bisymmetric matrices (full form)
    double *matfull1;
    matrixops_allocsq(&matfull1, dim);
    matrixops_bisym_randomfull(matfull1, dim);
    res = matrixops_bisym_isvalid(matfull1, dim);
    if(res == 0) printf("matfull1 is not bisymmetric\n");
    double *matfull2;
    matrixops_allocsq(&matfull2, dim);
    matrixops_bisym_randomfull(matfull2, dim);
    res = matrixops_bisym_isvalid(matfull2, dim);
    if(res == 0) printf("matfull2 is not bisymmetric\n");

    // Extract compressed forms
    double *matcomp1;
    matrixops_bisym_alloc(&matcomp1, dim);
    matrixops_bisym_extractcomp(matcomp1, matfull1, dim);
    double *matcomp2;
    matrixops_bisym_alloc(&matcomp2, dim);
    matrixops_bisym_extractcomp(matcomp2, matfull2, dim);


    // Perform the product in full form
    double *matfull3;
    matrixops_allocsq(&matfull3, dim);
    t1 = clock();
    matrixops_productsq(matfull3, matfull1, matfull2, dim);
    t2 = clock();
    tmean_product_full += t2 - t1;
    //printf("Matrix product in full form : %4.4e seconds\n",(t2 - t1) / (double)CLOCKS_PER_SEC);
    double *matcomp3;
    matrixops_bisym_alloc(&matcomp3, dim);
    matrixops_bisym_extractcomp(matcomp3, matfull3, dim);

    // Check that the full form is still centrosymmetric
    res = matrixops_bisym_isvalid(matfull3, dim);
    if(res == 0) printf("matfull3 is not bisymmetric\n");

    // Perform the product in compressed fonm
    double *matcomp4;
    matrixops_bisym_alloc(&matcomp4, dim);
    t1 = clock();
    matrixops_bisym_product(matcomp4, matcomp1, matcomp2, dim);
    t2 = clock();
    tmean_product_comp += t2 - t1;
    //printf("Matrix product in comp form : %4.4e seconds\n",(t2 - t1) / (double)CLOCKS_PER_SEC);

    // Extract and compare the two results
    res = matrixops_bisym_assertequal(matcomp4, matcomp3, dim);
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
  printf("Mean time for matrix product in full form : %4.4e seconds\n",tmean_product_full / (double)CLOCKS_PER_SEC);
  printf("Mean time for matrix product in comp form : %4.4e seconds\n",tmean_product_comp / (double)CLOCKS_PER_SEC);
  printf("> Acceleration factor : %f \n", (double)tmean_product_full / tmean_product_comp );
  printf("> Storage size factor : %f \n", ((double)dim*dim)/(dim*(dim+1)/2) );

  printf("----------------------------------------------\n");

  */

}


void test_legendrematrices(double thetacut, long nside)
{

  printf("\n==============================================\n");
  printf("Testing properties of legendre matrices\n");
  printf("----------------------------------------------\n");
  printf("Testing angular maps with thetacut = %f\n", thetacut);

  int res;
  long npixtot = 12 * nside * nside;
  // Creating mask with azimutal cut
  int *mask = (int*)calloc(npixtot, sizeof(mask));
  healpix_azimask(mask, nside, thetacut);
  long npix = count(mask, 12*nside*nside);

  // Create angular map
  double *angmapfull = (double*)calloc(npix * npix, sizeof(double));
  healpix_angmap_full(angmapfull, mask, nside);


  // Allocate legendre matrices (full and compressed)
  double *Plfull = (double*)calloc(npix * npix, sizeof(double));
  double *Plcomp;
  matrixops_centrosym_alloc(&Plcomp, npix);

  // Create some legendre matrix
  const int ell = 4;
  legendre_matrix_full(Plfull, mask, nside, ell);
  legendre_matrix_comp(Plcomp, mask, nside, ell);

  res = matrixops_centrosym_isvalid(Plfull, npix);
  if(res == 0) printf("matfull3 is not centrosymmetric\n");

  double *matcomp3;
  matrixops_centrosym_alloc(&matcomp3, npix);
  matrixops_centrosym_extractcomp(matcomp3, Plfull, npix);
  res = matrixops_centrosym_assertequal(Plcomp, matcomp3, npix);
  if(res == 0) printf("Plcomp is not equal to Plfull\n");

  free(Plcomp);
  free(Plfull);

}
 

void test_clestimators(thetacut, nside, ellmax){

  // Create mask
  long npixtot = 12*nside*nside;
  int *mask = (int*)calloc(npixtot, sizeof(mask));
  healpix_azimask(mask, nside, thetacut);
  long npix = count(mask, 12*nside*nside);

  // Create map
  double *map = (double*)calloc(npixtot, sizeof(double));
  matrixops_randomvector(map, npixtot);

  // Perform PCL estimation (full)
  double *cl_pcl_full = (double*)calloc(ellmax+1, sizeof(double));
  map2cl_PCL_full(cl_pcl_full, map, nside, ellmax, mask);
  matrixops_printvec(cl_pcl_full);

}


int main(int argc, char *argv[]) 
{

  const int dim = 100;
  const long nside = 4;
  const int NREPEAT = 5;
  const int ellmax = 8;
  const double thetacut = 60.0;

  printf("\n==============================================\n");
  printf("MELISSA TESTS\n");
  printf("----------------------------------------------\n");
  printf("\n");
  printf("HEALPIX nside parameter : %li\n", nside);
  printf("Maximum multipole elmax : %i\n", ellmax);
  printf("Size of square matrices : %i\n", dim);
  printf("Number of realisations  : %i\n", NREPEAT);

  // Testing centrosymmetric matrices
  test_matrixops_centrosym(NREPEAT, dim);

  // Testing bisymmetric matrices  
  test_legendrematrices(thetacut, nside);

  // Testing power spectrum estimators
  test_clestimators(thetacut, nside, ellmax);

  printf("\n");
  return 1;
}