#include "symtrx.h"

/*!
 * Create random vector.
 *
 * \param[out]  mat The vector.
 * \param[in]  dim Its dimension.
 * \retval none
 */
void vector_random(double *x, int dim)
{
  const int seed = (int)(10000.0*(double)clock()/(double)CLOCKS_PER_SEC);
  int i;
  for(i = 0; i < dim; i++)
    x[i] = ran2_dp(seed);
}