// Symmetrix (SMTRX) package
// Copyright (C) 2012
// Boris Leistedt

#include "symtrx.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

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
