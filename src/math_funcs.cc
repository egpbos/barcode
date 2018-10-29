/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <cmath> // sqrt, pow
#include <algorithm> // nth_element, max_element, min_element

#include "convenience.h" // copyArray

#include "math_funcs.h"

////////////////////////////////////////////////////////////////
// General scalar math functions
////////////////////////////////////////////////////////////////

// Check for odd-ness (instead of even-ness).
// Taken from http://forums.devshed.com/software-design-43/quick-algorithm-to-determine-odd-even-numbers-29843.html
int odd(int inputval) {
  return inputval & 1;
}

real_prec factorial(int term) {
  real_prec out = 1.0;

  if (term > 0) {
    for (int i = 1; i <= term; i++) {
      out *= static_cast<real_prec>(i);
    }
  } else {
    out = static_cast<real_prec>(1.);
  }

  return out;
}

real_prec power_mean(real_prec x, real_prec y, real_prec p) {
  real_prec mean;
  if (p == 0) {
    mean = std::sqrt(x*y);
  } else {
    mean = std::pow( ( std::pow(x, p) + std::pow(y, p) )/2, 1/p);
  }
  return mean;
}




////////////////////////////////////////////////////////////////
// Array functions / operators
// TODO: when the new array class (vector derived) is done, make
//       versions of these functions that take such a vector as
//       argument (by reference).
////////////////////////////////////////////////////////////////

real_prec min_arr ( ULONG factor, real_prec *in ) {
  real_prec *firstn=in;
  real_prec *lastn=in+factor;
  real_prec minn;
  real_prec *min = std::min_element(firstn, lastn);
  minn = *min;
  return minn;
}

real_prec max_arr ( ULONG factor, real_prec *in ) {
  real_prec *firstn=in;
  real_prec *lastn=in+factor;
  real_prec max;
  real_prec *maxn = std::max_element(firstn, lastn);
  max =*maxn;
  return max;
}

real_prec mean_arr ( ULONG size, const real_prec *in ) {
  real_prec mean = 0.0;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:mean)
#endif // MULTITHREAD
  for(ULONG i=0; i<size; i++) {
    mean += in[i];
  }
  return mean/static_cast<real_prec>(size);
}

real_prec median_arr (ULONG size, real_prec *in) {
  auto *copy = new real_prec[size]; // EGP: suggested by Johan; C++ way of doing it
  copyArray(in, copy, size);

  real_prec median;

  if (odd(static_cast<int>(size))) {
    std::nth_element(copy, copy + size/2, copy + size);
    median = copy[size/2];
  } else {
    std::nth_element(copy, copy + size/2 - 1, copy + size);
    median = copy[size/2-1]/2.;
    median += *std::min_element(copy + size/2, copy + size)/2.;
  }

  delete [] copy;
  return median;
}

real_prec std_arr (ULONG size, real_prec *in) {
  real_prec s_sq = 0.;
  real_prec mean = mean_arr(size, in);
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:s_sq)
#endif // MULTITHREAD
  for (ULONG i = 0; i < size; i++) {
    s_sq += std::pow(in[i] - mean, 2);
  }
  s_sq /= static_cast<real_prec>(size-1); // "corrected sample standard deviation" (squared)

  return std::sqrt(s_sq);
}

void complexfactor_mult (ULONG factor, real_prec in_a, complex_prec *in_b, complex_prec *out ) {
  for (ULONG i = 0; i < factor; i++) {
    re(out[i]) = in_a*re(in_b[i]);
    im(out[i]) = in_a*im(in_b[i]);
  }
}
