/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include "define_opt.h"
#include "fftw_array.h"

#include <cmath>

using namespace std;

real_prec rankorder_leclercq_ZA(real_prec delta_ZA)
{
  // threshold value between two approximation ranges
  real_prec delta_th = pow(0.610/0.371, 1./(1.752-1.424));
  real_prec delta_Nbody;
  if (delta_ZA < delta_th)
    delta_Nbody = 0.610 * pow(delta_ZA + 1, 1.424) - 1;
  else
    delta_Nbody = 0.371 * pow(delta_ZA + 1, 1.752) - 1;
  return delta_Nbody;
}

real_prec rankorder_leclercq_2LPT(real_prec delta_2LPT)
{
  // threshold value between two approximation ranges
  real_prec delta_th = pow(0.642/0.257, 1./(1.922-1.401));
  real_prec delta_Nbody;
  if (delta_2LPT < delta_th)
    delta_Nbody = 0.642 * pow(delta_2LPT + 1, 1.401) - 1;
  else
    delta_Nbody = 0.257 * pow(delta_2LPT + 1, 1.922) - 1;
  return delta_Nbody;
}

void rankorder_leclercq_ZA(real_prec *delta_ZA, real_prec *delta_Nbody, ULONG size)
{
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (long i = 0; i < static_cast<long>(size); i++)
    delta_Nbody[i] = rankorder_leclercq_ZA(delta_ZA[i]);
}

void rankorder_leclercq_2LPT(real_prec *delta_2LPT, real_prec *delta_Nbody, ULONG size)
{
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (long i = 0; i < static_cast<long>(size); i++)
    delta_Nbody[i] = rankorder_leclercq_2LPT(delta_2LPT[i]);
}
