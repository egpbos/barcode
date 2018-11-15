/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <gsl/gsl_math.h> // gsl_pow_2

#include "struct_hamil.h"

#include "hmc/likelihood/gaussian_random_field.hpp"

// Model: Gaussian random field (no evolution, just lagrangian field) //
// N.B.: grad_f_delta_x_comp should no longer be used. Only necessary for calc_h != 2, which is bogus!
void grf_likelihood_grad_f_delta_x_comp(struct HAMIL_DATA *, real_prec *, real_prec *, unsigned int)
{
}

void grf_likelihood_partial_f_delta_x_log_like(struct HAMIL_DATA *, real_prec *, real_prec *)
{
}

void grf_likelihood_grad_log_like(struct HAMIL_DATA *hd, real_prec *delta, real_prec *out)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for(ULONG i=0;i<n->N;i++)
    if (hd->window[i]>0.)
      out[i] = (delta[i] - hd->nobs[i]) / (hd->noise[i] * hd->noise[i]);
    else
      out[i] = 0;
}

real_prec grf_likelihood_log_like(struct HAMIL_DATA *hd, real_prec *delta)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;

  real_prec out=0.;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:out)
#endif // MULTITHREAD
  for(ULONG i=0;i<n->N;i++)
    if (hd->window[i]>0.)
      out += num_0_5 * gsl_pow_2((delta[i] - hd->nobs[i])/hd->noise[i]);

  return (out);
}

// END model Gaussian random field //
