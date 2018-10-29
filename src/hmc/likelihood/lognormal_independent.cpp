/*
 * Barcode
 * Copyright E.G.P. Bos
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <cmath>

#include "struct_hamil.h"
#include "fftw_array.h"
#include "gradient.hpp" // gradfindif
#include "Lag2Eul.h"
#include "SPH_kernel.hpp" // SPH_kernel_scale

#include "hmc/likelihood/lognormal_independent.hpp"

// Model: log-normal likelihood //

// TEST: with lambda function (c++11) for Lambda
//       This is also a useful structure for when/if we put these models in
//       full-fledged classes; one could then define the function for Lambda
//       as a class method and reuse it in the different functions.
//void lognormal_likelihood_partial_f_delta_x_log_like(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *dummy)
//{
//struct HAMIL_NUMERICAL *n = hd->numerical;
//// EGP: partial -log L/partial log(1+delta_x) = (Lambda-nobs)/sigma**2
//auto Lambda = [&] (ULONG index) -> real_prec {log(hd->rho_c * pow(num_1 + hd->biasP * deltaX[index], hd->biasE));};
//#ifdef MULTITHREAD
//#pragma omp parallel for
//#endif // MULTITHREAD
//for(ULONG i=0; i < n->N;i++)
//if (hd->window[i]>0.)
//dummy[i] = (hd->nobs[i] - Lambda(i))/(hd->noise[i] * hd->noise[i]);
//else
//dummy[i] = 0.0;
//}

void lognormal_likelihood_partial_f_delta_x_log_like(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *dummy)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  // EGP: partial -log L/partial log(1+delta_x) = (Lambda-nobs)/sigma**2
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for(ULONG i=0; i < n->N;i++)
  {
    real_prec Lambda = log(hd->rho_c * std::pow(num_1 + hd->biasP * deltaX[i], hd->biasE));
    if (hd->window[i]>0.)
      dummy[i]=(hd->nobs[i] - Lambda)/(hd->noise[i] * hd->noise[i]);
    else
      dummy[i]=0.0;
  }
}

real_prec lognormal_likelihood_f_delta_x_i_calc(real_prec rho_c, real_prec delta_min, real_prec deltaX_i)
{
  // keep above zero density
  if (deltaX_i < delta_min)
    deltaX_i = delta_min;

  return std::log(rho_c * (num_1 + deltaX_i));
}

real_prec lognormal_likelihood_f_delta_x_i(struct HAMIL_DATA *hd, real_prec deltaX_i)
{
  return lognormal_likelihood_f_delta_x_i_calc(hd->rho_c, hd->delta_min, deltaX_i);
}

void lognormal_likelihood_f_delta_x(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (ULONG i = 0; i < n->N; i++)
    out[i] = lognormal_likelihood_f_delta_x_i(hd, deltaX[i]);
}

// N.B.: grad_f_delta_x_comp should no longer be used. Only necessary for calc_h != 2, which is bogus!
void lognormal_likelihood_grad_f_delta_x_comp(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out,
                                              unsigned int component)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  fftw_array<real_prec> f_delta_x(n->N);

  lognormal_likelihood_f_delta_x(hd, deltaX, f_delta_x);

  gradfindif(n->N1, n->L1, f_delta_x, out, component);
}

real_prec lognormal_likelihood_log_like(struct HAMIL_DATA *hd, real_prec *delta)
{
  // -log L = 0.5*(Lambda-nobs)^2/sigma^2
  struct HAMIL_NUMERICAL *n = hd->numerical;

  // calculate deltaX and pos
  {
    unsigned facL=1;
    bool reggrid=true;
    gsl_rng *seed = nullptr; // empty: reggrid is true anyway
    real_prec kernel_scale = SPH_kernel_scale(hd);
    Lag2Eul(delta, hd->deltaX, hd->posx, hd->posy, hd->posz, n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1, n->d2,
            n->d3, n->min1, n->min2, n->min3, hd->D1, hd->D2, hd->ascale, hd->OM, hd->OL, hd->sfmodel, n->mk, n->kth, facL,
            reggrid, seed, "", kernel_scale, n->R2Cplan, n->C2Rplan);
  }

  // actual likelihood calculation
  real_prec out=0.;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:out)
#endif // MULTITHREAD
  for(ULONG i=0;i<n->N;i++)
  {
    real_prec Lambda = lognormal_likelihood_f_delta_x_i(hd, hd->deltaX[i]);
    if (hd->window[i]>0.)
    {
      real_prec resid=Lambda - hd->nobs[i];
      out+=num_0_5*resid*resid/(hd->noise[i]*hd->noise[i]);
    }
  }

  return (out);
}
// End model: log-normal likelihood //
