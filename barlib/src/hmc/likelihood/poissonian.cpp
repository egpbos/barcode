/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <cmath> // pow, log

#include "struct_hamil.h"
#include "gradient.hpp" // gradfindif
#include "Lag2Eul.h"
#include "SPH_kernel.hpp" // SPH_kernel_scale

#include "hmc/likelihood/poissonian.hpp"

// Model: Poissonian likelihood //
void poissonian_likelihood_partial_f_delta_x_log_like(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *dummy)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;

  for(ULONG i = 0; i < n->N; i++)
  {
    auto dens = static_cast<real_prec>(1. + hd->biasP * deltaX[i]);// must be positive!
    real_prec Lambda = hd->window[i] * hd->rho_c * std::pow(dens, hd->biasE);

    if ((hd->window[i] > 0.0) && (dens > 0.0))
      //dummy[i] = (Lambda - hd->nobs[i]) * (hd->biasE * hd->biasP / dens);
      dummy[i] = (1 - hd->nobs[i]/Lambda) * hd->rho_c * hd->biasE * hd->biasP * std::pow(dens, hd->biasE - 1);
    else
      dummy[i] = 0.0;
  }
}

// N.B.: grad_f_delta_x_comp should no longer be used. Only necessary for calc_h != 2, which is bogus!
void poissonian_likelihood_grad_f_delta_x_comp(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out,
                                               unsigned int component)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  gradfindif(n->N1, n->L1, deltaX, out, component);
}

real_prec poissonian_likelihood_log_like(struct HAMIL_DATA *hd, real_prec *delta)
{
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
  real_prec out = 0.;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:out)
#endif // MULTITHREAD
  for(ULONG i = 0; i < n->N; i++)
  {
    auto dens = static_cast<real_prec>(1. + hd->biasP * hd->deltaX[i]);// must be positive!
    real_prec Lambda = hd->window[i] * hd->rho_c * std::pow(dens, hd->biasE);

    if ((hd->window[i] > 0.) && (Lambda > 0.0))
      out += Lambda - hd->nobs[i] * std::log(Lambda);
  }

  return (out);
}

// END model Poissonian likelihood //
