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
#include "gradient.hpp" // gradfft
#include "convenience.h" // multiply_factor_array
#include "Lag2Eul.h"
#include "gsl/gsl_math.h" // gsl_pow_2
#include "SPH_kernel.hpp" // SPH_kernel_scale

#include "hmc/likelihood/gaussian_independent.hpp"

// NOTE:
// partial_f_delta_x_log_like is also used in Gaussian RSD likelihood!
// Take care when introducing rank ordering, then it might become different.
void gaussian_likelihood_partial_f_delta_x_log_like(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  // EGP: partial -log L/partial delta_x = (Lambda-nobs)/sigma**2
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for(ULONG i = 0; i < n->N; i++)
  {
    real_prec Lambda = hd->window[i] * hd->rho_c * std::pow(num_1 + hd->biasP * deltaX[i], hd->biasE);
    if ((hd->window[i] > 0.) && (Lambda > 0.0))
    {
      real_prec resid = hd->nobs[i] - Lambda;
      out[i] = resid / (hd->noise[i] * hd->noise[i]);
    }
    else
      out[i]=0.0;
  }
}
// N.B.: grad_f_delta_x_comp should no longer be used. Only necessary for calc_h != 2, which is bogus!
void gaussian_likelihood_grad_f_delta_x_comp(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out,
                                             unsigned int component)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  //gradfindif(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, deltaX, out, component);
  gradfft(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, deltaX, out, component);
}
real_prec gaussian_likelihood_log_like(struct HAMIL_DATA *hd, real_prec *deltaQ)
{
  // -log L = 0.5*(Lambda-nobs)^2/sigma^2
  struct HAMIL_NUMERICAL *n = hd->numerical;
  fftw_array<real_prec> delta_growing(n->N);

  // single out growing mode (Peebles) -> delta_growing
  multiply_factor_array(n->deltaQ_factor, deltaQ, delta_growing, n->N);

  // calculate deltaX and pos
  {
    unsigned facL = 1;
    bool reggrid = true;
    gsl_rng *seed = nullptr; // empty: reggrid is true anyway
    real_prec kernel_scale = SPH_kernel_scale(hd);
    if (hd->rsd_model)
      Lag2Eul_rsd_zeldovich(delta_growing, hd->deltaX, hd->posx, hd->posy, hd->posz, n->N1, n->N2, n->N3, n->L1, n->L2,
                            n->L3, n->d1, n->d2, n->d3, n->min1, n->min2, n->min3, hd->D1, hd->ascale, hd->OM, hd->OL,
                            n->mk, facL,
                            reggrid, seed, kernel_scale, n->xobs, n->yobs, n->zobs, n->planepar, n->periodic,
                            n->R2Cplan, n->C2Rplan);
    else
      Lag2Eul(delta_growing, hd->deltaX, hd->posx, hd->posy, hd->posz, n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1,
              n->d2, n->d3, n->min1, n->min2, n->min3, hd->D1, hd->D2, hd->ascale, hd->OM, hd->OL, hd->sfmodel, n->mk, n->kth,
              facL,
              reggrid, seed, "", kernel_scale, n->R2Cplan, n->C2Rplan);
  }

  // actual likelihood calculation
  real_prec out=0.;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:out)
#endif // MULTITHREAD
  for(ULONG i=0;i<n->N;i++)
  {
    real_prec Lambda = hd->window[i] * hd->rho_c * std::pow(num_1+hd->biasP*hd->deltaX[i],hd->biasE);
    if ((hd->window[i]>0.) && (Lambda > 0.0))
      out+=num_0_5 * gsl_pow_2( (Lambda - hd->nobs[i]) / hd->noise[i] );
  }

  return (out);
}
