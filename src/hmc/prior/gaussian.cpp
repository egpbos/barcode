/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include "struct_hamil.h"
#include "fftw_array.h"
#include "HMC_help.h" // convolveInvCorrFuncWithSignal

#include "hmc/prior/gaussian.hpp"

void prior_gaussian_grad_log_prior(struct HAMIL_DATA *hd, real_prec *signal, real_prec *out) {
  // dPsiP/ds = S^-1 * s
  convolveInvCorrFuncWithSignal(hd, signal, out, hd->signal_PS);
}

real_prec prior_gaussian_log_prior(struct HAMIL_DATA *hd, real_prec *signal) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  fftw_array<real_prec> dummy(n->N);
  convolveInvCorrFuncWithSignal(hd, signal, dummy, hd->signal_PS);

  // Psi_prior = 1/2 [s_k] * IFT[ 1/P*FT[s_k] ]

  real_prec psi_prior=0.;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:psi_prior)
#endif // MULTITHREAD
  for(ULONG i=0;i<n->N;i++)
    psi_prior += num_0_5 * signal[i] * dummy[i];

  return(psi_prior);
}
