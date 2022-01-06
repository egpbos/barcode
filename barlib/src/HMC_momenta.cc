/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include <fftw3.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>

#include <fstream>  //NOLINT

#include "struct_main.h"
#include "struct_hamil.h"
#include "fftw_array.h"
#include "math_funcs.h"
#include "convenience.h"
#include "random.hpp"

#include "debug.h"

using namespace std;

#ifdef MASKING
bool masking2 = true;
#else
bool masking2 = false;
#endif  // ifdef MASKING




// Momenta: Gaussian momentum distribution //
void draw_masked_momenta(struct HAMIL_DATA *hd, gsl_rng *seed,
  real_prec *momenta);  // helper function defined below
void draw_real_space_momenta(struct HAMIL_DATA *hd, gsl_rng *seed,
  real_prec *momenta);  // helper function defined below
void draw_momenta(struct HAMIL_DATA *hd, gsl_rng *seed, real_prec *momenta,
#ifdef DEBUG
    struct DATA *data
#else
    struct DATA *
#endif // DEBUG
      ) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  // p ~ exp[-p^2/(2*M)]

  if (n->mass_fs) {
    create_GARFIELD(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, momenta,
      hd->mass_f, seed);

    // N.B.: als we masking weer mee willen nemen moeten we onderstaande functie
    //       hebben
    // draw_masked_momenta(hd, seed, momenta);
  } else {
    fillZero(momenta, n->N);  // prepare for addition to real part
  }

  if (n->mass_rs) {
    fftw_array<real_prec> dummy(n->N);
    draw_real_space_momenta(hd, seed, dummy);
    add_to_array(dummy, momenta, n->N);
  }

#ifdef DEBUG
  debug_array_statistics(momenta, n->N, "momenta");
  debug_scalar_dump(momenta, n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->N_bin,
    data->numerical->dir + string("momenta"));
#endif // DEBUG
}

void draw_real_space_momenta(struct HAMIL_DATA *hd, gsl_rng *seed,
  real_prec *momenta) {
  struct HAMIL_NUMERICAL *n = hd->numerical;

#ifdef MULTITHREAD_RNG
#pragma omp parallel for
#endif  // MULTITHREAD_RNG
  for (unsigned i = 0 ; i < n->N1; i++)
    for (unsigned j = 0 ; j < n->N2; j++)
      for (unsigned k = 0 ; k < n->N3; k++) {
        ULONG iind = k + n->N3 * (j + n->N2 * i);
        real_prec sigma = sqrt(hd->mass_r[iind]);  // real_prec(N));  // care

        momenta[iind] = static_cast<real_prec>(sigma *
          static_cast<real_prec>(GR_NUM(seed, 1., 0)));
      }
}


void draw_masked_momenta(struct HAMIL_DATA *hd, gsl_rng *seed,
  real_prec *momenta) {
  struct HAMIL_NUMERICAL *n = hd->numerical;

#ifdef MULTITHREAD_RNG
#pragma omp parallel for
#endif  // MULTITHREAD_RNG
  for (ULONG i = 0; i < n->N; i++)
    momenta[i] = static_cast<real_prec>(GR_NUM(seed, 1., 0));

  // fftw_array<complex_prec> AUX(n->N), dummyC(n->N);
  fftw_array<complex_prec> momentaC(n->Nhalf);

  if (masking2) {
    fftw_array<real_prec> dummy(n->N);
    // X^{1/2}p
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
    for (ULONG i = 0; i < n->N; i++) {
      real_prec corr = hd->corrf[i];
      real_prec winprime = sqrt(num_1+corr);
      dummy[i] = momenta[i]*winprime;
    }
    // complexify_array(dummy, dummyC, n->N);
    // FFT3d(n->N1, n->N2, n->N3, true, dummyC, AUX);
    fftR2C(n->N1, n->N2, n->N3, dummy, momentaC);
  } else {
    // complexify_array(momenta, dummyC, n->N);
    // FFT3d(n->N1, n->N2, n->N3, true, dummyC, AUX);
    fftR2C(n->N1, n->N2, n->N3, momenta, momentaC);
  }

  unsigned N3half = n->N3/2 + 1;

#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (unsigned i = 0; i < n->N1; i++)
    for (unsigned j = 0; j < n->N2; j++)
      for (unsigned k = 0; k < N3half; k++) {
        ULONG ix = k + n->N3 * (j + n->N2 * i);
        ULONG ix_C = k + N3half * (j + n->N2 * i);

        real_prec mass = 0.0;
#ifdef FOURIER_DEF_1
        mass = hd->mass_f[ix];  // testing
#endif
#ifdef FOURIER_DEF_2
        mass = hd->mass_f[ix] * static_cast<real_prec>(n->N);  // care
        throw runtime_error("Look into normalization of draw_masked_momenta! Seems wrong now.");
#endif

        real_prec sigma = sqrt(mass);

        re(momentaC[ix_C]) *= sigma;
        im(momentaC[ix_C]) *= sigma;
      }

  // FFT3d(n->N1, n->N2, n->N3, false, AUX, dummyC);
  // real_part_array(dummyC, momenta, n->N);
  fftC2R(n->N1, n->N2, n->N3, momentaC, momenta);
}
