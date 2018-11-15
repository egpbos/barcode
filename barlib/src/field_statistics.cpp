/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include "field_statistics.h"
#include "fftw_array.h"
#include "convenience.h"
#include "fftwrapper.h"
#include <cmath>
//#include "math_funcs.h"
#include "scale_space.hpp"

// Used in barcoderunner.cc and HMC_mass.cc:
// Calculating power spectrum of the density field
void measure_spectrum(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec L2, real_prec L3,
                      real_prec *signal, real_prec *kmode, real_prec *power, ULONG N_bin)
{
  fftw_array<ULONG> nmode(N_bin);

  auto N = (static_cast<ULONG>(N1) * N2) * N3;

  /* Initialize the arrays */
  fillZero(power, N_bin);
  fillZero(kmode, N_bin);
  fillZero(nmode, N_bin);
  fftw_array<complex_prec> Signal(N);
  fillZero(Signal, N);

  FFT3dR2C (N1,N2,N3,signal,Signal);  // keep this one 3D, otherwise factor 2 in loop!

  /* measure the greatest |k| in the box */
  real_prec kmax=sqrt(k_squared(N1/2,N2/2,N3/2,L1,L2,L3,N1,N2,N3));
  /* bin width in k-space */
  real_prec dk=kmax/real_prec(N_bin);

  /* Compute the power-spectrum : P(k) */
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for(unsigned i=0;i<N1;i++)
    for(unsigned j=0;j<N2;j++)
      for(unsigned k=0;k<N3;k++) {
        real_prec ktot=sqrt(k_squared(i,j,k,L1,L2,L3,N1,N2,N3));
        auto nbin=static_cast<ULONG>(ktot/dk);
        // N.B.: we leave the kmax contribution out; it's only one mode, so very
        //       noisy anyway, but also otherwise we have to sacrifice room for
        //       the better modes, because we have to keep the arrays N_bin size
        if (nbin < N_bin) {
          real_prec akl=re(Signal[k+N3*(j+N2*i)]);
          real_prec bkl=im(Signal[k+N3*(j+N2*i)]);
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
          kmode[nbin]+=1*ktot;
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
          power[nbin]+=(akl*akl+bkl*bkl);
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
          nmode[nbin]+=1;
        }
      }

  // Transform from discrete Fourier coefficients to actual power spectrum (see Martel 2005)
  auto NORM=static_cast<real_prec>(1.0);
#ifdef	  FOURIER_DEF_1
  NORM=static_cast<real_prec>(L1*L2*L3);//4./M_PI);
#endif
#ifdef	  FOURIER_DEF_2
  NORM=static_cast<real_prec>(L1*L2*L3/real_prec(N)/real_prec(N));///4./M_PI);
#endif

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (ULONG l=0;l<N_bin;l++)
    if(nmode[l]>0)
    {
      kmode[l]=kmode[l]/static_cast<real_prec>(nmode[l]);
      power[l]=power[l]/static_cast<real_prec>(nmode[l])*NORM;
    }

}
