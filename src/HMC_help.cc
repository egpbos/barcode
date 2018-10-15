/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include "define_opt.h"
#include "struct_hamil.h"

using namespace std;

void convolveInvCorrFuncWithSignal(struct HAMIL_DATA *hd, real_prec *signal,
                                   real_prec *out, real_prec *corrFunc) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  unsigned N3half = (n->N3/2 + 1);
  // fftw_array<complex_prec> signal_C(n->Nhalf);

#ifdef FOURIER_DEF_1
  real_prec normFS = n->vol;
#endif // FOURIER_DEF_1
#ifdef FOURIER_DEF_2
  real_prec normFS = n->vol/static_cast<real_prec>(n->N);
#endif // FOURIER_DEF_2
#ifdef FOURIER_DEF_2_20151021
  real_prec normFS = static_cast<real_prec>(n->N)/n->vol;
#endif  // FOURIER_DEF_2_20151021
    
  // 1) FT[signal] 
  // fftR2C(n->N1, n->N2, n->N3, signal, signal_C);
  fftR2Cplanned(signal, n->R2Cplan->C, n->R2Cplan);
  // fftR2Cplanned(signal, signal_C, n->R2Cplan);

  // 2) 1/C*FT[signal]
#ifdef MULTITHREAD
#pragma omp parallel for 
#endif // MULTITHREAD
  for (unsigned i=0 ; i<n->N1;i++)
    for (unsigned j=0 ; j<n->N2;j++)
      for (unsigned k = 0 ; k < N3half; ++k) {
        ULONG ix = k + n->N3 * (j + n->N2*i);
        ULONG ix_C = k + N3half * (j + n->N2*i);

        // merged invC and normFS into one variable (saves a multiplication)
        real_prec invC_normFS;
        if (corrFunc[ix]>0.0)
          invC_normFS = normFS / corrFunc[ix];
        else
          invC_normFS = 0.;

        // re(signal_C[ix_C]) *= invC_normFS;
        // im(signal_C[ix_C]) *= invC_normFS;
        re(n->R2Cplan->C[ix_C]) *= invC_normFS;
        im(n->R2Cplan->C[ix_C]) *= invC_normFS;
      }

  // 3) IFT[ 1/C*FT[signal] ]
  // fftC2R(n->N1, n->N2, n->N3, signal_C, out);
  fftC2Rplanned(n->R2Cplan->C, out, n->C2Rplan);
  // fftC2Rplanned(signal_C, out, n->C2Rplan);
}
