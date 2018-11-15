/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#pragma once
#include <climits>  // INT_MAX
#include <exception>  // std::exception
#include <stdexcept>  // std::logic_error

#include <fftw3.h>
#include "define_opt.h"

void fftC2R(unsigned int N1, unsigned int N2, unsigned int N3, complex_prec *in, real_prec *out);

void fftR2C(unsigned int N1, unsigned int N2, unsigned int N3, real_prec *in, complex_prec *out);

void FFT3d(unsigned int N1, unsigned int N2, unsigned int N3, bool direction, complex_prec *in, complex_prec *out);

void FFT3dR2C(unsigned int N1, unsigned int N2, unsigned int N3, real_prec *in, complex_prec *out);

void FFT3dC2R(unsigned int N1, unsigned int N2, unsigned int N3, complex_prec *in, real_prec *out);

struct plan_pkg {
  unsigned N1, N2, N3;
  ULONG N;
  ULONG Nhalf;
  complex_prec *C;
  real_prec *R;
  #ifdef SINGLE_PREC
  fftwf_plan plan;
  #endif
  #ifdef DOUBLE_PREC
  fftw_plan plan;
  #endif

  // R2C
  plan_pkg(unsigned int _N1, unsigned int _N2, unsigned int _N3, real_prec *in, complex_prec *out);
  // C2R
  plan_pkg(unsigned int _N1, unsigned int _N2, unsigned int _N3, complex_prec *in, real_prec *out);
  ~plan_pkg();
};

void fftC2Rplanned(complex_prec *in, real_prec *out, struct plan_pkg *plan);
void fftR2Cplanned(real_prec *in, complex_prec *out, struct plan_pkg *plan);
