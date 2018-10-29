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
#include <string>

#include <fftw3.h>
#include "define_opt.h"
#include <iostream>

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
  plan_pkg(unsigned int _N1, unsigned int _N2, unsigned int _N3, real_prec *in, complex_prec *out) {
    std::cout << "plan_pkg R2C" << std::endl;

    if (_N1 > INT_MAX || _N2 > INT_MAX || _N3 > INT_MAX) {
      std::string msg = "FFT dimensions must not be larger than INT_MAX!";
      throw std::logic_error(msg);
    }

    N1 = _N1, N2 = _N2, N3 = _N3;
    N = (static_cast<ULONG>(N1) * N2) * N3;
    Nhalf = (static_cast<ULONG>(N1) * N2) * (N3/2 + 1);
    R = in;
    C = out;

    #ifdef SINGLE_PREC
    plan = fftwf_plan_dft_r2c_3d(static_cast<int>(N1), static_cast<int>(N2), static_cast<int>(N3), in, out, FFTW_PATIENT);
    #endif
    #ifdef DOUBLE_PREC
    plan = fftw_plan_dft_r2c_3d(static_cast<int>(N1), static_cast<int>(N2), static_cast<int>(N3), in, out, FFTW_PATIENT);
    #endif
  }

  // C2R
  plan_pkg(unsigned int _N1, unsigned int _N2, unsigned int _N3, complex_prec *in, real_prec *out) {
    std::cout << "plan_pkg C2R" << std::endl;

    if (_N1 > INT_MAX || _N2 > INT_MAX || _N3 > INT_MAX) {
      std::string msg = "FFT dimensions must not be larger than INT_MAX!";
      throw std::logic_error(msg);
    }

    N1 = _N1, N2 = _N2, N3 = _N3;
    N = static_cast<ULONG>(N1*N2*N3);
    Nhalf = static_cast<ULONG>(N1 * N2 * (N3/2 + 1));
    C = in;
    R = out;

    #ifdef SINGLE_PREC
    plan = fftwf_plan_dft_c2r_3d(static_cast<int>(N1), static_cast<int>(N2), static_cast<int>(N3), in, out, FFTW_PATIENT);
    #endif
    #ifdef DOUBLE_PREC
    plan = fftw_plan_dft_c2r_3d(static_cast<int>(N1), static_cast<int>(N2), static_cast<int>(N3), in, out, FFTW_PATIENT);
    #endif
  }

  ~plan_pkg() {
    #ifdef SINGLE_PREC
    fftwf_destroy_plan(plan);
    #endif
    #ifdef DOUBLE_PREC
    fftw_destroy_plan(plan);
    #endif
  }
};

void fftC2Rplanned(complex_prec *in, real_prec *out, struct plan_pkg *plan);
void fftR2Cplanned(real_prec *in, complex_prec *out, struct plan_pkg *plan);
