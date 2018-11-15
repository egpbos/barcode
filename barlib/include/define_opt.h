/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#pragma once
#include <fftw3.h>

// nullptr isn't implemented in GCC before version 4.6. This template is a
// temporary solution.
#ifndef __INTEL_COMPILER  // icc also defines GNUC, but we don't want this there
#ifndef __clang__  // same story for clang
#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 \
    + __GNUC_MINOR__ * 100 \
    + __GNUC_PATCHLEVEL__)
// Test for GCC < 4.6
#if GCC_VERSION < 40600
const                        // this is a const object...
class {
public:
  template<class T>          // convertible to any type
    operator T*() const      // of null non-member
    { return 0; }            // pointer...
  template<class C, class T> // or any type of null
    operator T C::*() const  // member pointer...
    { return 0; }
private:
  void operator&() const;    // whose address can't be taken
} nullptr = {};              // and whose name is nullptr
#endif  // GCC_VERSION
#endif  // __GNUC__
#endif  // __clang__
#endif  // __INTEL_COMPILER

// EGP: some often used variables and functions
#define gsl_real double
#define ULONG unsigned long
#define re(c) ((c)[0])
#define im(c) ((c)[1])

// EGP: behaviour of FFTW routines
#define FFTW_OPTION FFTW_ESTIMATE
#define FORWARD FFTW_FORWARD 
#define BACKWARD FFTW_BACKWARD
#ifdef SINGLE_PREC
#define fftwf_real float
#define real_prec fftwf_real
#define complex_prec fftwf_complex
#endif
#ifdef DOUBLE_PREC
#define fftw_real double
#define real_prec fftw_real
#define complex_prec fftw_complex
#endif

// EGP: more variables (must come after definition of real_prec)
#define num_1 static_cast<real_prec>(1.)
#define num_2 static_cast<real_prec>(2.)
#define num_0_1 static_cast<real_prec>(0.1)
#define num_0_5 static_cast<real_prec>(0.5)

// EGP: no idea what these things do:
#define fwd true
//#define inv false // does nothing in whole code AND causes problems with gsl_linalg (needed in mass_fourier_vs_matrix unittest)!
#define fourier_space true
#define real_space false
#define to_Fspace true
#define to_Rspace false

// EGP: units and physical constants. Not in Makefile, calculations.
#define cgs_Mpc 1. // h^{-1}
#define cgs_sec 1.
#define cgs_km (0.3240779290e-19 * cgs_Mpc)
#define cgs_clight (0.9715611892e-14 * cgs_Mpc/cgs_sec)

// EGP: numerical constants. Not in Makefile, because whatever.
#define eps 1.e-14
#define epsint 1.0e-6 /* numerical accuracy for integrations */
#define NEVAL 1000    /* numerical integration a2com */


