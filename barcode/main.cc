/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

/*
BAyesian Reconstruction of COsimc DEnsity fields

purpose: reconstruct initial and final fields compatible with a dark matter distribution (in redshift space)

input  : cosmology, random seed
output : initial and final density and velocity field on a cartesian grid
*/

// c headers
#ifdef WITH_OPENMP
#include <omp.h>
#else
#include <thread>
#endif  // WITH_OPENMP
#include <fftw3.h>
#include <ncurses.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_randist.h>

// c++ headers
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept> // runtime_error
#ifndef RESTART_FILE
#endif  // RESTART_FILE

// barcode headers
#include "struct_main.h"
#include "fftw_array.h"
#include "init_par.h"
#include "IOfunctionsGen.h"
#include "calc_power.h"
#include "barcoderunner.h"

// EGP: NaN detection (part 1)
#ifdef NAN_DETECTION
#ifdef __linux__
#ifdef __GNUC__
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <cfenv>  // FE_INVALID, FE_OVERFLOW
#endif  // __GNUC__
#elif defined(__APPLE__)
#ifdef __llvm__
#include <xmmintrin.h>  // _MM_*
#endif  // __llvm__
#endif  // OSes
#endif  // NAN_DETECTION
// EGP: end NaN detection (part 1)

using namespace std;

int main(int argc, char *argv[]) {
  // EGP: NaN detection (part 2)
#ifdef NAN_DETECTION
#ifdef __linux__
#ifdef __GNUC__
  feenableexcept(FE_INVALID | FE_OVERFLOW);
#endif  // __GNUC__
#elif defined(__APPLE__)
#ifdef __llvm__
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~static_cast<unsigned int>(_MM_MASK_INVALID));
#endif  // __llvm__
#endif  // OSes
#endif  // NAN_DETECTION
  // EGP: end NaN detection (part 2)

#ifdef WITH_OPENMP
  unsigned N_threads = omp_get_max_threads();
#else
  unsigned N_threads = std::thread::hardware_concurrency();
#endif

#ifdef MULTITHREAD
  printf("\nCompiled with MULTITHREAD. Using %d threads.\n*** Densities will "
         "usually be slightly different every multi-core run, due to floating "
         "point addition order! ***\n", N_threads);
#endif  // MULTITHREAD

#ifdef MULTITHREAD_FFTW
#ifdef SINGLE_PREC
  fftwf_init_threads();
  fftwf_plan_with_nthreads(N_threads);
#endif
#ifdef DOUBLE_PREC
  fftw_init_threads();
  fftw_plan_with_nthreads(static_cast<int>(N_threads));
#endif
  printf("Compiled with MULTITHREAD_FFTW support, with %dthreads\n",
         N_threads);
#endif

  // memory allocation before try/catch
  auto data = new DATA;
  gsl_rng *gBaseRand = gsl_rng_alloc(gsl_rng_mt19937);

  try {
    if (data == nullptr) {
      throw runtime_error("main: error allocating memory....");
    }

  /* read parameter file */
    INIT_PARAMS(data);
    struct NUMERICAL *n = data->numerical;
    fftw_array<real_prec> in_r2c(n->N), out_c2r(n->N);
    fftw_array<complex_prec> in_c2r(n->Nhalf), out_r2c(n->Nhalf);
    INIT_FFTW(data, in_r2c, out_c2r, in_c2r, out_r2c);

#ifndef RESTART_FILE
    if (argc > 1) {
      stringstream str;
      unsigned restart_step;
      str << argv[1];
      str >> restart_step;
      data->numerical->start_at = restart_step;
    } else {
      data->numerical->start_at = 0;
    }
#endif  // RESTART_FILE

    // generate random number
    ULONG seed = n->seed;

    // EGP  const gsl_rng_type *T;
    /*
    if (seed==0) {
      srand(time(NULL)); // initialization for rand()
      seed = time(NULL); // changes random chain based on time  
    }
    */
    gsl_rng_set(gBaseRand, seed);
    // end generate random number

    unsigned N1 = n->N1, N2 = n->N2, N3 = n->N3;

    INIT_COSMOLOGY(data->cosmology, n->codename);

    ULONG N = N1*N2*N3;

    fftw_array<real_prec> A(N), B(N), C(N), D(N), E(N), F(N), G(N);

    INIT_OBSERVATIONAL(data, A, B, C, D, E, F, G);

    string fname = n->dir + string("powerero");
    string fnameps = fname+string(".dat");

    ifstream inStream;
    inStream.open(fnameps.data());

//    if (n->readPS == true) {
      cout << endl;
      cout << "reading power-spectrum from table!" << endl;

      readtab(data);

      dump_scalar(data->observational->Power, N1, N2, N3, fname);
//    } else {
//      if ( inStream.is_open() == true ) {
//        cout << endl;
//        cout << "power-spectrum was computed before!" << endl;
//
//        get_scalar(fname, data->observational->Power, N1, N2, N3);
//      } else {
//        cout << endl;
//        cout << "computing power-spectrum from fitting formulae!" << endl;
//
//        initialize_pow_spec(data);
//
//        dump_scalar(data->observational->Power, N1, N2, N3, fname);
//      }
//    }

    // initialize curses
    CURSES_STRUCT curses("BARCODE", "sample cand eps   Neps K_i   pr_i  li_i  "
                         "K_f   pr_f  li_f  dH     P(a)  a a_N_a");
    data->curses = &curses;

    barcoderunner(data, gBaseRand);

    delete n->R2Cplan;
    delete n->C2Rplan;
  }
  catch(std::runtime_error &caught) {
    cout << endl << "ERROR: " << caught.what() << endl;
  }

  // clean up manually defined and c (gsl & fftw) stuff
  gsl_rng_free(gBaseRand);

  data->numerical->performance_log.close();

  // delete data->numerical->R2Cplan;
  // delete data->numerical->C2Rplan;

  string codename = data->numerical->codename;
  cout << "\n >>> " << codename << "  finished.\n";

  delete data;

#ifdef MULTITHREAD_FFTW
#ifdef SINGLE_PREC
  fftwf_cleanup_threads();
#endif
#ifdef DOUBLE_PREC
  fftw_cleanup_threads();
#endif
#endif
}
