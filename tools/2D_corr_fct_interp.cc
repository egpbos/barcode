/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#if defined(MULTITHREAD) | defined(MULTITHREAD_FFTW) | defined(MULTITHREAD_RNG)
#include <omp.h>
#endif  // openmp headers
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>  // min

#include "../struct_main.h"
#include "../fftw_array.h"

#include "../IOfunctionsGen.h"
#include "../math_funcs.h"  // interpolate_CIC
#include "../convenience.h"

// using namespace std;  // forbidden by Google C++ Style Guide
using std::cout;
using std::cerr;
using std::endl;
using std::string;


void interp_field(real_prec *input, unsigned N1, unsigned N2, unsigned N3, real_prec L1,
                  real_prec L2, real_prec L3, unsigned N1_out, unsigned N2_out,
                  unsigned N3_out, real_prec *output) {
  real_prec dx = L1 / static_cast<real_prec>(N1);
  real_prec dy = L2 / static_cast<real_prec>(N2);
  real_prec dz = L3 / static_cast<real_prec>(N3);
  real_prec dx_out = L1 / static_cast<real_prec>(N1_out);
  real_prec dy_out = L2 / static_cast<real_prec>(N2_out);
  real_prec dz_out = L3 / static_cast<real_prec>(N3_out);

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif  // MULTITHREAD
  for (unsigned i = 0; i < N1_out; ++i)
    for (unsigned j = 0; j < N2_out; ++j)
      for (unsigned k = 0; k < N3_out; ++k) {
        ULONG ix = k + N3_out*(j+N2_out*i);

        // +0.5, assuming we associate grid values with values at center of grid
        // cells
        real_prec posx = dx_out * (0.5 + static_cast<real_prec>(i));
        real_prec posy = dy_out * (0.5 + static_cast<real_prec>(j));
        real_prec posz = dz_out * (0.5 + static_cast<real_prec>(k));

        output[ix] = interpolate_CIC(N1, N2, N3, L1, L2, L3, dx, dy, dz,
                                     posx, posy, posz, input);
      }
}


real_prec pacman_center_on_origin(unsigned ix, unsigned Ni, real_prec di) {
  if (ix <= Ni/2)
    return di * static_cast<real_prec>(ix);
  else
    return -di * static_cast<real_prec>(Ni - ix);
}


// This version is destructive on signal!
void measure_corr2D(unsigned N1, unsigned N2, unsigned N3, real_prec L1, /*real_prec L2,
                    real_prec L3,*/ real_prec d1, real_prec d2, real_prec d3,
                    ULONG N_bin, real_prec *signal, real_prec *rmode,
                    ULONG *nmode, real_prec *corr, real_prec L_max,
                    bool planepar = true) {
  cout << "N.B.: function measure_corr2D in this program destroys the contents of signal! Use measure_corr2D from 2D_corr_fct.cc for non-destructive version." << endl;

  ULONG N = N1 * N2 * N3;
  ULONG Nhalf = (N1/2 + 1) * N2 * N3;
  ULONG N_bin_sq = N_bin * N_bin;

  cout << "... Fourier transforming signal" << endl;
  fftw_array<complex_prec> Signal(Nhalf);
  fftR2C(N1, N2, N3, signal, Signal);

  // Compute the Power Spectrum : P(k)
  cout << "... computing power spectrum" << endl;
  absolute_squared_array(Signal, Signal, Nhalf);  // P(k) -> Signal (real part)

  cout << "... computing correlation function" << endl;
  fftC2R(N1, N2, N3, Signal, signal);  // 3D correlation function -> signal

  // Initialize the arrays
  fillZero(rmode, N_bin_sq);
  fillZero(corr, N_bin_sq);
  fillZero(nmode, N_bin_sq);

  // measure the greatest |r| in the box
  real_prec rmax = L1/2 * sqrt(3);
  // For periodic box, there are no distances longer than sqrt(3)*L1/2.

  // bin width in r-space
  real_prec dr = rmax / static_cast<real_prec>(N_bin);

  // Linearize the 3D corr function: sum over shells
//  double progress = 0.0;
  cout << "... linearizing the 3D correlation function; progress:" << endl;
  #pragma omp parallel for
  for (unsigned i = 0; i < N1; i++)
    for (unsigned j = 0; j < N2; j++) {
      // progress = static_cast<double>(N3*(j+N2*i))/static_cast<double>(N)*100.;
      // cout << progress << " %  \r" << flush;
      for (unsigned k = 0; k < N3; k++) {
        real_prec xpos = pacman_center_on_origin(i, N1, d1);
        real_prec ypos = pacman_center_on_origin(j, N2, d2);
        real_prec zpos = pacman_center_on_origin(k, N3, d3);

        real_prec rpar, rperp;

        if (planepar == true) {
          rpar  = sqrt(zpos*zpos);
          rperp = sqrt(xpos*xpos + ypos*ypos);
        } else {
          // work in progress
          throw std::runtime_error("non-plane-parallel option not yet implemented");
        }

        if (rpar < L_max && rperp < L_max) {
          ULONG nbin_perp = static_cast<ULONG>(rperp/dr);
          ULONG nbin_par  = static_cast<ULONG>(rpar/dr);

          if (nbin_perp < N_bin && nbin_par < N_bin) {
            // par on "x axis", perp on "y axis" (reverse of regular sigma-pi
            // plots in literature!); par is LOS (RSD direction):
            // ULONG ii = N_bin * nbin_par + nbin_perp;
            // perp on "x axis", par on "y axis" (regular order in sigma-pi plots
            // in literature):
            ULONG ii = nbin_par + N_bin * nbin_perp;
            // EGP, 29 July 2016:
            // Right, this was wrong! I was using matplotlib.pcolormesh for
            // plotting the image, but was using x and y the wrong way around
            // (see http://matplotlib.org/examples/pylab_examples/pcolor_demo.html
            // , there you see that pcolorshow has to use y, x = np.mgrid...).
            // Because of this, the plot seemed right, but was actually wrong.
            // With imshow the axes were again flipped, which they shouldn't be.
            // The expected behaviour is that when loading an array with numpy
            // memmap, and reshaping it into a square 2D array, that the first
            // Nx numbers in the (1D) array are plotted along the bottom (with
            // origin='bottom'), from left to right, the next Nx above that, etc.
            // In our case, the horizontal axis is r_perp, so r_perp should vary
            // the fastest. In ULONG ii above this is not the case.
            // So, long story short: this tool does not produce the file I
            // expect.
            // Now that I know, I changed the loading function in Python to flip
            // the axes once again, since I already produced all the 2Dcf files
            // and don't want to do that again :)
            real_prec rtot = sqrt(xpos*xpos + ypos*ypos + zpos*zpos);
#pragma omp atomic
            rmode[ii] += rtot;
#pragma omp atomic
            corr[ii] += signal[k+N3*(j+N2*i)];
#pragma omp atomic
            nmode[ii] += 1;
          }
        }
      }
    }

  // And finally normalize the shells by the volume of the shells
  cout << "... finally, normalizing correlation function shells" << endl;
#pragma omp parallel for
  for (ULONG l = 0; l < N_bin_sq; ++l) {
    if (nmode[l] > 0) {
      rmode[l] /= static_cast<real_prec>(nmode[l]);
      corr[l] /= static_cast<real_prec>(static_cast<real_prec>(nmode[l]) *
                                        static_cast<real_prec>(N));
    }
  }
}


void measure_corr2D_FFTzeropad(real_prec *signal, unsigned N1, unsigned N2, unsigned N3,
                               unsigned N1_interp, unsigned N2_interp, unsigned N3_interp,
                               real_prec L1, real_prec L2, real_prec L3,
                               ULONG N_bin, real_prec *rmode, ULONG *nmode, real_prec *corr,
                               bool planepar = true) {
  ULONG N_interp = N1_interp * N2_interp * N3_interp;
  unsigned N3half = N3/2 + 1;
  ULONG Nhalf = (N1/2 + 1) * N2 * N3;
  unsigned N3half_interp = N3_interp/2 + 1;
  ULONG Nhalf_interp = (N1_interp/2 + 1) * N2_interp * N3_interp;
  ULONG N_bin_sq = N_bin * N_bin;

  real_prec d1_interp = L1 / static_cast<real_prec>(N1_interp);
  real_prec d2_interp = L2 / static_cast<real_prec>(N2_interp);
  real_prec d3_interp = L3 / static_cast<real_prec>(N3_interp);

  cout << "... Fourier transforming signal" << endl;
  fftw_array<complex_prec> Signal(Nhalf);
  fftR2C(N1, N2, N3, signal, Signal);

  // Compute the Power Spectrum : P(k)
  cout << "... computing power spectrum" << endl;
  absolute_squared_array(Signal, Signal, Nhalf);  // P(k) -> Signal (real part)

  // transfer to larger, zeropadded array
  cout << "... copying to larger, zero-padded Fourier space array" << endl;
  fftw_array<complex_prec> Power_Interp(Nhalf_interp);
  fillZero(Power_Interp, Nhalf_interp);
  for (unsigned i = 0; i < N1; i++) {
    unsigned I = i;
    if (i >= N1/2) {
      I = N1_interp - (N1 - i);
    }
    for (unsigned j = 0; j < N2; j++) {
      unsigned J = j;
      if (j >= N2/2) {
        J = N2_interp - (N2 - j);
      }
      for (unsigned k = 0; k < N3half; ++k) {
        ULONG ix_input = k + N3half * (j + N2 * i);
        ULONG ix_interp = k + N3half_interp * (J + N2_interp * I);
        re(Power_Interp[ix_interp]) = re(Signal[ix_input]);
        im(Power_Interp[ix_interp]) = im(Signal[ix_input]);
      }
    }
  }

  cout << "... computing correlation function" << endl;
  fftw_array<real_prec> power_interp(N_interp);
  fftC2R(N1_interp, N3_interp, N3_interp, Power_Interp, power_interp);  // 3D correlation function

  // Initialize the arrays
  fillZero(rmode, N_bin_sq);
  fillZero(corr, N_bin_sq);
  fillZero(nmode, N_bin_sq);

  // measure the greatest |r| in the box
  real_prec rmax = L1/2 * sqrt(3);
  // For periodic box, there are no distances longer than sqrt(3)*L1/2.

  // bin width in r-space
  real_prec dr = rmax / static_cast<real_prec>(N_bin);

  // Linearize the 3D corr function: sum over shells
//  double progress = 0.0;
  cout << "... linearizing the 3D correlation function; progress:" << endl;
#pragma omp parallel for
  for (unsigned i = 0; i < N1_interp; i++)
    for (unsigned j = 0; j < N2_interp; j++) {
      // progress = static_cast<double>(N3*(j+N2*i))/static_cast<double>(N)*100.;
      // cout << progress << " %  \r" << flush;
      for (unsigned k = 0; k < N3_interp; k++) {
        real_prec xpos = pacman_center_on_origin(i, N1_interp, d1_interp);
        real_prec ypos = pacman_center_on_origin(j, N2_interp, d2_interp);
        real_prec zpos = pacman_center_on_origin(k, N3_interp, d3_interp);

        real_prec rpar, rperp;

        if (planepar == true) {
          rpar  = sqrt(zpos*zpos);
          rperp = sqrt(xpos*xpos + ypos*ypos);
        } else {
          // work in progress
          throw std::runtime_error("non-plane-parallel option not yet implemented");
        }

        ULONG nbin_perp = static_cast<ULONG>(rperp/dr);
        ULONG nbin_par  = static_cast<ULONG>(rpar/dr);

        if (nbin_perp < N_bin && nbin_par < N_bin) {
          // par on "x axis", perp on "y axis" (reverse of regular sigma-pi
          // plots in literature!); par is LOS (RSD direction):
          // ULONG ii = N_bin * nbin_par + nbin_perp;
          // perp on "x axis", par on "y axis" (regular order in sigma-pi plots
          // in literature):
          ULONG ii = nbin_par + N_bin * nbin_perp;
          // EGP, 29 July 2016:
          // Right, this was wrong! I was using matplotlib.pcolormesh for
          // plotting the image, but was using x and y the wrong way around
          // (see http://matplotlib.org/examples/pylab_examples/pcolor_demo.html
          // , there you see that pcolorshow has to use y, x = np.mgrid...).
          // Because of this, the plot seemed right, but was actually wrong.
          // With imshow the axes were again flipped, which they shouldn't be.
          // The expected behaviour is that when loading an array with numpy
          // memmap, and reshaping it into a square 2D array, that the first
          // Nx numbers in the (1D) array are plotted along the bottom (with
          // origin='bottom'), from left to right, the next Nx above that, etc.
          // In our case, the horizontal axis is r_perp, so r_perp should vary
          // the fastest. In ULONG ii above this is not the case.
          // So, long story short: this tool does not produce the file I
          // expect.
          // Now that I know, I changed the loading function in Python to flip
          // the axes once again, since I already produced all the 2Dcf files
          // and don't want to do that again :)
          real_prec rtot = sqrt(xpos*xpos + ypos*ypos + zpos*zpos);
#pragma omp atomic
          rmode[ii] += rtot;
#pragma omp atomic
          corr[ii] += power_interp[k+N3_interp*(j+N2_interp*i)];
#pragma omp atomic
          nmode[ii] += 1;
        }
      }
    }

  // And finally normalize the shells by the volume of the shells
  cout << "... finally, normalizing correlation function shells" << endl;
#pragma omp parallel for
  for (ULONG l = 0; l < N_bin_sq; ++l) {
    if (nmode[l] > 0) {
      rmode[l] /= static_cast<real_prec>(nmode[l]);
      corr[l] /= static_cast<real_prec>(static_cast<real_prec>(nmode[l]) *
                                        static_cast<real_prec>(N_interp));
    }
  }
}


void load_arguments(int argc, char *argv[], string &fname_in, unsigned &N1,
                    real_prec &L1, unsigned &N1_interp, unsigned &N_bin,
                    unsigned &interp_mode, real_prec &L_max, string &fname_out) {
  stringstream N1_arg, L1_arg, N1_interp_arg, N_bin_arg, interp_mode_arg, L_max_arg;

  if (argc >= 6) {
    fname_in = string(argv[1]);
    N1_arg << argv[2];
    N1_arg >> N1;
    L1_arg << argv[3];
    L1_arg >> L1;
    N1_interp_arg << argv[4];
    N1_interp_arg >> N1_interp;
    N_bin_arg << argv[5];
    N_bin_arg >> N_bin;
    interp_mode_arg << argv[6];
    interp_mode_arg >> interp_mode;
    L_max_arg << argv[7];
    L_max_arg >> L_max;

    if (argc >= 9) {
      fname_out = string(argv[8]);
    } else {
      fname_out = fname_in + string("_interpCIC") + N1_interp_arg.str() + string("_corr2D");
    }
  } else {
    std::cerr  << "Need 7 parameters (file in, N1, L1, N1_interp, N_bin, interp_mode, L_max)!\n"
               << "N.B.: filenames must be given without extension (which must be .dat).\n"
               << "- Set N_bin to 0 to automatically get the largest number of bins without gaps.\n"
               << "- Interp_mode can be 0 (CIC) or 1 (Fourier zero padding).\n"
               << "- Optional: fname_out (default: fname_in+'_interpCIC[N1_interp]_corr2D')."
               << std::endl;
    exit(1);
  }
}


int main(int argc, char *argv[]) {
#ifdef MULTITHREAD_FFTW
  fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
  printf("Compiled with MULTITHREAD_FFTW support, with %dthreads\n",
         omp_get_max_threads());
#endif

  cout << "Note: plane-parallel RSDs only in this version!" << endl;
  bool planepar = true;

  unsigned N1, N1_interp, N_bin, interp_mode;
  real_prec L1, L_max;
  string fname_in, fname_out;
  load_arguments(argc, argv, fname_in, N1, L1, N1_interp, N_bin, interp_mode, L_max, fname_out);

  ULONG N = N1*N1*N1;
  ULONG N_interp = N1_interp*N1_interp*N1_interp;

  fftw_array<real_prec> input(N);

  // get data from input file
  get_scalar(fname_in, input, N1, N1, N1);

  // corr2D preparations
  if (N_bin > N1_interp)
    cout << "Warning: N_bin larger than N1_interp will cause gaps in the correlation "
            "function!" << endl;

  real_prec d1_interp = L1 / static_cast<real_prec>(N1_interp);

  // determine number of bins automatically
  if (N_bin == 0) {
    real_prec rmax = L1/2 * sqrt(3);
    // real_prec dmin = min(d1, min(d2, d3));
    real_prec dmin = d1_interp;
    N_bin = static_cast<unsigned>(ceil(rmax/dmin));
    stringstream N_bin_ss;
    N_bin_ss << N_bin;
    fname_out += string("_Nbin") + N_bin_ss.str();
  }

  unsigned N_bin_sq = N_bin * N_bin;

  fftw_array<real_prec> rmode(N_bin_sq), corr(N_bin_sq);
  fftw_array<ULONG> nmode(N_bin_sq);

  // interpolate
  switch (interp_mode) {
    case 0: {
      cout << "... interpolate" << endl;
      fftw_array<real_prec> grid(N_interp);
      interp_field(input, N1, N1, N1, L1, L1, L1, N1_interp, N1_interp, N1_interp, grid);
      cout << "... do measurement" << endl;
      measure_corr2D(N1_interp, N1_interp, N1_interp, L1, /*L1, L1,*/ d1_interp, d1_interp, d1_interp, N_bin, grid, rmode, nmode,
                     corr, L_max, planepar);
      break;
    }

    case 1: {
      measure_corr2D_FFTzeropad(input, N1, N1, N1,
          N1_interp, N1_interp, N1_interp,
          L1, L1, L1,
          N_bin, rmode, nmode, corr, planepar);
      break;
    }
  }

  // output results
  dump_scalar(rmode, N_bin_sq, 1, 1, fname_out + string("_r"));
  dump_scalar(corr, N_bin_sq, 1, 1, fname_out + string("_eta"));

#ifdef MULTITHREAD_FFTW
  fftw_cleanup_threads();
#endif
}
