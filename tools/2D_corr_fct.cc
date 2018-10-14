/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <cmath>
#include <string>
#include <sstream>
#include <iostream>

#include "../struct_main.h"
#include "../fftw_array.h"

#include "../IOfunctionsGen.h"
#include "../convenience.h"

using namespace std;

real_prec pacman_center_on_origin(unsigned ix, unsigned Ni, real_prec di) {
  if (ix <= Ni/2)
    return di * static_cast<real_prec>(ix);
  else
    return -di * static_cast<real_prec>(Ni - ix);
}

void measure_corr2D(unsigned N1, unsigned N2, unsigned N3, real_prec L1,/* real_prec L2,
                    real_prec L3,*/ real_prec d1, real_prec d2, real_prec d3,
                    ULONG N_bin, real_prec *signal, real_prec *rmode,
                    ULONG *nmode, real_prec *corr, bool planepar = true) {
  ULONG N = N1 * N2 * N3;
  ULONG N_bin_sq = N_bin * N_bin;

  cout << "... Fourier transforming signal" << endl;
  fftw_array<complex_prec> Signal(N);
  FFT3dR2C(N1, N2, N3, signal, Signal);

  // measure the greatest |r| in the box
  real_prec rmax = L1/2 * sqrt(3);
  // For periodic box, there are no distances longer than sqrt(3)*L1/2.

  // bin width in r-space
  real_prec dr = rmax / static_cast<real_prec>(N_bin);

  // Compute the Power Spectrum : P(k)
  cout << "... computing power spectrum" << endl;
  absolute_squared_array(Signal, Signal, N);  // P(k) -> Signal (real part)

  cout << "... computing correlation function" << endl;
  fftw_array<real_prec> dummy(N);
  FFT3dC2R(N1, N2, N3, Signal, dummy);  // 3D correlation function -> dummy

  // Initialize the arrays
  fillZero(rmode, N_bin_sq);
  fillZero(corr, N_bin_sq);
  fillZero(nmode, N_bin_sq);

  // Linearize the 3D corr function: sum over shells
  double progress = 0.0;
  cout << "... linearizing the 3D correlation function; progress:" << endl;
  for (unsigned i = 0; i < N1; i++)
    for (unsigned j = 0; j < N2; j++) {
      progress = static_cast<double>(N3*(j+N2*i))/static_cast<double>(N)*100.;
      cout << progress << " %  \r" << flush;
      for (unsigned k = 0; k < N3; k++) {
        real_prec xpos = pacman_center_on_origin(i, N1, d1);
        real_prec ypos = pacman_center_on_origin(j, N2, d2);
        real_prec zpos = pacman_center_on_origin(k, N3, d3);

        real_prec rtot = sqrt(xpos*xpos + ypos*ypos + zpos*zpos);

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

          rmode[ii] += rtot;
          corr[ii] += dummy[k+N3*(j+N2*i)];
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
                                        static_cast<real_prec>(N));
    }
  }
}


// This version is destructive on signal!
void measure_corr2D_destructive(unsigned N1, unsigned N2, unsigned N3, real_prec L1,/* real_prec L2,
                    real_prec L3,*/ real_prec d1, real_prec d2, real_prec d3,
                    ULONG N_bin, real_prec *signal, real_prec *rmode,
                    ULONG *nmode, real_prec *corr, bool planepar = true) {
  cout << "N.B.: this function destroys the contents of signal! Use measure_corr2D for non-destructive." << endl;

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
  double progress = 0.0;
  cout << "... linearizing the 3D correlation function; progress:" << endl;
  for (unsigned i = 0; i < N1; i++)
    for (unsigned j = 0; j < N2; j++) {
      progress = static_cast<double>(N3*(j+N2*i))/static_cast<double>(N)*100.;
      cout << progress << " %  \r" << flush;
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
          rmode[ii] += rtot;
          corr[ii] += signal[k+N3*(j+N2*i)];
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
                                        static_cast<real_prec>(N));
    }
  }
}


void load_arguments(int argc, char *argv[], string &fname_in, unsigned &N1,
                    real_prec &L1, unsigned &N_bin, string &fname_out) {
  stringstream N1_arg, L1_arg, N_bin_arg;

  if (argc >= 5) {
    fname_in = string(argv[1]);
    N1_arg << argv[2];
    N1_arg >> N1;
    L1_arg << argv[3];
    L1_arg >> L1;
    N_bin_arg << argv[4];
    N_bin_arg >> N_bin;

    if (argc >= 6) {
      fname_out = string(argv[5]);
    } else {
      fname_out = fname_in + string("_corr2D");
    }
  } else {
    cerr << "Need 4 parameters (file in, N1, L1, N_bin)!" << endl
         << "N.B.: filenames must be given without extension (which must be "
            ".dat)." << endl
         << "Set N_bin to 0 to automatically get the largest number of bins "
            "without gaps." << endl
         << "Optional: fname_out (default: fname_in+'_corr2D')." << endl;
    exit(1);
  }
}

int main(int argc, char *argv[]) {
  cout << "Note: plane-parallel RSDs only in this version!" << endl;
  bool planepar = true;

  unsigned N1, N_bin;
  real_prec L1;
  string fname_in, fname_out;
  load_arguments(argc, argv, fname_in, N1, L1, N_bin, fname_out);

  if (N_bin > static_cast<unsigned int>(N1))
    cout << "Warning: N_bin larger than N1 will cause gaps in the correlation "
            "function!" << endl;

  real_prec d1 = L1 / static_cast<real_prec>(N1);

  // determine number of bins automatically
  if (N_bin == 0) {
    real_prec rmax = L1/2 * sqrt(3);
    // real_prec dmin = min(d1, min(d2, d3));
    real_prec dmin = d1;
    N_bin = static_cast<unsigned>(ceil(rmax/dmin));
    stringstream N_bin_ss;
    N_bin_ss << N_bin;
    fname_out += string("_Nbin") + N_bin_ss.str();
  }

  unsigned N_bin_sq = N_bin * N_bin;
  ULONG N = N1*N1*N1;

  fftw_array<real_prec> grid(N), rmode(N_bin_sq), corr(N_bin_sq);
  fftw_array<ULONG> nmode(N_bin_sq);

  // get data from input file
  get_scalar(fname_in, grid, N1, N1, N1);

  // do measurement
  measure_corr2D_destructive(N1, N1, N1, L1, /*L1, L1,*/ d1, d1, d1, N_bin, grid, rmode, nmode,
                 corr, planepar);

  // output results
  dump_scalar(rmode, N_bin_sq, 1, 1, fname_out + string("_r"));
  dump_scalar(corr, N_bin_sq, 1, 1, fname_out + string("_eta"));

  return(0);
}
