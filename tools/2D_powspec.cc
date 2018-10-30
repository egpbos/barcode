/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <iostream>
#include "define_opt.h"
#include "struct_main.h"
#include "fftw_array.h"
#include <cmath>
#include <sstream>
#include "IOfunctionsGen.h"
//#include "math_funcs.h"
#include "scale_space.hpp"

#include "convenience.h"

using namespace std;

// N.B.: planepar = false is not yet properly implemented! Only use planepar = true!
// This means that the redshift distortions must be applied in the z-direction.
void measure_spec2D(unsigned N1, unsigned N2, unsigned N3,
                    real_prec L1, real_prec L2, real_prec L3,
                    real_prec *signal, real_prec *kmode, real_prec *power,
                    ULONG N_bin, bool planepar = true) {
  ULONG N = N1*N2*N3;
  ULONG N_bin_sq = N_bin * N_bin;

  real_prec NORM = L1*L2*L3 / static_cast<real_prec>(4.*M_PI) / static_cast<real_prec>(N*N);
  //
  // TODO: CHECK IF THIS NORMALIZATION IS RIGHT!!
  //

  fftw_array<complex_prec> Signal(N);
  FFT3dR2C(N1, N2, N3, signal, Signal);

  // measure the greatest |k| in the box
  real_prec kmax = sqrt(k_squared(N1/2, N2/2, N3/2, L1, L2, L3, N1, N2, N3));
  // bin width in k-space
  real_prec dk = kmax / real_prec(N_bin-1);  

  //// Compute the Power Spectrum : P(k)

  // Initialize the arrays
  fillZero(power, N_bin_sq);
  fillZero(kmode, N_bin_sq);
  fftw_array<ULONG> nmode(N_bin_sq);
  fillZero(nmode, N_bin_sq);

  for(unsigned i = 0; i < N1; i++)
    for(unsigned j=0;j<N2;j++)
      for(unsigned k=0;k<N3;k++)
      {           
        real_prec ktot = sqrt(k_squared(i, j, k, L1, L2, L3, N1, N2, N3));
        real_prec kx = calc_kx(i, L1, N1);
        real_prec ky = calc_kz(j, L2, N2);
        real_prec kz = calc_kz(k, L3, N3);

        real_prec kperp, kpar;

        if (planepar == true)
        {
          kpar  = sqrt(kz*kz);
          kperp = sqrt(kx*kx+ky*ky);
        }
        else
        {
          throw std::runtime_error("non-plane-parallel option not yet implemented");
          //real_prec pos1[3],pos2[3];

          //work in progress...!!!
          //pos1[0]=kx;
          //pos1[1]=ky;
          //pos1[2]=kz;

          //pos2[0]=kx;
          //pos2[1]=ky;
          //pos2[2]=kz;

          // EGP: these two were already commented out, the others above I commented out
          //kpar =calc_rpar (pos1,pos2);
          //kperp=calc_rperp(pos1,pos2);
        }

        auto nbin_perp = static_cast<ULONG>(kperp/dk);
        auto nbin_par  = static_cast<ULONG>(kpar/dk);

        if (nbin_perp < N_bin && nbin_par < N_bin )
        {
          //ULONG ii = N_bin * nbin_par + nbin_perp; // => par on "x axis", perp on "y axis". See 2D_corr_fct for why this is wrong (sigma-pi plots).
          ULONG ii = nbin_par + N_bin * nbin_perp; // => perp on "x axis", par on "y axis"
          kmode[ii] += 1 * ktot;
          power[ii] += absolute_squared( Signal[k + N3*(j + N2*i)] );
          nmode[ii] += 1;
        }
      }

#pragma omp parallel for
  for (ULONG l = 0; l < N_bin_sq; ++l)
  {
    if (nmode[l] > 0)
    {     
      kmode[l] =        kmode[l] / static_cast<real_prec>(nmode[l]);
      power[l] = NORM * power[l] / static_cast<real_prec>(nmode[l]);
    }
  }
}


void load_arguments(int argc, char *argv[], string &fname_in, unsigned &N1, real_prec &L1, unsigned &N_bin, string &fname_out)
{
  stringstream N1_arg, L1_arg, N_bin_arg;

  if (argc >= 5)
  {
    fname_in = string(argv[1]);
    N1_arg << argv[2];
    N1_arg >> N1;
    L1_arg << argv[3];
    L1_arg >> L1;
    N_bin_arg << argv[4];
    N_bin_arg >> N_bin;

    if (argc >= 6)
      fname_out = string(argv[5]);
    else
      fname_out = fname_in + string("_pow2D");
  }
  else
  {
    cerr << "Need 4 parameters (file in, N1, L1, N_bin)! N.B.: filenames must be given without extension (which must be .dat). Optional: fname_out (default: fname_in+'_pow2D')." << endl;
    exit(1);
  }
}

int main(int argc, char *argv[])
{
  cout << "Note: plane-parallel RSDs only in this version!" << endl;
  bool planepar = true;

  unsigned N1, N_bin;
  real_prec L1;
  string fname_in, fname_out;
  load_arguments(argc, argv, fname_in, N1, L1, N_bin, fname_out);

  ULONG N = N1*N1*N1;
  unsigned N_bin_sq = N_bin * N_bin;

  fftw_array<real_prec> power2D(N_bin_sq), kmode2D(N_bin_sq);
  fftw_array<real_prec> grid(N);

  // get data from input file
  get_scalar(fname_in, grid, N1, N1, N1);

  // do measurement
  measure_spec2D(N1, N1, N1, L1, L1, L1, grid, kmode2D, power2D, N_bin, planepar);

  // output results
  dump_scalar(kmode2D, N_bin_sq, 1, 1, fname_out + string("_k"));
  dump_scalar(power2D, N_bin_sq, 1, 1, fname_out + string("_P"));

  return(0);
}
