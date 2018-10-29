/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <cmath>

#include "struct_main.h"
#include "fftw_array.h"
#include "pacman.hpp" // pacman_center_on_origin
#include "convenience.h"
#include "IOfunctionsGen.h"

using namespace std;


void measure_corr_grid(unsigned N1, unsigned N2, unsigned N3, real_prec L1, /*real_prec L2, real_prec L3,*/ real_prec d1, real_prec d2, real_prec d3, ULONG N_bin, real_prec *signal, real_prec *rmode, ULONG *nmode, real_prec *corr)
{
  ULONG N = N1 * N2 * N3;
  fftw_array<complex_prec> Signal(N);
  FFT3dR2C(N1, N2, N3, signal, Signal);

  // measure the greatest |r| in the box
  real_prec rmax = L1/2 * sqrt(3);
  // For periodic box, there are no distances longer than sqrt(3)*L1/2.

  // bin width in r-space
  real_prec dr = rmax / static_cast<real_prec>(N_bin);

  // Compute the Power Spectrum : P(k) 
  fftw_array<real_prec> dummy(N);
  absolute_squared_array(Signal, dummy, N); // P(k) -> dummy

  fftw_array<complex_prec> AUX(N);
  complexify_array(dummy, AUX, N); // P(k) -> AUX

  FFT3dC2R (N1, N2, N3, AUX, dummy); // 3D correlation function -> dummy

  // Initialize the arrays
  fillZero(rmode, N_bin);
  fillZero(corr, N_bin);
  fillZero(nmode, N_bin);

  // Linearize the 3D corr function: sum over shells
  for(unsigned i=0;i<N1;i++)
    for(unsigned j=0;j<N2;j++)
      for(unsigned k=0;k<N3;k++)
      { 
        real_prec xpos = pacman_center_on_origin(i, N1, d1);
        real_prec ypos = pacman_center_on_origin(j, N2, d2);
        real_prec zpos = pacman_center_on_origin(k, N3, d3);

        real_prec rtot = sqrt(xpos*xpos + ypos*ypos + zpos*zpos);

        auto nbin = static_cast<ULONG>(rtot/dr);

        //if (rtot <= rmax) // this condition is useless, it will never happen with the given rmax
        //{
          rmode[nbin] += rtot;
          corr[nbin] += dummy[k+N3*(j+N2*i)];
          nmode[nbin] += 1;
        //}
        //else
          //cout << "rtot " << rtot << " is groter dan rmax " << rmax << "! Verhip!" << endl;
      }

  // And finally normalize the shells by the volume of the shells
#pragma omp parallel for
  for (ULONG l = 0; l < N_bin; ++l)
  {
    if(nmode[l] > 0)
    {     
      rmode[l] /= static_cast<real_prec>(nmode[l]);
      corr[l] /= static_cast<real_prec>(static_cast<real_prec>(nmode[l]) * static_cast<real_prec>(N));
    }
  }
}

int main(/*int argc, char *argv[]*/)
{
  //int N1 = 64;
  unsigned N1 = 256;
  real_prec L1 = 1000.;
  unsigned N_bin = 200;

  ULONG N = N1*N1*N1;

  real_prec d1 = L1 / static_cast<real_prec>(N1);

  fftw_array<real_prec> grid(N), rmode(N_bin), corr(N_bin);
  fftw_array<ULONG> nmode(N_bin);

  // First the old regular field (or the 64 cubed one)
  //get_scalar(string("/Users/patrick/Downloads/MD_clusters_Mvir_weighted_histogram"), grid, N1, N1, N1);
  get_scalar(string("/Users/patrick/Downloads/MD_clusters_Mvir_weighted_delta"), grid, N1, N1, N1);

  measure_corr_grid(N1, N1, N1, L1, /*L1, L1,*/ d1, d1, d1, N_bin, grid, rmode, nmode, corr);

  //dump_scalar(rmode, static_cast<int>(N_bin), 1, 1, L1, 0., 0., 0, string("corr_fct_data/corr_r_64_NGP"));
  //dump_scalar(corr, static_cast<int>(N_bin), 1, 1, L1, 0., 0., 0, string("corr_fct_data/corr_eta_64_NGP"));
  dump_scalar(rmode, N_bin, 1, 1, string("corr_fct_data/corr_r_NGP"));
  dump_scalar(corr, N_bin, 1, 1, string("corr_fct_data/corr_eta_NGP"));

  // Then the new regular field (with CIC instead of NGP)
  string dir("/Users/patrick/barcode/testcodes/MD_cluster_density_data/");

  get_scalar(dir + string("delta_regCIC"), grid, N1, N1, N1);

  measure_corr_grid(N1, N1, N1, L1,/* L1, L1,*/ d1, d1, d1, N_bin, grid, rmode, nmode, corr);

  dump_scalar(rmode, N_bin, 1, 1, string("corr_fct_data/corr_r_CIC"));
  dump_scalar(corr, N_bin, 1, 1, string("corr_fct_data/corr_eta_CIC"));

  // Then the redshift-distorted field
  get_scalar(dir + string("delta_rsdCIC"), grid, N1, N1, N1);

  measure_corr_grid(N1, N1, N1, L1, /*L1, L1,*/ d1, d1, d1, N_bin, grid, rmode, nmode, corr);

  dump_scalar(rmode, N_bin, 1, 1, string("corr_fct_data/corr_r_rsd_CIC"));
  dump_scalar(corr, N_bin, 1, 1, string("corr_fct_data/corr_eta_rsd_CIC"));

  return(0);
}
