/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include "define_opt.h"
#include "struct_hamil.h"
#include "HMC.h"

#include "fftw_array.h"


void call_hamil(struct DATA *data, gsl_rng * seed) {
  struct NUMERICAL *dn = data->numerical;

  unsigned N1 = dn->N1;
  unsigned N2 = dn->N2;
  unsigned N3 = dn->N3;

  ULONG  N = N1*N2*N3;

  // // prepare SPH kernel stuff
  // int particle_kernel_type = data->numerical->particle_kernel;
  // real_prec particle_kernel_h = data->numerical->particle_kernel_h;
  // real_prec d1 = data->numerical->d1;
  // real_prec d2 = data->numerical->d2;
  // real_prec d3 = data->numerical->d3;
  // int N_cells = SPH_kernel_3D_cells_count(particle_kernel_type,
  //                                         particle_kernel_h, d1, d2, d3);
  // std::vector<int> cells_i, cells_j, cells_k;
  // SPH_kernel_3D_cells(particle_kernel_type, particle_kernel_h, d1, d2, d3,
  //                     cells_i, cells_j, cells_k);

  // set arrays
  fftw_array<real_prec> A(N), B_f(N), B_r(N), C(N), D(N), E(N);

  struct HAMIL_DATA *hamil_data;
  // memory allocation & initialization (== construction)
  hamil_data = new HAMIL_DATA(data, A, B_f, B_r, C, D, E, dn->kernel_cells_i,
                              dn->kernel_cells_j, dn->kernel_cells_k,
                              dn->N_cells);

  if (hamil_data == nullptr) {
    delete hamil_data;
    throw runtime_error("Hamiltonian: error allocating memory....");
  }


  HamiltonianMC(hamil_data, seed, data);

  // hand over inversion success variable
  dn->INV_SUCCESS = hamil_data->numerical->INV_SUCCESS;

  // EGP: hand over rejections
  dn->rejections += hamil_data->numerical->rejections;

  // hand back possibly change eps_fac
  dn->eps_fac = hamil_data->numerical->eps_fac;

  delete hamil_data;
}
