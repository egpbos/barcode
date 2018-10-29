/*
 * Barcode
 * Copyright E.G.P. Bos
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#ifndef BARCODE_SPH_KERNEL_HPP
#define BARCODE_SPH_KERNEL_HPP

#include <vector>
#include <utility> // pair

#include "define_opt.h"
#include "struct_hamil.h"

real_prec SPH_kernel_radius(int particle_kernel_type,
                            real_prec particle_kernel_h);

int SPH_kernel_3D_cells_count(int particle_kernel_type,
                              real_prec particle_kernel_h,
                              real_prec d1, real_prec d2, real_prec d3);
void SPH_kernel_3D_cells(int particle_kernel_type, real_prec particle_kernel_h,
                         real_prec d1, real_prec d2, real_prec d3,
                         std::vector<int> &out_i,
                         std::vector<int> &out_j, std::vector<int> &out_k);
void SPH_kernel_3D_cells_hull_1(const std::vector<int> &i, const std::vector<int> &j, const std::vector<int> &k,
                                std::vector< std::pair<int, int> > &ij_out,
                                std::vector<int> &k_begin, std::vector<int> &k_last);
void grad_SPH_kernel_3D_h_units(real_prec x_h, real_prec y_h, real_prec z_h,
                                real_prec norm,
                                real_prec &out_x, real_prec &out_y, real_prec &out_z);

template <class T = struct HAMIL_DATA>
real_prec SPH_kernel_scale(T *d)
{
  real_prec scale;
  switch (d->numerical->particle_kernel)
  {
    case 0: // SPH spline kernel
      scale = d->numerical->particle_kernel_h;
      break;
    default:
      scale = 0;
  }
  return(scale);
}

#endif //BARCODE_SPH_KERNEL_HPP
