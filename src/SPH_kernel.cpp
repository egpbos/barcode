/*
 * Barcode
 * Copyright E.G.P. Bos
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <stdexcept>
#include <cmath> // sqrt, abs
#include <algorithm> // find
#include <iterator>  // distance

#include "SPH_kernel.hpp"

real_prec SPH_kernel_radius(int particle_kernel_type,
                            real_prec particle_kernel_h) {
  real_prec radius;
  switch (particle_kernel_type)
  {
    case 0: // SPH kernel
      radius = particle_kernel_h * 2;
      break;
    default:
      throw std::runtime_error("\"Eat my shorts!\" -- B. J. Simpson");
  }
  return(radius);
}

int SPH_kernel_3D_cells_count(int particle_kernel_type,
                              real_prec particle_kernel_h,
                              real_prec d1, real_prec d2, real_prec d3) {

  // real_prec kernel_reach = SPH_kernel_radius(hd);
  real_prec kernel_reach = SPH_kernel_radius(particle_kernel_type,
                                             particle_kernel_h);

  int reach1 = static_cast<int>(kernel_reach/d1) + 1;
  int reach2 = static_cast<int>(kernel_reach/d2) + 1;
  int reach3 = static_cast<int>(kernel_reach/d3) + 1;
  real_prec kernel_reach_sq = kernel_reach * kernel_reach;

  int out = 0;

  for (int i1 = -reach1; i1 <= reach1; ++i1)
    for (int i2 = -reach2; i2 <= reach2; ++i2)
      for (int i3 = -reach3; i3 <= reach3; ++i3) {
        real_prec dx = (std::abs(static_cast<real_prec>(i1)) - 0.5)*d1;
        real_prec dy = (std::abs(static_cast<real_prec>(i2)) - 0.5)*d2;
        real_prec dz = (std::abs(static_cast<real_prec>(i3)) - 0.5)*d3;
        real_prec r_sq = dx*dx + dy*dy + dz*dz;  // squared for efficiency

        if (r_sq <= kernel_reach_sq) {
          ++out;
        }
      }

  return(out);
}


void SPH_kernel_3D_cells(int particle_kernel_type, real_prec particle_kernel_h,
                         real_prec d1, real_prec d2, real_prec d3,
                         std::vector<int> &out_i,
                         std::vector<int> &out_j, std::vector<int> &out_k) {
  // First determine the reach of the kernel, which determines how many cells
  // every particle needs to loop over to check for contributions. This is
  // based on the kernel radius (e.g. 2*kernel_h for SPH with splines).
  real_prec kernel_reach = SPH_kernel_radius(particle_kernel_type,
                                             particle_kernel_h);
  int reach1 = static_cast<int>(kernel_reach/d1) + 1;
  int reach2 = static_cast<int>(kernel_reach/d2) + 1;
  int reach3 = static_cast<int>(kernel_reach/d3) + 1;
  // Determine sphere to loop over (+ 0.5dx, because we need to accomodate all
  // possible positions in a cell; note that this is implemented as -0.5dx
  // instead of just increasing the kernel reach by 0.5dx, because in principle
  // dx can be different in each direction, so that can only be incorporated if
  // done in each direction independently).
  real_prec kernel_reach_sq = kernel_reach * kernel_reach;
  // std::vector<int> kernel_cells;
  // int ix_cell = 0;

  for(int i1 = -reach1; i1 <= reach1; ++i1)
    for(int i2 = -reach2; i2 <= reach2; ++i2)
      for(int i3 = -reach3; i3 <= reach3; ++i3)
      {
        real_prec dx = (std::abs(static_cast<real_prec>(i1)) - 0.5)*d1; // abs, otherwise the "- 0.5" doesn't work for negative indices
        real_prec dy = (std::abs(static_cast<real_prec>(i2)) - 0.5)*d2;
        real_prec dz = (std::abs(static_cast<real_prec>(i3)) - 0.5)*d3;
        real_prec r_sq = dx*dx + dy*dy + dz*dz; // squared for efficiency

        if (r_sq <= kernel_reach_sq)
        {
          out_i.push_back(i1);
          out_j.push_back(i2);
          out_k.push_back(i3);
          // ++ix_cell;
        }
      }

  // return(kernel_cells);
}



// twee mogelijkheden:
// 1. 3 vectors met i&j, k_begin en k_last (niet k_end, inclusive!)
// 2. 2 2D vectors met k_begin en k_end (vierkant) waar indices i en j zijn
// optie 1:
void SPH_kernel_3D_cells_hull_1(const std::vector<int> &i, const std::vector<int> &j,
                                const std::vector<int> &k, std::vector< std::pair<int, int> > &ij_out,
                                std::vector<int> &k_begin, std::vector<int> &k_last) {
  auto N_cells = i.size();

  std::pair<int, int> ij0(i[0], j[0]);
  ij_out.push_back(ij0);
  k_begin.push_back(k[0]);
  k_last.push_back(k[0]);

  for (unsigned int ix = 1; ix < N_cells; ++ix) {
    // cout << ix << endl;
    // cout << i[ix] << " " << j[ix] << endl;
    std::pair<int, int> ij(i[ix], j[ix]);
    auto ij_duplicate = std::find(ij_out.begin(), ij_out.end(), ij);
    if (ij_duplicate == ij_out.end()) {
      // cout << "nieuw" << endl << endl;
      ij_out.push_back(ij);
      k_begin.push_back(k[ix]);
      k_last.push_back(k[ix]);
    } else {
      // cout << "duplicaat" << endl << endl;
      ULONG ix_duplicate = static_cast<ULONG>(std::distance(ij_out.begin(), ij_duplicate));
      if (k_begin[ix_duplicate] > k[ix])
        k_begin[ix_duplicate] = k[ix];
      if (k_last[ix_duplicate] < k[ix])
        k_last[ix_duplicate] = k[ix];
    }
  }
}

// optie 2:
// void SPH_kernel_3D_cells_hull_2(vector<int> &i, vector<int> &j, vector<int> &k,
//                                 vector<int> &k_begin, vector<int> &k_last) {
// auto N_cells = i.size();
// for (int ix = 0; ix < N_cells; ++ix) {}
// }

void grad_SPH_kernel_3D_h_units(real_prec x_h, real_prec y_h, real_prec z_h,
                                real_prec norm,
                                real_prec &out_x, real_prec &out_y, real_prec &out_z) {
  // N.B.: when using this for density estimation, you need to normalize the
  // result afterwards with V/N! See e.g. getDensity_SPH. Same goes for
  // SPH_kernel_3D.

  // norm must be 1/(PI * h^4)

  // derivative of Monaghan kernel W_4
  // real_prec r_sq = x*x + y*y + z*z;
  real_prec q_sq = x_h*x_h + y_h*y_h + z_h*z_h;
  // __m128 d_h = _mm_setr_ps(static_cast<float>(x_h), static_cast<float>(y_h),
  // static_cast<float>(z_h), 0.0);
  // real_prec q_sq = static_cast<real_prec>(_mm_cvtss_f32(_mm_dp_ps(d_h, d_h, 0x71)));

  real_prec partial;
  if (q_sq > 4)
    partial = 0.;
  else if (q_sq > 1) {
    // real_prec r = sqrt(r_sq);
    // real_prec q = r * h_inv;
    // real_prec qmin2 = q - 2;
    // partial = -0.75*qmin2*qmin2 * norm / r;
    real_prec q = std::sqrt(q_sq);
    real_prec qmin2 = q - 2;
    partial = -0.75*qmin2*qmin2 * norm / q;
  } else {
    // real_prec r = sqrt(r_sq);
    // partial = (2.25*r*h_inv*h_inv - 3*h_inv) * norm;
    real_prec q = std::sqrt(q_sq);
    partial = (2.25*q - 3) * norm;
  }
  // int q_sq_i = static_cast<int>(q_sq);
  // switch (q_sq_i) {
  //   case 0:
  //     {
  //       real_prec q = sqrt(q_sq);
  //       partial = (2.25*q - 3) * norm;
  //     }
  //     break;
  //   case 1:
  //   case 2:
  //   case 3:
  //     {
  //       real_prec q = sqrt(q_sq);
  //       // real_prec q = r * h_inv;
  //       real_prec qmin2 = q - 2;
  //       // partial = -0.75*qmin2*qmin2 * norm / r;
  //       partial = -0.75*qmin2*qmin2 * norm / q;
  //     }
  //     break;
  //   default:
  //     partial = 0.;
  //     break;
  // }

  out_x = partial * x_h;
  out_y = partial * y_h;
  out_z = partial * z_h;
}
