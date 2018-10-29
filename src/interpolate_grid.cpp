/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <cmath> // abs

#include <gsl/gsl_math.h> // gsl_pow_2

#include "interpolate_grid.hpp"
#include "pacman.hpp"

////////////////////////////////////////////////////////////////
// CIC: helper stuff and interpolation
////////////////////////////////////////////////////////////////


// EGP: convenience function for getting the 8 nearest cells to a coordinate.
// Since these 8 form a cube, only two cells need to be returned, those that
// form a diagonal of the cube. cell_index1 is the one whose cell center
// coordinates xc_1 are lower than the coordinate x's, i.e. xc_1 < x, etc., and
// cell_index2 has the higher cell center coordinates xc_2, xc_2 > x.
// N.B.: cell_index arrays must have length 3!
void getCICcells(ULONG N1, ULONG N2, ULONG N3, real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec x, real_prec y, real_prec z, ULONG *cell_index1, ULONG *cell_index2)
{
  real_prec xpos, ypos, zpos;

  xpos = x - num_0_5*d1;
  ypos = y - num_0_5*d2;
  zpos = z - num_0_5*d3;
  // periodic boundary conditions:
  pacman_coordinate(&xpos, L1);
  pacman_coordinate(&ypos, L2);
  pacman_coordinate(&zpos, L3);

  cell_index1[0] = static_cast<ULONG>(xpos/d1); // indices of the cell of the particle
  cell_index1[1] = static_cast<ULONG>(ypos/d2);
  cell_index1[2] = static_cast<ULONG>(zpos/d3);

  cell_index1[0] = (cell_index1[0] + N1) % N1; // periodic boundary conditions
  cell_index1[1] = (cell_index1[1] + N2) % N2;
  cell_index1[2] = (cell_index1[2] + N3) % N3;

  cell_index2[0] = (cell_index1[0] + 1) % N1;
  cell_index2[1] = (cell_index1[1] + 1) % N2;
  cell_index2[2] = (cell_index1[2] + 1) % N3;
}

// EGP: CIC weights, used in the CIC density estimator and CIC interpolator.
// tx is just 1-dx. Note that these arrays must have length 3, and cell_index1
// as well. The latter can be calculated using getCICcells.
void getCICweights(real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec x, real_prec y, real_prec z,
                   const ULONG *cell_index1, real_prec *dx, real_prec *tx)
{
  real_prec xpos, ypos, zpos, xc, yc, zc;

  xpos = x - num_0_5*d1;
  ypos = y - num_0_5*d2;
  zpos = z - num_0_5*d3;
// periodic boundary conditions:
  pacman_coordinate(&xpos, L1);
  pacman_coordinate(&ypos, L2);
  pacman_coordinate(&zpos, L3);

  xc = static_cast<real_prec>(cell_index1[0]);
  yc = static_cast<real_prec>(cell_index1[1]);
  zc = static_cast<real_prec>(cell_index1[2]);

  dx[0] = xpos/d1 - xc;
  dx[1] = ypos/d2 - yc;
  dx[2] = zpos/d3 - zc;

  tx[0] = num_1 - dx[0];
  tx[1] = num_1 - dx[1];
  tx[2] = num_1 - dx[2];
}

// Single coordinate interpolator without all the array mess.
real_prec interpolate_CIC(ULONG N1, ULONG N2, ULONG N3, real_prec L1,
                          real_prec L2, real_prec L3, real_prec d1,
                          real_prec d2, real_prec d3, real_prec xp,
                          real_prec yp, real_prec zp, const real_prec *input) {
  real_prec dx[3], tx[3];
  ULONG i[3], ii[3];
  getCICcells(N1, N2, N3, L1, L2, L3, d1, d2, d3, xp, yp, zp, i, ii);
  getCICweights(L1, L2, L3, d1, d2, d3, xp, yp, zp, i, dx, tx);

#define FIELD(i, j, k) input[k[2]+N3*(j[1]+N2*i[0])]
  real_prec output = FIELD(i, i, i)    * tx[0]*tx[1]*tx[2] +
                     FIELD(ii, i, i)   * dx[0]*tx[1]*tx[2] +
                     FIELD(i, ii, i)   * tx[0]*dx[1]*tx[2] +
                     FIELD(i, i, ii)   * tx[0]*tx[1]*dx[2] +
                     FIELD(ii, ii, i)  * dx[0]*dx[1]*tx[2] +
                     FIELD(ii, i, ii)  * dx[0]*tx[1]*dx[2] +
                     FIELD(i, ii, ii)  * tx[0]*dx[1]*dx[2] +
                     FIELD(ii, ii, ii) * dx[0]*dx[1]*dx[2];
#undef FIELD

  return output;
}


// EGP: this CIC interpolator was added because getDensity_CIC was used to do
// this, which obviously doesn't work.
void interpolate_CIC(ULONG N1, ULONG N2, ULONG N3, real_prec L1, real_prec L2,
                     real_prec L3, real_prec d1, real_prec d2, real_prec d3,
                     const real_prec *xp, const real_prec *yp, const real_prec *zp,
                     const real_prec *field_grid, ULONG N_OBJ,
                     real_prec *interpolation) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG n = 0; n < N_OBJ; n++) {
    interpolation[n] = interpolate_CIC(N1, N2, N3, L1, L2, L3, d1, d2, d3,
                                       xp[n], yp[n], zp[n], field_grid);
  }
}


////////////////////////////////////////////////////////////////
// TSC: interpolation
////////////////////////////////////////////////////////////////
// Interpolating with TSC is in a sense even easier than CIC,
// since the 27 nearest cells are trivially found; they are the
// cells around the cell that contains the particle.
// Kernel values taken from p. 8 of:
// http://www-hpcc.astro.washington.edu/simulations/DARK_MATTER/adap.ps
////////////////////////////////////////////////////////////////

// Single coordinate interpolator without all the array mess.
real_prec interpolate_TSC(ULONG N1_ul, ULONG N2_ul, ULONG N3_ul, real_prec d1, real_prec d2, real_prec d3, real_prec xp,
                          real_prec yp, real_prec zp, const real_prec *field) {
  auto N1 = static_cast<unsigned>(N1_ul);
  auto N2 = static_cast<unsigned>(N2_ul);
  auto N3 = static_cast<unsigned>(N3_ul);

  real_prec output = 0;
  // particle coordinates in "kernel units" (divided by d1/2/3)
  real_prec xp_kern = xp/d1;
  real_prec yp_kern = yp/d2;
  real_prec zp_kern = zp/d3;
  // indices of the cell of the particle
  auto ix_x_cell = static_cast<unsigned>(xp_kern);
  auto ix_y_cell = static_cast<unsigned>(yp_kern);
  auto ix_z_cell = static_cast<unsigned>(zp_kern);
  // cell_center coordinates (kernel units)
  real_prec x_cell = static_cast<real_prec>(ix_x_cell) + 0.5;
  real_prec y_cell = static_cast<real_prec>(ix_y_cell) + 0.5;
  real_prec z_cell = static_cast<real_prec>(ix_z_cell) + 0.5;
  // distance of particle from cell center
  real_prec dx_p_cell = xp_kern - x_cell;
  real_prec dy_p_cell = yp_kern - y_cell;
  real_prec dz_p_cell = zp_kern - z_cell;

  // weights per dimension
  real_prec wx[3], wy[3], wz[3];
  wx[1] = 0.75 - dx_p_cell*dx_p_cell;
  wy[1] = 0.75 - dy_p_cell*dy_p_cell;
  wz[1] = 0.75 - dz_p_cell*dz_p_cell;
  wx[0] = 0.5 * gsl_pow_2(1.5 - std::abs(dx_p_cell+1));
  wy[0] = 0.5 * gsl_pow_2(1.5 - std::abs(dy_p_cell+1));
  wz[0] = 0.5 * gsl_pow_2(1.5 - std::abs(dz_p_cell+1));
  wx[2] = 0.5 * gsl_pow_2(1.5 - std::abs(dz_p_cell-1));
  wy[2] = 0.5 * gsl_pow_2(1.5 - std::abs(dz_p_cell-1));
  wz[2] = 0.5 * gsl_pow_2(1.5 - std::abs(dz_p_cell-1));

  // indices per dimension, periodic box-ized
  unsigned ix_x[3], ix_y[3], ix_z[3];
  ix_x[1] = ix_x_cell;
  ix_y[1] = ix_y_cell;
  ix_z[1] = ix_z_cell;
  ix_x[0] = (ix_x_cell - 1 + N1) % N1;
  ix_y[0] = (ix_y_cell - 1 + N2) % N2;
  ix_z[0] = (ix_z_cell - 1 + N3) % N3;
  ix_x[2] = (ix_x_cell + 1) % N1;
  ix_y[2] = (ix_y_cell + 1) % N2;
  ix_z[2] = (ix_z_cell + 1) % N3;

  // add the 27 contributions
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k) {
        unsigned ix_f = (ix_x[i]*N2 + ix_y[j])*N3 + ix_z[k];
        output  +=    wx[i]    *   wy[j]     *   wz[k] * field[ix_f] ;
      }

  return output;
}


void interpolate_TSC(ULONG N1, ULONG N2, ULONG N3, real_prec d1, real_prec d2, real_prec d3, const real_prec *xp, const real_prec *yp,
                     const real_prec *zp, const real_prec *field_grid, ULONG N_OBJ, real_prec *interpolation) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG n = 0; n < N_OBJ; n++) {
    interpolation[n] = interpolate_TSC(N1, N2, N3, d1, d2, d3, xp[n], yp[n], zp[n], field_grid);
  }
}


// multi-field variant
// Single coordinate interpolator without all the array mess.
void interpolate_TSC_multi(ULONG N1_ul, ULONG N2_ul, ULONG N3_ul, real_prec d1, real_prec d2, real_prec d3, real_prec xp,
                           real_prec yp, real_prec zp, real_prec **fields, int N_fields, real_prec *out) {
  auto N1 = static_cast<unsigned>(N1_ul);
  auto N2 = static_cast<unsigned>(N2_ul);
  auto N3 = static_cast<unsigned>(N3_ul);

  for (int ix_fields = 0; ix_fields < N_fields; ++ix_fields) {
    out[ix_fields] = 0;
  }
  // particle coordinates in "kernel units" (divided by d1/2/3)
  real_prec xp_kern = xp/d1;
  real_prec yp_kern = yp/d2;
  real_prec zp_kern = zp/d3;
  // indices of the cell of the particle
  auto ix_x_cell = static_cast<unsigned>(xp_kern);
  auto ix_y_cell = static_cast<unsigned>(yp_kern);
  auto ix_z_cell = static_cast<unsigned>(zp_kern);
  // cell_center coordinates (kernel units)
  real_prec x_cell = static_cast<real_prec>(ix_x_cell) + 0.5;
  real_prec y_cell = static_cast<real_prec>(ix_y_cell) + 0.5;
  real_prec z_cell = static_cast<real_prec>(ix_z_cell) + 0.5;
  // distance of particle from cell center
  real_prec dx_p_cell = xp_kern - x_cell;
  real_prec dy_p_cell = yp_kern - y_cell;
  real_prec dz_p_cell = zp_kern - z_cell;

  // weights per dimension
  real_prec wx[3], wy[3], wz[3];
  wx[1] = 0.75 - dx_p_cell*dx_p_cell;
  wy[1] = 0.75 - dy_p_cell*dy_p_cell;
  wz[1] = 0.75 - dz_p_cell*dz_p_cell;
  wx[0] = 0.5 * gsl_pow_2(1.5 - abs(dx_p_cell+1));
  wy[0] = 0.5 * gsl_pow_2(1.5 - abs(dy_p_cell+1));
  wz[0] = 0.5 * gsl_pow_2(1.5 - abs(dz_p_cell+1));
  wx[2] = 0.5 * gsl_pow_2(1.5 - abs(dz_p_cell-1));
  wy[2] = 0.5 * gsl_pow_2(1.5 - abs(dz_p_cell-1));
  wz[2] = 0.5 * gsl_pow_2(1.5 - abs(dz_p_cell-1));

  // indices per dimension, periodic box-ized
  unsigned ix_x[3], ix_y[3], ix_z[3];
  ix_x[1] = ix_x_cell;
  ix_y[1] = ix_y_cell;
  ix_z[1] = ix_z_cell;
  ix_x[0] = (ix_x_cell - 1 + N1) % N1;
  ix_y[0] = (ix_y_cell - 1 + N2) % N2;
  ix_z[0] = (ix_z_cell - 1 + N3) % N3;
  ix_x[2] = (ix_x_cell + 1) % N1;
  ix_y[2] = (ix_y_cell + 1) % N2;
  ix_z[2] = (ix_z_cell + 1) % N3;

  // add the 27 contributions
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k) {
        unsigned ix_f =         (ix_x[i]*N2 + ix_y[j])*N3 + ix_z[k];
        // ...for each field
        for (int ix_fields = 0; ix_fields < N_fields; ++ix_fields) {
          real_prec *field = fields[ix_fields];
          out[ix_fields]  +=    wx[i]    *   wy[j]     *   wz[k] * field[ix_f];
        }
      }
}


void interpolate_TSC_multi(ULONG N1, ULONG N2, ULONG N3, real_prec d1, real_prec d2, real_prec d3, real_prec *xp,
                           real_prec *yp, real_prec *zp, ULONG N_OBJ, real_prec **field_grids,
                           real_prec **interpolations, int N_fields) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG n = 0; n < N_OBJ; n++) {
    real_prec out[N_fields];
    // interpolations[n] = interpolate_TSC(N1, N2, N3, L1, L2, L3, d1, d2, d3,
    //                                    xp[n], yp[n], zp[n], field_grid);
    interpolate_TSC_multi(N1, N2, N3, d1, d2, d3, xp[n], yp[n], zp[n], field_grids, N_fields, out);
    for (int ix_fields = 0; ix_fields < N_fields; ++ix_fields) {
      interpolations[ix_fields][n] = out[ix_fields];
    }
  }
}
