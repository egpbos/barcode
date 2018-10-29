/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include "struct_main.h"
#include "struct_hamil.h"

#include <cmath>
#include <iomanip>
#include <algorithm>  // std::max_element, std::find
#include <iterator>  // std::distance
#include <vector>
#include <utility> // pair

#include <gsl/gsl_integration.h>

#include "fftw_array.h"

#include "interpolate_grid.hpp" // interpolate_CIC/_TSC
#include "scale_space.hpp" // calc_k*
#include "convolution.hpp" // convolve
#include "gradient.hpp" // grad*
#include "pacman.hpp" // pad_array_pacman

#include "Lag2Eul.h"
#include "cosmo.h"

#include "HMC_help.h"

#include "convenience.h"

#include "HMC_models_testing.hpp" // likelihood_calc_h, likelihood_calc_V_SPH_fourier_TSC
#include "SPH_kernel.hpp"

#include "HMC_models.h"

using namespace std;


// Note: in HMC_{momenta,mass}.cc the mass is tuned (if mass_fs)
// to the Gaussian likelihood model!




// new SPH method stuff //


// double inline __attribute__((fastcall)) sqrt14(double n) {
//   __asm__ ("fld qword ptr [esp+4];"
//            "fsqrt;"
//            "ret 8");
//         // _asm fld qword ptr [esp+4]
//         // _asm fsqrt
//         // _asm ret 8
// }


// #include <ia32intrin.h>






// TODO: check changed signedness warnings in below code, but make sure the changes don't make it slower, since this is an arduously optimized piece of code! Original version of the function put in comments below.

// TODO: make sure somehow (using new type?) that ix/y/z are never larger than INT_MAX! Otherwise casts below can go wrong.
// TODO: same goes for index_xy_part below, which should not be larger than LONG_MAX.
/*
 * N2pad and N3pad are ULONG because we need that to promote to ULONG in the calculation of index_xy_part, but in fact
 * they are not expected to be higher than uint_max! I.e. conceptually they are unsigned int, not ulong.
 */
void _likelihood_calc_V_SPH_kernel_loop_h_units(ULONG N2pad, ULONG N3pad, real_prec d1_h, real_prec d2_h, real_prec d3_h,
                                                vector<pair<int, int> > &ij, vector<int> &k_begin, vector<int> &k_last,
                                                unsigned ix, unsigned iy, unsigned iz, real_prec dpcx_h, real_prec dpcy_h,
                                                real_prec dpcz_h, const real_prec *part_like_padded, unsigned padding,
                                                real_prec grad_SPH_kernel_norm, real_prec &_out_x_j, real_prec &_out_y_j,
                                                real_prec &_out_z_j) {
  real_prec out_x_j = 0., out_y_j = 0., out_z_j = 0.; // to avoid having to add to out_#[j]'s inside the loop, which is expensive, because these arrays are not at all contiguous => factor 0.82 timesaving
  unsigned ix_pad = ix + padding;
  unsigned iy_pad = iy + padding;
  unsigned iz_pad = iz + padding;
  // Note: don't use unsigned int as index below! Slower than signed. See:
  // http://stackoverflow.com/a/2044021/1199693.
  // Can't seem to be able to remove warning (unsigned vs signed comparison).
  // Tried converting the ij.size() to long and int, but is again slower...
  for (vector< pair<int, int> >::size_type ij_ix = 0;
       ij_ix < ij.size(); ++ij_ix) {
    int i1 = ij[ij_ix].first;
    int i2 = ij[ij_ix].second;
    auto kx = static_cast<unsigned>(static_cast<int>(ix_pad) + i1);
    auto ky = static_cast<unsigned>(static_cast<int>(iy_pad) + i2);
    // Cell position (relative to the central cell):
    real_prec diff_x_h = dpcx_h - static_cast<real_prec>(i1)*d1_h;
    real_prec diff_y_h = dpcy_h - static_cast<real_prec>(i2)*d2_h;
    // magic
    ULONG index_xy_part = N3pad*(ky + N2pad*kx);  // N2/3pad will promote kx/y to ULONG
    int kz_begin = k_begin[ij_ix];
    int kz_last = k_last[ij_ix];
    auto index_begin = static_cast<ULONG>(kz_begin + static_cast<long>(index_xy_part + iz_pad));
    ULONG index_end = index_begin + static_cast<ULONG>(kz_last - kz_begin);
    real_prec diff_z_h = dpcz_h - static_cast<real_prec>(kz_begin) * d3_h;
    // the actual loop
//    for (int i3 = kz_begin; i3 <= kz_last; ++i3) {
    for (ULONG index = index_begin; index <= index_end; ++index) {
      real_prec common_part = part_like_padded[index];

      real_prec grad_kernel_x, grad_kernel_y, grad_kernel_z;
      grad_SPH_kernel_3D_h_units(diff_x_h, diff_y_h, diff_z_h, grad_SPH_kernel_norm,
                                 grad_kernel_x, grad_kernel_y, grad_kernel_z);

      out_x_j += common_part * grad_kernel_x;
      out_y_j += common_part * grad_kernel_y;
      out_z_j += common_part * grad_kernel_z;

      // for next iteration:
//      ++index;
      diff_z_h -= d3_h;
    }
  }
  _out_x_j = out_x_j;
  _out_y_j = out_y_j;
  _out_z_j = out_z_j;
}


// NOTE:
// Below the original version of the above function. The original version was optimized by hand.
// However, the compiler complained about implicit conversions, which are now removed from the above
// version. TODO: check at some point which version is faster!
//void _likelihood_calc_V_SPH_kernel_loop_h_units(int N2pad, int N3pad, real_prec d1_h, real_prec d2_h, real_prec d3_h,
//                                               vector<pair<int, int> > &ij, vector<int> &k_begin, vector<int> &k_last,
//                                               int ix, int iy, int iz, real_prec dpcx_h, real_prec dpcy_h,
//                                               real_prec dpcz_h, real_prec *part_like_padded, int padding,
//                                               real_prec grad_SPH_kernel_norm, real_prec &_out_x_j, real_prec &_out_y_j,
//                                               real_prec &_out_z_j) {
//  real_prec out_x_j = 0., out_y_j = 0., out_z_j = 0.; // to avoid having to add to out_#[j]'s inside the loop, which is expensive, because these arrays are not at all contiguous => factor 0.82 timesaving
//  int ix_pad = ix + padding;
//  int iy_pad = iy + padding;
//  int iz_pad = iz + padding;
//  // Note: don't use unsigned int as index below! Slower than signed. See:
//  // http://stackoverflow.com/a/2044021/1199693.
//  // Can't seem to be able to remove warning (unsigned vs signed comparison).
//  // Tried converting the ij.size() to long and int, but is again slower...
//  for (vector< pair<int, int> >::size_type ij_ix = 0;
//       ij_ix < ij.size(); ++ij_ix) {
//    int i1 = ij[ij_ix].first;
//    int i2 = ij[ij_ix].second;
//    int kx = ix_pad + i1;
//    int ky = iy_pad + i2;
//    // Cell position (relative to the central cell):
//    real_prec diff_x_h = dpcx_h - static_cast<real_prec>(i1)*d1_h;
//    real_prec diff_y_h = dpcy_h - static_cast<real_prec>(i2)*d2_h;
//    // magic
//    ULONG index_xy_part = N3pad*(ky + N2pad*kx);
//    int kz_begin = k_begin[ij_ix];
//    int kz_last = k_last[ij_ix];
//    ULONG index = kz_begin + index_xy_part + iz_pad;
//    real_prec diff_z_h = dpcz_h - static_cast<real_prec>(kz_begin) * d3_h;
//    // the actual loop
//    for (int i3 = kz_begin; i3 <= kz_last; ++i3) {
//      real_prec common_part = part_like_padded[index];
//
//      real_prec grad_kernel_x, grad_kernel_y, grad_kernel_z;
//      // grad_SPH_kernel_3D(diff_x_h, diff_y_h, diff_z_h, particle_kernel_h,
//      //                    h_sq, h_inv, h_sq_inv, grad_SPH_kernel_norm,
//      //                    grad_kernel_x, grad_kernel_y, grad_kernel_z);
//      grad_SPH_kernel_3D_h_units(diff_x_h, diff_y_h, diff_z_h, grad_SPH_kernel_norm,
//                                 grad_kernel_x, grad_kernel_y, grad_kernel_z);
//
//      out_x_j += common_part * grad_kernel_x;
//      out_y_j += common_part * grad_kernel_y;
//      out_z_j += common_part * grad_kernel_z;
//
//      // for next iteration:
//      ++index;
//      diff_z_h -= d3_h;
//    }
//  }
//  _out_x_j = out_x_j;
//  _out_y_j = out_y_j;
//  _out_z_j = out_z_j;
//}



// TODO: COMPARE THE OUTPUT OF THIS FUNCTION FROM AFTER 6 JUNE 2017 TO VERSION FROM BEFORE!!!
// TODO: COMPARE THE OUTPUT OF THIS FUNCTION FROM AFTER 6 JUNE 2017 TO VERSION FROM BEFORE!!!
// TODO: COMPARE THE OUTPUT OF THIS FUNCTION FROM AFTER 6 JUNE 2017 TO VERSION FROM BEFORE!!!
// TODO: COMPARE THE OUTPUT OF THIS FUNCTION FROM AFTER 6 JUNE 2017 TO VERSION FROM BEFORE!!!
// TODO: COMPARE THE OUTPUT OF THIS FUNCTION FROM AFTER 6 JUNE 2017 TO VERSION FROM BEFORE!!!
// TODO: COMPARE THE OUTPUT OF THIS FUNCTION FROM AFTER 6 JUNE 2017 TO VERSION FROM BEFORE!!!
// TODO: COMPARE THE OUTPUT OF THIS FUNCTION FROM AFTER 6 JUNE 2017 TO VERSION FROM BEFORE!!!
// TODO: also pad_array_pacman

void likelihood_calc_V_SPH(struct HAMIL_DATA *hd, real_prec *part_like, real_prec *posx, real_prec *posy, real_prec *posz, real_prec *out_x, real_prec *out_y, real_prec *out_z)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;

  // int N_cells = SPH_kernel_3D_cells_count(hd);
  // // const int kernel_cells_i[N_cells], kernel_cells_j[N_cells], kernel_cells_k[N_cells];
  // vector<int> kernel_cells_i, kernel_cells_j, kernel_cells_k;
  // SPH_kernel_3D_cells(hd, kernel_cells_i, kernel_cells_j, kernel_cells_k);
  const vector<int> kernel_cells_i = hd->kernel_cells_i;
  const vector<int> kernel_cells_j = hd->kernel_cells_j;
  const vector<int> kernel_cells_k = hd->kernel_cells_k;
  // const int N_cells = hd->N_cells;

  std::vector< std::pair<int, int> > ij;
  std::vector<int> k_begin, k_last;
  SPH_kernel_3D_cells_hull_1(kernel_cells_i, kernel_cells_j, kernel_cells_k,
                             ij, k_begin, k_last);

  unsigned padding = static_cast<unsigned>(*max_element(kernel_cells_i.begin(), kernel_cells_i.end()));
  fftw_array<real_prec> part_like_padded((n->N1+2*padding)*(n->N2+2*padding)*
                                         (n->N3+2*padding));
  pad_array_pacman(part_like, n->N1, part_like_padded, padding);

  // Normalization given that we want the integral over density to be V (left-
  // hand side) and the integral over all particle kernels to be N (right-hand
  // side).
  // Basically, this is the mass of all particles given that we want the
  // density to be rho_c.
  real_prec normalize = hd->rho_c * n->L1*n->L2*n->L3/static_cast<real_prec>(n->N1*n->N2*n->N3);

  real_prec f1;
  bool rsd_model = hd->rsd_model;
  bool planepar = n->planepar;
  if (rsd_model)
    f1 = fgrow(hd->ascale, hd->OM, hd->OL, 1); // first order growth factor (Zel'dovich only!)

  real_prec h = n->particle_kernel_h;
  real_prec h_sq = h*h;
  real_prec h_inv = 1. / h;
  // real_prec h_sq_inv = 1. / h_sq;
  real_prec grad_SPH_kernel_norm = 1. / (M_PI*h_sq*h_sq);
  unsigned N2 = n->N2, N3 = n->N3;
  unsigned N3pad = N3 + 2*padding;
  unsigned N2pad = N2 + 2*padding;
  real_prec d1 = n->d1, d2 = n->d2, d3 = n->d3;
  real_prec d1_h = d1 * h_inv, d2_h = d2 * h_inv, d3_h = d3 * h_inv;

  // Initialize this thing outside the for-loop for performance
  // std::vector<real_prec> grad_kernel(3);
  // Firstprivate below initializes grad_kernel in each thread to a copy
  // of the value outside the loop. Private only initializes an empty
  // vector, which then does not have the correct length.
  #ifdef MULTITHREAD
  // #pragma omp parallel for firstprivate(grad_kernel)
  #pragma omp parallel for
  #endif // MULTITHREAD 
  for (ULONG j = 0; j < n->N; j++)
  {
    // Load particle position
    real_prec px = posx[j], py = posy[j], pz=posz[j];
    // Determine central cell index where particle resides
    // TODO: this calculation will go wrong when posx,y,z are not already put
    // in periodic box coordinates; might be negative then.
    // Addition 6 Jun 2017: TODO: guarantee that positions are in periodic box coords by making special type for that.
    auto ix = static_cast<int>(px/d1);
    auto iy = static_cast<int>(py/d2);
    auto iz = static_cast<int>(pz/d3);
    // Central cell position:
    // real_prec ccx = (static_cast<real_prec>(ix) + 0.5)*d1;
    // real_prec ccy = (static_cast<real_prec>(iy) + 0.5)*d2;
    // real_prec ccz = (static_cast<real_prec>(iz) + 0.5)*d3;
    real_prec ccx_h = (static_cast<real_prec>(ix) + 0.5)*d1_h;
    real_prec ccy_h = (static_cast<real_prec>(iy) + 0.5)*d2_h;
    real_prec ccz_h = (static_cast<real_prec>(iz) + 0.5)*d3_h;
    // Optimization: precalculate diff. pos and cell
    // real_prec dpcx = px - ccx;
    // real_prec dpcy = py - ccy;
    // real_prec dpcz = pz - ccz;
    real_prec dpcx_h = px * h_inv - ccx_h;
    real_prec dpcy_h = py * h_inv - ccy_h;
    real_prec dpcz_h = pz * h_inv - ccz_h;

    real_prec out_x_j, out_y_j, out_z_j;

    // 6 Jun 2017: don't convert ix/y/z to unsigned above, since conversion between uint and double (for ccx_h) is
    //             slower than between int and double!
    _likelihood_calc_V_SPH_kernel_loop_h_units(N2pad, N3pad, d1_h, d2_h, d3_h, ij, k_begin, k_last,
                                               static_cast<unsigned>(ix), static_cast<unsigned>(iy),
                                               static_cast<unsigned>(iz), dpcx_h, dpcy_h, dpcz_h, part_like_padded,
                                               padding, grad_SPH_kernel_norm, out_x_j, out_y_j, out_z_j);

    out_x[j] = normalize * out_x_j;
    out_y[j] = normalize * out_y_j;
    out_z[j] = normalize * out_z_j;

    if (rsd_model) {
      if (!planepar) {
        throw runtime_error("Non-plane-parallel RSD model is not yet implemented in calc_V! Use planepar = true.");
      } else {
        out_z[j] += f1 * out_z[j];
      }
    }
  }
}








void likelihood_calc_h_SPH(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out) {
  // fastest if out == n->C2Rplan->R
  struct HAMIL_NUMERICAL *n = hd->numerical;

  if (!(n->mk == 3)) {
    throw runtime_error("Must use SPH mass kernel (masskernel = 3) when "
                           "using likelihood_calc_h_SPH (calc_h = 2 or 3)!");
  }

  bool rfft = true;

  // fftw_array<real_prec> part_like(n->N), V_x(n->N), V_y(n->N), V_z(n->N);
  // fftw_array<real_prec> V_x(n->N), V_y(n->N), V_z(n->N);
  fftw_array<real_prec> V_y(n->N), V_z(n->N);
  // hd->partial_f_delta_x_log_like(hd, deltaX, part_like);
  hd->partial_f_delta_x_log_like(hd, deltaX, n->R2Cplan->R);

  // The first fft below this switch works fastest if V_x == n->R2Cplan->R.
  // This is allowed in likelihood_calc_V_SPH_fourier_TSC, but it overwrites
  // part_like, because we also use n->R2Cplan->R for that!
  switch (n->calc_h) {
    case 2:
      // likelihood_calc_V_SPH(hd, part_like, hd->posx, hd->posy, hd->posz,
      //                       V_x, V_y, V_z);
      likelihood_calc_V_SPH(hd, n->R2Cplan->R, hd->posx, hd->posy, hd->posz,
                            n->R2Cplan->R, V_y, V_z);
      break;
    case 3:
      // likelihood_calc_V_SPH_fourier_TSC(hd, part_like, hd->posx, hd->posy,
      //                                   hd->posz, V_x, V_y, V_z);
      likelihood_calc_V_SPH_fourier_TSC(hd, n->R2Cplan->R, n->R2Cplan->R, V_y, V_z);
      break;
  }

  // fftw_array<complex_prec> outC(n->Nhalf), dummyC(n->Nhalf);

  // use outC as dummy in first run, saves one addition and one fillZero
  // fftR2C(n->N1, n->N2, n->N3, V_x, outC);
  fftR2Cplanned(n->R2Cplan->R, n->C2Rplan->C, n->R2Cplan);  // use n->C2Rplan->C as outC

  // grad_inv_lap_FS(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, outC, outC, 1, rfft); // multiply by -ik_i/k^2
  grad_inv_lap_FS(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, n->C2Rplan->C, n->C2Rplan->C, 1, rfft); // multiply by -ik_i/k^2

  // fftR2C(n->N1, n->N2, n->N3, V_y, dummyC);
  fftR2Cplanned(V_y, n->R2Cplan->C, n->R2Cplan);  // use n->R2Cplan->C as dummyC
  // grad_inv_lap_FS(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, dummyC, dummyC, 2, rfft); // multiply by -ik_i/k^2
  grad_inv_lap_FS(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, n->R2Cplan->C, n->R2Cplan->C, 2, rfft); // multiply by -ik_i/k^2
  // add_to_array(dummyC, outC, n->Nhalf); // sum the resulting h(k) term to outC
  add_to_array(n->R2Cplan->C, n->C2Rplan->C, n->Nhalf); // sum the resulting h(k) term to outC

  // fftR2C(n->N1, n->N2, n->N3, V_z, dummyC);
  fftR2Cplanned(V_z, n->R2Cplan->C, n->R2Cplan);  // use n->R2Cplan->C as dummyC
  // grad_inv_lap_FS(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, dummyC, dummyC, 3, rfft); // multiply by -ik_i/k^2
  grad_inv_lap_FS(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, n->R2Cplan->C, n->R2Cplan->C, 3, rfft); // multiply by -ik_i/k^2
  // add_to_array(dummyC, outC, n->Nhalf); // sum the resulting h(k) term to outC
  add_to_array(n->R2Cplan->C, n->C2Rplan->C, n->Nhalf); // sum the resulting h(k) term to outC

  // transform h(k) to real space -> h(q)
  // fftC2R(n->N1, n->N2, n->N3, outC, out);
  fftC2Rplanned(n->C2Rplan->C, out, n->C2Rplan);  // use n->R2Cplan->C as dummyC
}
// end SPH method stuff //


// Model: Gaussian likelihood (or any other that is dependent on delta_q only through delta_x) //
void likelihood_grad_log_like(struct HAMIL_DATA *hd, real_prec *delta, real_prec *out)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  // fftw_array<real_prec> dummy(n->N);  // use n->C2Rplan->R as dummy
  
  // single out growing mode (Peebles) -> n->C2Rplan->R (dummy)
  if (n->deltaQ_factor != 1.) {
    multiply_factor_array(n->deltaQ_factor, delta, n->C2Rplan->R, n->N);
  } else {
    copyArray(delta, n->C2Rplan->R, n->N);
  }

  // calculate deltaX and pos
  {
    unsigned facL=1;
    bool reggrid=true;
    gsl_rng *seed = nullptr; // empty: reggrid is true anyway
    real_prec kernel_scale = SPH_kernel_scale(hd);
    if (hd->rsd_model)
      Lag2Eul_rsd_zeldovich(n->C2Rplan->R, hd->deltaX, hd->posx, hd->posy, hd->posz, n->N1, n->N2, n->N3, n->L1, n->L2,
                            n->L3, n->d1, n->d2, n->d3, n->min1, n->min2, n->min3, hd->D1, hd->ascale, hd->OM, hd->OL,
                            n->mk, facL,
                            reggrid, seed, kernel_scale, n->xobs, n->yobs, n->zobs, n->planepar, n->periodic,
                            n->R2Cplan, n->C2Rplan);
    else
      Lag2Eul(n->C2Rplan->R, hd->deltaX, hd->posx, hd->posy, hd->posz, n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1,
              n->d2, n->d3, n->min1, n->min2, n->min3, hd->D1, hd->D2, hd->ascale, hd->OM, hd->OL, hd->sfmodel, n->mk, n->kth,
              facL,
              reggrid, seed, "", kernel_scale, n->R2Cplan, n->C2Rplan);
  }

  // h -> n->C2Rplan->R (dummy)
  switch (n->calc_h)
  {
    case 0:
      likelihood_calc_h(hd, hd->deltaX, n->C2Rplan->R);
      break;
    case 1:
      hd->partial_f_delta_x_log_like(hd, hd->deltaX, n->C2Rplan->R);
      break;
    case 2:
    case 3:
      likelihood_calc_h_SPH(hd, hd->deltaX, n->C2Rplan->R);
      break;
  }

  // now for some normalization terms
  real_prec norm = 1.;

  // ************************* WARNING **************************
  // No longer using heuristic factor, it was caused by a mismatch between mass
  // kernel used for density estimation and the one used for calc_h! As long as
  // the two match, the correspondence between fin. diff. and calculated h is
  // nearly perfect!
  // ************************* WARNING **************************

  // heuristically determined correction factor
  // EGP: I have no idea where this came from, but it seems necessary. See tests
  //      in ipynb "unittest - calc_V_fourier"
  //      NOTE: I'm using n->d1 now; this is of course only correct if all d's
  //            are the same!
  // real_prec heuristic_correction = 1.;
  // switch (n->calc_h)
  // {
  //   case 2:
  //     heuristic_correction = 0.52/0.33;
  //     break;
  //   case 3:
  //     if (n->d1 != n->d2 || n->d2 != n->d3) {
  //       throw runtime_error("likelihood_grad_log_like: d1, d2 and d3 are not"
  //                              " equal, so heuristic_correction is unknown! "
  //                              "Aborting.");
  //     } else {
  //       heuristic_correction = 0.52/0.31/pow(n->d1, 1.04);
  //     }
  //     break;
  // }
  // norm *= heuristic_correction;

  // ************************* WARNING ************************** (see above)

  // EGP: for Zel'dovich, -h is the entire term. For other models more has to be done.
  // FIXME: hier moet nog een factor D1 bij als het signaal alleen delta^(1) is
  real_prec zeldovich_norm = -1.;  // EGP: Zel'dovich: -gradLogLike = -h
  norm *= zeldovich_norm;

  // TODO: check dit theoretisch!
  norm *= n->deltaQ_factor;

  if (n->correct_delta) {
    norm *= hd->D1;
  }

  multiply_factor_array(norm, n->C2Rplan->R, out, n->N);
}
