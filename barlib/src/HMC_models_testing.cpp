/*
 * Barcode
 * Copyright E.G.P. Bos
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <cmath> // sqrt, cos, sin
#include <stdexcept> // runtime_error

#include "struct_hamil.h"
#include "fftw_array.h"
#include "fftwrapper.h" // fftR2Cplanned
#include "scale_space.hpp" // calc_k*
#include "interpolate_grid.hpp" // interpolate_TSC
#include "cosmo.h" // fgrow
#include "convenience.h" // fillZero, add_to_array, multiplyArrays
#include "gradient.hpp" // grad_inv_lap_FS

#include "HMC_models_testing.hpp"


// Model: Gaussian likelihood (or any other that is dependent on delta_q only through delta_x) //
void likelihood_calc_h(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  fftw_array<real_prec> partLike(n->N), dummy(n->N);
  fftw_array<complex_prec> outC(n->Nhalf), dummyC(n->Nhalf);

  fillZero(outC, n->Nhalf);

  hd->partial_f_delta_x_log_like(hd, deltaX, partLike);

  for (unsigned i = 1; i <= 3; i++){
    // compute i-th component of x-gradient of f(deltaX) -> dummy
    hd->grad_f_delta_x_comp(hd, deltaX, dummy, i);
    // multiply with partLike -> dummy
    multiplyArrays(partLike, dummy, dummy, n->N);
    // transform to FS => g_i(k) -> dummyC
    fftR2C(n->N1, n->N2, n->N3, dummy, dummyC);
    // multiply by -ik_i/k^2
    bool rfft = true;
    grad_inv_lap_FS(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, dummyC, dummyC, i, rfft);
    // sum the resulting h(k) term to outC
    add_to_array(dummyC, outC, n->Nhalf);
  }
  // transform h(k) to real space -> h(q)
  fftC2R(n->N1, n->N2, n->N3, outC, out);
}


// new Fourier based version with TSC interpolation
void likelihood_calc_V_SPH_fourier_TSC(struct HAMIL_DATA *hd, real_prec *part_like, real_prec *out_x, real_prec *out_y,
                                       real_prec *out_z) {
  // Works fastest if part_like == n->R2Cplan->R.
  // NOTE: n->R2Cplan->R is used in calc_h_SPH for posx as well, so part_like
  //       will be destroyed if it is also n->R2Cplan->R!

  // 0. initialization
  struct HAMIL_NUMERICAL *n = hd->numerical;
  real_prec h = n->particle_kernel_h;

  // kernel normalization
  real_prec norm_kernel = 24./(h*h*h);

  // Normalization given that we want the integral over density to be V (left-
  // hand side) and the integral over all particle kernels to be N (right-hand
  // side).
  // Basically, this is the mass of all particles given that we want the
  // density to be rho_c.
  real_prec norm_density = hd->rho_c * n->L1*n->L2*n->L3/static_cast<real_prec>(n->N1*n->N2*n->N3);

  real_prec norm = norm_kernel * norm_density;

  // 1. calculate convolution of part_like and kernel
  // 1.a. calculate fourier transform of part_like
  // fftw_array<complex_prec> part_like_F(n->Nhalf), conv_x_F(n->Nhalf),
  //                          conv_y_F(n->Nhalf);
  // fftw_array<complex_prec> conv_x_F(n->Nhalf), conv_y_F(n->Nhalf);
  fftw_array<complex_prec> conv_y_F(n->Nhalf);  // use n->C2Rplan->C as conv_x_F

  // fftR2C(n->N1, n->N2, n->N3, part_like, part_like_F);
  fftR2Cplanned(part_like, n->R2Cplan->C, n->R2Cplan);  // use n->R2Cplan->C as part_like_F

  real_prec L1 = n->L1, L2 = n->L2, L3 = n->L3;
  unsigned N1 = n->N1, N2 = n->N2, N3 = n->N3;
  unsigned N3half = N3/2 + 1;

#ifdef MULTITHREAD
#pragma omp parallel for //schedule(static,chunk) //collapse(3)
#endif // MULTITHREAD
  for (unsigned i = 0; i < N1; ++i) {
    real_prec kx = calc_kx(i, L1, N1);
    for (unsigned j = 0; j < N2; ++j) {
      real_prec ky = calc_ky(j, L2, N2);
      for (unsigned k = 0; k < N3half; ++k) {
        real_prec kz = calc_kz(k, L3, N3);
        real_prec k_sq = kx*kx + ky*ky + kz*kz;
        // 1.b. calculate fourier transform of SPH kernel (it's real, so don't
        //      need complex number)
        real_prec SPH_kernel_F;
        if (k_sq == 0.) {
          SPH_kernel_F = 1./(h*h*h);
        } else {
          real_prec kk = std::sqrt(k_sq);
          real_prec ksink = kk * std::sin(kk);

          SPH_kernel_F = norm * (3 + std::cos(2*kk) - ksink + std::cos(kk) * (ksink - 4))
                         / (k_sq*k_sq*k_sq);
        }

        // 1.c. multiply fourier transforms of part_like and SPH kernel and take
        //      derivative (multiply by ik)
        // Also multiply with h to correct for wrong coordinate in derivative
        // (x(q_i) instead of q=(x(q_i)-x)/h which is the Fourier dual of k)
        ULONG ix = k + N3half*(j + static_cast<ULONG>(N2)*i);
        // x
        // re(conv_x_F[ix]) = h * kx * -im(part_like_F[ix]) * SPH_kernel_F;
        // im(conv_x_F[ix]) = h * kx *  re(part_like_F[ix]) * SPH_kernel_F;
        re(n->C2Rplan->C[ix]) = h * kx * -im(n->R2Cplan->C[ix]) * SPH_kernel_F;
        im(n->C2Rplan->C[ix]) = h * kx *  re(n->R2Cplan->C[ix]) * SPH_kernel_F;
        // y
        // re(conv_y_F[ix]) = h * ky * -im(part_like_F[ix]) * SPH_kernel_F;
        // im(conv_y_F[ix]) = h * ky *  re(part_like_F[ix]) * SPH_kernel_F;
        re(conv_y_F[ix]) = h * ky * -im(n->R2Cplan->C[ix]) * SPH_kernel_F;
        im(conv_y_F[ix]) = h * ky *  re(n->R2Cplan->C[ix]) * SPH_kernel_F;
        // z -> part_like_C itself, to save memory
        // real_prec dummy = re(part_like_F[ix]);
        // re(part_like_F[ix]) = h * kz * -im(part_like_F[ix]) * SPH_kernel_F;
        real_prec dummy = re(n->R2Cplan->C[ix]);
        re(n->R2Cplan->C[ix]) = h * kz * -im(n->R2Cplan->C[ix]) * SPH_kernel_F;
        im(n->R2Cplan->C[ix]) = h * kz * dummy                  * SPH_kernel_F;
      }
    }
  }

  ULONG N_part = n->N;

  // 1.d. transform back to real space => convolution done
  // fftC2R(n->N1, n->N2, n->N3, conv_x_F, out_z);  // out_z is dummy
  fftC2Rplanned(n->C2Rplan->C, n->C2Rplan->R, n->C2Rplan);  // use n->C2Rplan->R as dummy
  // 2. interpolate convolution of part_like and kernel to particle positions
  // interpolate_TSC(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1, n->d2, n->d3,
  //                 hd->posx, hd->posy, hd->posz, out_z, N_part, out_x);
  interpolate_TSC(n->N1, n->N2, n->N3, n->d1, n->d2, n->d3, hd->posx, hd->posy, hd->posz, n->C2Rplan->R, N_part, out_x);

  // fftC2R(n->N1, n->N2, n->N3, conv_y_F, out_z);  // again out_z dummy
  fftC2Rplanned(conv_y_F, n->C2Rplan->R, n->C2Rplan);  // again n->C2Rplan->R as dummy
  // interpolate_TSC(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1, n->d2, n->d3,
  //                 hd->posx, hd->posy, hd->posz, out_z, N_part, out_y);
  interpolate_TSC(n->N1, n->N2, n->N3, n->d1, n->d2, n->d3, hd->posx, hd->posy, hd->posz, n->C2Rplan->R, N_part, out_y);

  // fftw_array<real_prec> conv_x(n->N), conv_y(n->N);  // multi version
  // fftC2R(n->N1, n->N2, n->N3, conv_x_F, conv_x);     // multi version
  // fftC2R(n->N1, n->N2, n->N3, conv_y_F, conv_y);     // multi version
  // here part_like_F is still dummy and we also use part_like as dummy
  // fftC2R(n->N1, n->N2, n->N3, part_like_F, part_like);  // ALSO for multi!
  // fftC2R(n->N1, n->N2, n->N3, n->R2Cplan->C, part_like);  // ALSO for multi!
  fftC2Rplanned(n->R2Cplan->C, n->C2Rplan->R, n->C2Rplan);  // again n->C2Rplan->R as dummy  // ALSO for multi!
  // interpolate_TSC(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1, n->d2, n->d3,
  //                 hd->posx, hd->posy, hd->posz, part_like, N_part, out_z);
  interpolate_TSC(n->N1, n->N2, n->N3, n->d1, n->d2, n->d3, hd->posx, hd->posy, hd->posz, n->C2Rplan->R, N_part, out_z);

  // real_prec *input_fields[3]  = {conv_x, conv_y, part_like}; // multi version
  // real_prec *output_fields[3] = {out_x, out_y, out_z};       // multi version
  // interpolate_TSC_multi(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1,
  //                       n->d2, n->d3, hd->posx, hd->posy, hd->posz, N_part,
  //                       input_fields, output_fields, 3);     // multi version

  // RSD stuff
  if (hd->rsd_model) {
    if (!n->planepar) {
      throw std::runtime_error("Non-plane-parallel RSD model is not yet "
                          "implemented in calc_V! Use planepar = true.");
    } else {
      // first order growth factor (Zel'dovich only!)
      real_prec f1 = fgrow(hd->ascale, hd->OM, hd->OL, 1);
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
      for (unsigned int ix = 0; ix < n->N; ++ix) {
        out_z[ix] += f1 * out_z[ix];
      }
    }
  }

}
