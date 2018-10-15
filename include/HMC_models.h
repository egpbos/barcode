/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#pragma once
#include "define_opt.h"
#include "struct_hamil.h"
#include <utility>  // std::pair
#include <vector>

// Prior //
void prior_gaussian_grad_log_prior(struct HAMIL_DATA *hd, real_prec *signal, real_prec *out);
real_prec prior_gaussian_log_prior(struct HAMIL_DATA *hd, real_prec *signal);


// Likelihood //
void grf_likelihood_grad_f_delta_x_comp(struct HAMIL_DATA *, real_prec *, real_prec *, unsigned int);
void grf_likelihood_partial_f_delta_x_log_like(struct HAMIL_DATA *, real_prec *, real_prec *);
real_prec grf_likelihood_log_like(struct HAMIL_DATA *hamil_data, real_prec *delta);

void poissonian_likelihood_partial_f_delta_x_log_like(struct HAMIL_DATA *hamil_data, real_prec *deltaX, real_prec *dummy);
real_prec poissonian_likelihood_log_like(struct HAMIL_DATA *hamil_data, real_prec *delta);
void poissonian_likelihood_grad_f_delta_x_comp(struct HAMIL_DATA *hamil_data, real_prec *deltaX, real_prec *out,
                                               unsigned int component);

// old testing stuff (interpolate):
void gaussian_likelihood_partial_f_delta_x_log_like_interpolate(struct HAMIL_DATA *hamil_data, real_prec *deltaX, real_prec *dummy);
void gaussian_likelihood_grad_f_delta_x_comp_interpolate(struct HAMIL_DATA *hamil_data, real_prec *deltaX,
                                                         real_prec *out, unsigned int component);

void gaussian_likelihood_partial_f_delta_x_log_like(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *dummy);
real_prec gaussian_likelihood_log_like(struct HAMIL_DATA *hamil_data, real_prec *delta);
void gaussian_likelihood_grad_f_delta_x_comp(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out,
                                             unsigned int component);

real_prec lognormal_likelihood_f_delta_x_i_calc(real_prec rho_c, real_prec delta_min, real_prec deltaX_i);
void lognormal_likelihood_partial_f_delta_x_log_like(struct HAMIL_DATA *hamil_data, real_prec *deltaX, real_prec *dummy);
real_prec lognormal_likelihood_log_like(struct HAMIL_DATA *hamil_data, real_prec *delta);
void lognormal_likelihood_grad_f_delta_x_comp(struct HAMIL_DATA *hamil_data, real_prec *deltaX, real_prec *out,
                                              unsigned int component);

void likelihood_grad_log_like(struct HAMIL_DATA *hamil_data, real_prec *delta, real_prec *dummy);
void likelihood_grad_log_like_old(struct HAMIL_DATA *hamil_data, real_prec *delta, real_prec *dummy);

void likelihood_calc_h(struct HAMIL_DATA *hamil_data, real_prec *deltaX, real_prec *out);
void likelihood_calc_h_old(struct HAMIL_DATA *hamil_data, real_prec *deltaX, real_prec *out);
void likelihood_calc_h_simple(struct HAMIL_DATA *hamil_data, real_prec *deltaX, real_prec *out);

// special for grf:
void grf_likelihood_grad_log_like(struct HAMIL_DATA *hd, real_prec *delta, real_prec *out);

real_prec gaussian_exponent_convolution(struct HAMIL_DATA *hamil_data, real_prec *signal, real_prec *var);//EGP,real_prec normFS);

// SPH stuff
//real_prec SPH_kernel_scale(struct HAMIL_DATA *hd);
real_prec SPH_kernel_radius(struct HAMIL_DATA *hd);
void grad_SPH_kernel_3D(real_prec x, real_prec y, real_prec z, real_prec h, std::vector<real_prec> &out);
std::vector<int> SPH_kernel_3D_cells(struct HAMIL_DATA *hd);
void likelihood_calc_V_SPH(struct HAMIL_DATA *hd, real_prec *part_like, real_prec *posx, real_prec *posy, real_prec *posz, real_prec *out_x, real_prec *out_y, real_prec *out_z);
void likelihood_calc_h_SPH(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out);

// testing
void SPH_kernel_3D_cells(struct HAMIL_DATA *hd, std::vector<int> &out_i,
                         std::vector<int> &out_j, std::vector<int> &out_k);
int SPH_kernel_3D_cells_count(struct HAMIL_DATA *hd);
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
void grad_SPH_kernel_3D(real_prec x, real_prec y, real_prec z, real_prec h_inv, real_prec h_sq_inv, real_prec norm,
                        real_prec &out_x, real_prec &out_y, real_prec &out_z);


void pad_array_pacman(real_prec *input, unsigned int N1_in, real_prec *out, unsigned int padding);


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
