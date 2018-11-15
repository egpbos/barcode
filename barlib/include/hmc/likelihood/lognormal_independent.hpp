/*
 * Barcode
 * Copyright E.G.P. Bos
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#ifndef BARCODE_LOGNORMAL_INDEPENDENT_HPP
#define BARCODE_LOGNORMAL_INDEPENDENT_HPP

#include "define_opt.h"

real_prec lognormal_likelihood_f_delta_x_i_calc(real_prec rho_c, real_prec delta_min, real_prec deltaX_i);
void lognormal_likelihood_partial_f_delta_x_log_like(struct HAMIL_DATA *hamil_data, real_prec *deltaX, real_prec *dummy);
real_prec lognormal_likelihood_log_like(struct HAMIL_DATA *hamil_data, real_prec *delta);
void lognormal_likelihood_grad_f_delta_x_comp(struct HAMIL_DATA *hamil_data, real_prec *deltaX, real_prec *out,
                                              unsigned int component);

#endif //BARCODE_LOGNORMAL_INDEPENDENT_HPP
