/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#ifndef BARCODE_GAUSSIAN_RANDOM_FIELD_HPP
#define BARCODE_GAUSSIAN_RANDOM_FIELD_HPP

#include "define_opt.h"

void grf_likelihood_grad_f_delta_x_comp(struct HAMIL_DATA *, real_prec *, real_prec *, unsigned int);
void grf_likelihood_partial_f_delta_x_log_like(struct HAMIL_DATA *, real_prec *, real_prec *);
real_prec grf_likelihood_log_like(struct HAMIL_DATA *hamil_data, real_prec *delta);

#endif //BARCODE_GAUSSIAN_RANDOM_FIELD_HPP
