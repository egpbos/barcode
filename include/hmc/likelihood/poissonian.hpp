/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#ifndef BARCODE_POISSONIAN_HPP
#define BARCODE_POISSONIAN_HPP

#include "define_opt.h"

void poissonian_likelihood_partial_f_delta_x_log_like(struct HAMIL_DATA *hamil_data, real_prec *deltaX, real_prec *dummy);
real_prec poissonian_likelihood_log_like(struct HAMIL_DATA *hamil_data, real_prec *delta);
void poissonian_likelihood_grad_f_delta_x_comp(struct HAMIL_DATA *hamil_data, real_prec *deltaX, real_prec *out,
                                               unsigned int component);

#endif //BARCODE_POISSONIAN_HPP
