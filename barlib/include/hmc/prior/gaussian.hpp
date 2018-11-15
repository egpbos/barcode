/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#ifndef BARCODE_GAUSSIAN_HPP
#define BARCODE_GAUSSIAN_HPP

#include "define_opt.h"

void prior_gaussian_grad_log_prior(struct HAMIL_DATA *hd, real_prec *signal, real_prec *out);
real_prec prior_gaussian_log_prior(struct HAMIL_DATA *hd, real_prec *signal);

#endif //BARCODE_GAUSSIAN_HPP
