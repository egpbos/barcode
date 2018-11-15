/*
 * Barcode
 * Copyright E.G.P. Bos
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#ifndef BARCODE_GRADIENT_HPP
#define BARCODE_GRADIENT_HPP

#include "define_opt.h"

void gradfft(unsigned N1, unsigned N2, unsigned N3, real_prec L1, real_prec L2, real_prec L3, real_prec *in,
             real_prec *out, unsigned int dim);

void gradfindif(unsigned N1, real_prec L1, const real_prec *in, real_prec *out, unsigned int dim);

void grad_inv_lap_FS(unsigned N1, unsigned N2, unsigned N3, real_prec L1, real_prec L2, real_prec L3, complex_prec *in,
                     complex_prec *out, unsigned int index, bool rfft = false);

#endif //BARCODE_GRADIENT_HPP
