/*
 * Barcode
 * Copyright E.G.P. Bos
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#ifndef BARCODE_SCALE_SPACE_HPP
#define BARCODE_SCALE_SPACE_HPP

#include "define_opt.h"

real_prec k_squared(unsigned int i, unsigned int j, unsigned int k, real_prec L1, real_prec L2, real_prec L3,
                    unsigned int N1, unsigned int N2,
                    unsigned int N3);

real_prec calc_kx(unsigned int i, real_prec L1, unsigned int N1);

real_prec calc_ky(unsigned int j, real_prec L2, unsigned int N2);

real_prec calc_kz(unsigned int k, real_prec L3, unsigned int N3);

#endif //BARCODE_SCALE_SPACE_HPP
