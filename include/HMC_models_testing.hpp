/*
 * Barcode
 * Copyright E.G.P. Bos
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#ifndef BARCODE_HMC_MODELS_TESTING_HPP
#define BARCODE_HMC_MODELS_TESTING_HPP

#include "define_opt.h"

void likelihood_calc_h(struct HAMIL_DATA *hamil_data, real_prec *deltaX, real_prec *out);
void likelihood_calc_V_SPH_fourier_TSC(struct HAMIL_DATA *hd, real_prec *part_like, real_prec *out_x, real_prec *out_y,
                                       real_prec *out_z);

#endif //BARCODE_HMC_MODELS_TESTING_HPP
