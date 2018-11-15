/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#pragma once
#include "define_opt.h"

void likelihood_grad_log_like(struct HAMIL_DATA *hamil_data, real_prec *delta, real_prec *dummy);
void likelihood_calc_V_SPH(struct HAMIL_DATA *hd, real_prec *part_like, real_prec *posx, real_prec *posy, real_prec *posz, real_prec *out_x, real_prec *out_y, real_prec *out_z);
void likelihood_calc_h_SPH(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out);
