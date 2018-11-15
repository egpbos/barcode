/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#pragma once

#include "struct_main.h"
#include <gsl/gsl_rng.h>
#include "define_opt.h"

void barcoderunner(struct DATA *data,gsl_rng * gBaseRand);
void load_initial_fields(struct DATA *data, int resnum, gsl_rng *gBaseRand);
unsigned int initial_iteration_number(struct DATA *data);
void setup_random_test(struct DATA *data, real_prec *delta_lag, real_prec *delta_eul, unsigned int facL,
                       real_prec *posx, real_prec *posy, real_prec *posz, gsl_rng *gBaseRand);
void make_initial_guess(struct DATA *data, gsl_rng *gBaseRand);
