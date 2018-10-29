/*
 * Barcode
 * Copyright E.G.P. Bos
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#ifndef BARCODE_PACMAN_HPP
#define BARCODE_PACMAN_HPP

#include "define_opt.h"

void pacman_coordinate(real_prec *x, real_prec L);
void pacman_difference(real_prec *d_x, real_prec L);
real_prec pacman_d_x_from_d_x(real_prec d_x, real_prec L);
int pacman_d_ix_from_d_ix(int d_ix, int N);
real_prec pacman_center_on_origin(unsigned ix, unsigned Ni, real_prec di);
void pad_array_pacman(real_prec *input, unsigned int N1_in, real_prec *out,
                      unsigned int padding);

#endif //BARCODE_PACMAN_HPP
