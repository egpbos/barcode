/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#pragma once

#include "define_opt.h"

// simple math functions
real_prec factorial(int term);
int odd(int inputval);
real_prec power_mean(real_prec x, real_prec y, real_prec p);

// array functions
real_prec min_arr ( ULONG factor, real_prec *in );
real_prec max_arr ( ULONG factor, real_prec *in );
real_prec mean_arr ( ULONG size, const real_prec *in );
real_prec median_arr (ULONG size, real_prec *in);
real_prec std_arr (ULONG size, real_prec *in);
void complexfactor_mult (ULONG factor, real_prec in_a, complex_prec *in_b, complex_prec *out );
