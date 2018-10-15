/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#pragma once
#include "define_opt.h"

real_prec calc_dcom(real_prec scale_factor, real_prec Omega_M, real_prec Omega_L, real_prec hconst);

real_prec D_growth(real_prec scale_factor, real_prec Omega_M, real_prec Omega_L, real_prec hconst);

//real_prec c_pecvel(real_prec a, real_prec Omega_M, real_prec Omega_L, real_prec H0);
real_prec c_pecvel(real_prec a, real_prec Omega_M, real_prec Omega_L, int term);

real_prec fgrow(real_prec a, real_prec Omega_M, real_prec Omega_L, int term);


