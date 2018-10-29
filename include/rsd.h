/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include "define_opt.h"
#include "fftw_array.h"

void calc_pos_rsd(ULONG Npart, real_prec L3, real_prec xobs, real_prec yobs, real_prec zobs, const real_prec *x, const real_prec *y,
                  const real_prec *z, const real_prec *vx, const real_prec *vy, const real_prec *vz, real_prec *xr, real_prec *yr,
                  real_prec *zr, real_prec ascale, real_prec Omega_M, real_prec Omega_L, bool planepar, bool periodic);
