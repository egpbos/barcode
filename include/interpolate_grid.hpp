/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#ifndef BARCODE_INTERPOLATE_GRID_HPP
#define BARCODE_INTERPOLATE_GRID_HPP

#include "define_opt.h"

void interpolate_CIC(ULONG N1, ULONG N2, ULONG N3, real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, const real_prec *xp, const real_prec *yp, const real_prec *zp, const real_prec *field_grid, ULONG N_OBJ, real_prec *interpolation);
real_prec interpolate_CIC(ULONG N1, ULONG N2, ULONG N3, real_prec L1,
                          real_prec L2, real_prec L3, real_prec d1,
                          real_prec d2, real_prec d3, real_prec xp,
                          real_prec yp, real_prec zp, const real_prec *input);

real_prec interpolate_TSC(ULONG N1, ULONG N2, ULONG N3, real_prec d1, real_prec d2, real_prec d3, real_prec xp, real_prec yp,
                          real_prec zp, const real_prec *field);
void interpolate_TSC(ULONG N1, ULONG N2, ULONG N3, real_prec d1, real_prec d2, real_prec d3, const real_prec *xp, const real_prec *yp,
                     const real_prec *zp, const real_prec *field_grid, ULONG N_OBJ, real_prec *interpolation);

void interpolate_TSC_multi(ULONG N1, ULONG N2, ULONG N3, real_prec d1, real_prec d2, real_prec d3, real_prec *xp,
                           real_prec *yp, real_prec *zp, ULONG N_OBJ, real_prec **field_grids,
                           real_prec **interpolations, int N_fields);

void getCICcells(ULONG N1, ULONG N2, ULONG N3, real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec x, real_prec y, real_prec z, ULONG *cell_index1, ULONG *cell_index2);

void getCICweights(real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec x, real_prec y, real_prec z, const ULONG *cell_index1, real_prec *dx, real_prec *tx);

#endif //BARCODE_INTERPOLATE_GRID_HPP
