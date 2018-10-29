/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#pragma once

void cellbound(unsigned int N1, unsigned int N2, unsigned int N3, real_prec *v1, real_prec *v2, real_prec *v3);

void cellboundcomp(unsigned int N1, unsigned int N2, unsigned int N3, real_prec *vi);

void getDensity_NGP(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec L2, real_prec L3,
                    real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,
                    const real_prec *xp, const real_prec *yp, const real_prec *zp, const real_prec *Par_mass, ULONG N_OBJ, real_prec *delta);

void getDensity_CIC(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3, const real_prec *xp, const real_prec *yp, const real_prec *zp, const real_prec *Par_mass, ULONG N_OBJ, real_prec *delta, bool weightmass);
void getDensity_CIC_old(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec L2, real_prec L3,
                        real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,
                        const real_prec *xp, const real_prec *yp, const real_prec *zp, const real_prec *Par_mass, ULONG N_OBJ, real_prec *delta);

void getDensity_TSC(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec L2, real_prec L3,
                    real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,
                    const real_prec *xp, const real_prec *yp, const real_prec *zp, const real_prec *Par_mass, ULONG N_OBJ, real_prec *delta);

real_prec SPH_kernel_3D(real_prec r, real_prec h);
void getDensity_SPH(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3, const real_prec *xp, const real_prec *yp, const real_prec *zp, const real_prec *Par_mass, ULONG N_OBJ, real_prec *delta, bool weightmass, real_prec kernel_h);
void getDensity_SPH(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3, const real_prec *xp, const real_prec *yp, const real_prec *zp, const real_prec *Par_mass, ULONG N_OBJ, real_prec *delta, bool weightmass, real_prec kernel_h, bool old_cell_index);

void overdens(unsigned int N1, unsigned int N2, unsigned int N3, const real_prec *in, real_prec *out);
