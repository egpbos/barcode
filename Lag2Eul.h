/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#pragma once
void Lag2Eul(real_prec *in, real_prec *dummy, real_prec *posx, real_prec *posy, real_prec *posz, unsigned int N1,
             unsigned int N2, unsigned int N3, real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2,
             real_prec d3, real_prec min1, real_prec min2, real_prec min3, real_prec D1, real_prec D2, real_prec scale,
             real_prec Omega_M, real_prec Omega_L, int sfmodel, int masskernel, real_prec kth, unsigned int facL,
             bool reggrid, gsl_rng *gBaseRand, string dir, real_prec kernel_scale_factor, plan_pkg *R2Cplan,
             plan_pkg *C2Rplan);

void Lag2Eul_zeldovich(real_prec *in, real_prec *dummy, real_prec *posx, real_prec *posy, real_prec *posz,
                       unsigned int N1, unsigned int N2,
                       unsigned int N3, real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2,
                       real_prec d3,
                       real_prec min1, real_prec min2, real_prec min3, real_prec D1, real_prec scale, real_prec Omega_M,
                       real_prec Omega_L, int masskernel, unsigned int facL, bool reggrid, gsl_rng *gBaseRand,
                       real_prec kernel_scale_factor, plan_pkg *R2Cplan, plan_pkg *C2Rplan);

void Lag2Eul_rsd_zeldovich(real_prec *in, real_prec *out, real_prec *posx, real_prec *posy, real_prec *posz,
                           unsigned int N1, unsigned int N2,
                           unsigned int N3, real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2,
                           real_prec d3,
                           real_prec min1, real_prec min2, real_prec min3, real_prec D1, real_prec scale,
                           real_prec Omega_M,
                           real_prec Omega_L, int masskernel, unsigned int facL, bool reggrid, gsl_rng *gBaseRand,
                           real_prec kernel_scale_factor, real_prec xobs, real_prec yobs, real_prec zobs, bool planepar,
                           bool periodic, plan_pkg *R2Cplan, plan_pkg *C2Rplan);
