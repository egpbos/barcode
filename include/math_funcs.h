/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#pragma once

#include "define_opt.h"
#include <gsl/gsl_randist.h>
#include <string>

void interpolate_CIC(ULONG N1, ULONG N2, ULONG N3, real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec *xp, real_prec *yp, real_prec *zp, real_prec *field_grid, ULONG N_OBJ, real_prec *interpolation);
real_prec interpolate_CIC(ULONG N1, ULONG N2, ULONG N3, real_prec L1,
                          real_prec L2, real_prec L3, real_prec d1,
                          real_prec d2, real_prec d3, real_prec xp,
                          real_prec yp, real_prec zp, real_prec *input);

real_prec
interpolate_TSC(ULONG N1, ULONG N2, ULONG N3, real_prec d1, real_prec d2, real_prec d3, real_prec xp, real_prec yp,
                real_prec zp, real_prec *field);
void interpolate_TSC(ULONG N1, ULONG N2, ULONG N3, real_prec d1, real_prec d2, real_prec d3, real_prec *xp, real_prec *yp,
                     real_prec *zp, real_prec *field_grid, ULONG N_OBJ, real_prec *interpolation);

void interpolate_TSC_multi(ULONG N1, ULONG N2, ULONG N3, real_prec d1, real_prec d2, real_prec d3, real_prec *xp,
                           real_prec *yp, real_prec *zp, ULONG N_OBJ, real_prec **field_grids,
                           real_prec **interpolations, int N_fields);

void getCICcells(ULONG N1, ULONG N2, ULONG N3, real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec x, real_prec y, real_prec z, ULONG *cell_index1, ULONG *cell_index2);

void getCICweights(real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec x, real_prec y, real_prec z, ULONG *cell_index1, real_prec *dx, real_prec *tx);

void pacman_coordinate(real_prec *x, real_prec L);
void pacman_difference(real_prec *d_x, real_prec L);
real_prec pacman_d_x_from_d_x(real_prec d_x, real_prec L);
int pacman_d_ix_from_d_ix(int d_ix, int N);


void complexfactor_mult (ULONG factor, real_prec in_a, complex_prec *in_b, complex_prec *out );

real_prec factorial(int term);

real_prec k_squared(unsigned int i, unsigned int j, unsigned int k, real_prec L1, real_prec L2, real_prec L3,
                    unsigned int N1, unsigned int N2,
                    unsigned int N3);

real_prec calc_kx(unsigned int i, real_prec L1, unsigned int N1);

real_prec calc_ky(unsigned int j, real_prec L2, unsigned int N2);

real_prec calc_kz(unsigned int k, real_prec L3, unsigned int N3);



void convolve(real_prec L1, real_prec L2, real_prec L3, unsigned N1, unsigned N2, unsigned N3, real_prec *in, real_prec *out,
              real_prec smol, bool zeropad, int filtertype);

void kernelcomp(real_prec L1, real_prec L2, real_prec L3, unsigned N1, unsigned N2, unsigned N3, real_prec smol, int filtertype,
                struct DATA *data);

void convcomp(unsigned N1, unsigned N2, unsigned N3, real_prec *in, real_prec *out, real_prec smol, std::string dir);

void convcompb(unsigned N1, unsigned N2, unsigned N3, real_prec *in, real_prec *out);

real_prec min_arr ( ULONG factor, real_prec *in );

real_prec max_arr ( ULONG factor, real_prec *in );

real_prec mean_arr ( ULONG size, real_prec *in );

void gradfft(unsigned N1, unsigned N2, unsigned N3, real_prec L1, real_prec L2, real_prec L3, real_prec *in,
             real_prec *out, unsigned int dim);

void gradfindif(unsigned N1, real_prec L1, real_prec *in, real_prec *out, unsigned int dim);





real_prec power_mean(real_prec x, real_prec y, real_prec p);



real_prec GR_NUM(gsl_rng * SEED, real_prec sigma ,int GR_METHOD);


void create_GARFIELDR2(int N1,int N2,int N3,real_prec *delta,real_prec * Power,gsl_rng * seed);

void create_GARFIELD_old(int N1,int N2,int N3, real_prec L1, real_prec L2, real_prec L3,real_prec *delta,real_prec * Power,gsl_rng * seed);

void create_GARFIELD(unsigned N1,unsigned N2,unsigned N3, real_prec L1, real_prec L2, real_prec L3,real_prec *delta,real_prec * Power,gsl_rng * seed);

void grad_inv_lap_FS(unsigned N1, unsigned N2, unsigned N3, real_prec L1, real_prec L2, real_prec L3, complex_prec *in,
                     complex_prec *out, unsigned int index, bool rfft = false);


real_prec median_arr (ULONG size, real_prec *in);

int odd(int inputval);

real_prec std_arr (ULONG size, real_prec *in);

