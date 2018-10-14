/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#pragma once
#include <string>
void INIT_PARAMS(struct DATA *data);

int INIT_COSMOLOGY(struct COSMOLOGY *c, std::string codename);

void set_likelihood_functions(struct DATA *data);

int INIT_OBSERVATIONAL(struct DATA *data, real_prec *POWER, real_prec *SIGNAL,
                       real_prec *SIGNALX, real_prec *WINDOW,
                       real_prec *NOISE_SF, real_prec *NOBS, real_prec *CORRF);

void INIT_FFTW(struct DATA *d, real_prec *in_r2c, real_prec *out_c2r,
               complex_prec *in_c2r, complex_prec *out_r2c);
