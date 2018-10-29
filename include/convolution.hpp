/*
 * Barcode
 * Copyright E.G.P. Bos
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#ifndef BARCODE_CONVOLUTION_HPP
#define BARCODE_CONVOLUTION_HPP

#include <string>
#include "define_opt.h"

void convolve(real_prec L1, real_prec L2, real_prec L3, unsigned N1, unsigned N2, unsigned N3, real_prec *in, real_prec *out,
              real_prec smol, bool zeropad, int filtertype);

void kernelcomp(real_prec L1, real_prec L2, real_prec L3, unsigned N1, unsigned N2, unsigned N3, real_prec smol, int filtertype,
                struct DATA *data);

void convcomp(unsigned N1, unsigned N2, unsigned N3, real_prec *in, real_prec *out, real_prec smol, const std::string & dir);

void convcompb(unsigned N1, unsigned N2, unsigned N3, real_prec *in, real_prec *out);

#endif //BARCODE_CONVOLUTION_HPP
