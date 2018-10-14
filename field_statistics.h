/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#ifndef BARCODE_FIELD_STATISTICS_H
#define BARCODE_FIELD_STATISTICS_H

#include "define_opt.h"

void measure_spectrum(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec L2, real_prec L3,
                      real_prec *signal, real_prec *kmode, real_prec *power, ULONG N_bin);

#endif //BARCODE_FIELD_STATISTICS_H
