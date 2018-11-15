/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#ifndef BARCODE_TRANSF_H
#define BARCODE_TRANSF_H

#include <string>
#include "define_opt.h"

void transflpt(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec L2, real_prec L3,
               int filtertype, std::string dir);

#endif //BARCODE_TRANSF_H
