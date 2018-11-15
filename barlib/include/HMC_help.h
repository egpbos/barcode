/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#pragma once
#include "define_opt.h"

void convolveInvCorrFuncWithSignal(struct HAMIL_DATA *hamil_data, real_prec *signal, real_prec *out, const real_prec *corrFunc);
