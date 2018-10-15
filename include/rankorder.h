/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include "define_opt.h"
#include "fftw_array.h"

real_prec rankorder_leclercq_ZA(real_prec delta_ZA);
real_prec rankorder_leclercq_2LPT(real_prec delta_2LPT);
void rankorder_leclercq_ZA(real_prec *delta_ZA, real_prec *delta_Nbody, ULONG size);
void rankorder_leclercq_2LPT(real_prec *delta_2LPT, real_prec *delta_Nbody, ULONG size);
