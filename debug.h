/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include "define_opt.h"
#include "struct_hamil.h"
#include "struct_main.h"
#include <string>

void debug_array_statistics(real_prec *array, ULONG size, std::string name);
void debug_scalar_dump(real_prec *array, unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec L2,
                       real_prec L3, ULONG N_bin, std::string fname);
void print_data(struct DATA *d);
void print_hamil_data(struct HAMIL_DATA *hd);
