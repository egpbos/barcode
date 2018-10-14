/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#pragma once
#include "struct_hamil.h"

void Hamiltonian_mass(struct HAMIL_DATA *hamil_data, real_prec *signal, struct DATA *data);
void dump_mass_spec(struct HAMIL_DATA *hd, struct DATA *data);
