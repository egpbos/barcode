/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#pragma once
#include "struct_hamil.h"
#include <gsl/gsl_rng.h>
#include "struct_main.h"

void HamiltonianMC(struct HAMIL_DATA *hamil_data, gsl_rng * seed, struct DATA *data);
