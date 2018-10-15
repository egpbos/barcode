/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#pragma once
void disp_part(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec L2, real_prec L3,
               real_prec d1, real_prec d2, real_prec d3, real_prec *posx, real_prec *posy, real_prec *posz,
               real_prec *psix, real_prec *psiy, real_prec *psiz, real_prec *dummyL, unsigned int facL, bool reggrid,
               gsl_rng *gBaseRand);
