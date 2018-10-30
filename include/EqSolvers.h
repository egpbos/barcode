/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#pragma once

#include "define_opt.h"
// forward declaration
struct plan_pkg;

real_prec linearvel3d(int index, real_prec kx, real_prec ky, real_prec kz, real_prec phi);

void PoissonSolver(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec L2, real_prec L3,
                   real_prec *delta, real_prec *Pot);

void calc_LapPhiv(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, complex_prec *philv, real_prec *LapPhiv,
                  int index1, int index2);

void theta2vel(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec L2, real_prec L3,
               real_prec scale, real_prec Omega_M,
               real_prec Omega_L, real_prec *delta, real_prec *vex, real_prec *vey, real_prec *vez, bool zeropad,
               bool norm, plan_pkg *R2Cplan, plan_pkg *C2Rplan);

void theta2velcomp(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec L2, real_prec L3,
                   real_prec scale, real_prec Omega_M,
                   real_prec Omega_L, real_prec *delta, real_prec *vei, bool zeropad, bool norm, int comp);

void calc_m2v_mem(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec *phiv, real_prec *m2v);


