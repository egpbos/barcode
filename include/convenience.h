/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#pragma once
#include <string>
#include <vector>
#include "./define_opt.h"

void copyArray(const real_prec *from, real_prec *to, ULONG size);
void copyArray(complex_prec *from, complex_prec *to, ULONG size);
void fillZero(real_prec *out, ULONG size);
void fill_one(real_prec *out, ULONG size);
void fill_one(complex_prec *out, ULONG size);
void fillZero(ULONG *out, ULONG size);
void fillZero(complex_prec *out, ULONG size);
void sum_arrays(const real_prec *array1, const real_prec *array2, real_prec *out,
                ULONG size);
void multiplyArrays(const real_prec *array1, const real_prec *array2, real_prec *out,
                    ULONG size);
void multiply_factor_array(real_prec factor, const real_prec *array, real_prec *out,
                           ULONG size);
void subtract_arrays(const real_prec *array1, const real_prec *array2, real_prec *out,
                     ULONG size);
void add_to_array(real_prec addition, real_prec *out, ULONG size);
void add_to_array(const real_prec *addition, real_prec *out, ULONG size);
void add_to_array(const complex_prec *addition, complex_prec *out, ULONG size);
void flip_sign(const complex_prec *in, complex_prec *out, ULONG size);
void flip_sign(const real_prec *in, real_prec *out, ULONG size);
void complexify_array(const real_prec *in, complex_prec *out, ULONG size);
void real_part_array(const complex_prec *in, real_prec *out, ULONG size);
void imaginary_part_array(const complex_prec *in, real_prec *out, ULONG size);

void conjugate_array(const complex_prec *in, complex_prec *out, ULONG size);
void times_i_array(const complex_prec *in, complex_prec *out, ULONG size);
real_prec absolute_squared(const complex_prec in);
void absolute_squared_array(const complex_prec *in, real_prec *out, ULONG size);
void absolute_squared_array(const complex_prec *in, complex_prec *out, ULONG size);

bool contains_nan(const real_prec *in, ULONG size);
bool contains_nan(const complex_prec *in, ULONG size);

std::vector<std::string> &split(const std::string &s, char delim,
                                std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);
