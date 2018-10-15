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

void copyArray(real_prec *from, real_prec *to, ULONG size);
void copyArray(complex_prec *from, complex_prec *to, ULONG size);
void fillZero(real_prec *out, ULONG size);
void fill_one(real_prec *out, ULONG size);
void fill_one(complex_prec *out, ULONG size);
void fillZero(ULONG *out, ULONG size);
void fillZero(complex_prec *out, ULONG size);
void sum_arrays(real_prec *array1, real_prec *array2, real_prec *out,
                ULONG size);
void multiplyArrays(real_prec *array1, real_prec *array2, real_prec *out,
                    ULONG size);
void multiply_factor_array(real_prec factor, real_prec *array, real_prec *out,
                           ULONG size);
void subtract_arrays(real_prec *array1, real_prec *array2, real_prec *out,
                     ULONG size);
void add_to_array(real_prec addition, real_prec *out, ULONG size);
void add_to_array(real_prec *addition, real_prec *out, ULONG size);
void add_to_array(complex_prec *addition, complex_prec *out, ULONG size);
void flip_sign(complex_prec *in, complex_prec *out, ULONG size);
void flip_sign(real_prec *in, real_prec *out, ULONG size);
void complexify_array(real_prec *in, complex_prec *out, ULONG size);
void real_part_array(complex_prec *in, real_prec *out, ULONG size);
void imaginary_part_array(complex_prec *in, real_prec *out, ULONG size);

void conjugate_array(complex_prec *in, complex_prec *out, ULONG size);
void times_i_array(complex_prec *in, complex_prec *out, ULONG size);
real_prec absolute_squared(complex_prec in);
void absolute_squared_array(complex_prec *in, real_prec *out, ULONG size);
void absolute_squared_array(complex_prec *in, complex_prec *out, ULONG size);

bool contains_nan(real_prec *in, ULONG size);
bool contains_nan(complex_prec *in, ULONG size);

std::vector<std::string> &split(const std::string &s, char delim,
                                std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);
