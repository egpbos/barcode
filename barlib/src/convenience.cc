/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

#include "./define_opt.h"

// using namespace std;

// TODO: almost all functions are in fact array functions. Rename file to array.cc and move other stuff out (absolute_squared and split functions)

void copyArray(const real_prec *from, real_prec *to, ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG i = 0; i < size; i++)
    to[i] = from[i];
}

void copyArray(complex_prec *from, complex_prec *to, ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG i = 0; i < size; i++) {
    re(to[i]) = re(from[i]);
    im(to[i]) = im(from[i]);
  }
}

void sum_arrays(const real_prec *array1, const real_prec *array2, real_prec *out,
                ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG i = 0; i < size; i++)
    out[i] = array1[i] + array2[i];
}

void multiplyArrays(const real_prec *array1, const real_prec *array2, real_prec *out,
                    ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG i = 0; i < size; i++)
    out[i] = array1[i] * array2[i];
}

void multiply_factor_array(real_prec factor, const real_prec *array, real_prec *out,
                           ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG i = 0; i < size; i++)
    out[i] = factor * array[i];
}

void subtract_arrays(const real_prec *array1, const real_prec *array2, real_prec *out,
                     ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG i = 0; i < size; i++)
    out[i] = array1[i] - array2[i];
}

void add_to_array(real_prec addition, real_prec *out, ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG i = 0; i < size; i++)
    out[i] += addition;
}

void add_to_array(const real_prec *addition, real_prec *out, ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG i = 0; i < size; i++)
    out[i] += addition[i];
}

void add_to_array(const complex_prec *addition, complex_prec *out, ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG j = 0; j < size; j++) {
    re(out[j]) += re(addition[j]);
    im(out[j]) += im(addition[j]);
  }
}

void fill_one(real_prec *out, ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG i = 0; i < size; i++)
    out[i] = 1.0;
}

void fill_one(complex_prec *out, ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG i = 0; i < size; i++) {
    re(out[i]) = 1.0;
    im(out[i]) = 1.0;
  }
}

void fillZero(real_prec *out, ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG i = 0; i < size; i++)
    out[i] = 0.0;
}

void fillZero(ULONG *out, ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG i = 0; i < size; i++)
    out[i] = 0;
}

void fillZero(complex_prec *out, ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG i = 0; i < size; i++) {
    re(out[i])=0.0;
    im(out[i])=0.0;
  }
}

void flip_sign(const real_prec *in, real_prec *out, ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG j = 0; j < size; j++)
    out[j] = -in[j];
}

void flip_sign(const complex_prec *in, complex_prec *out, ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG j = 0; j < size; j++) {
    re(out[j]) = -re(in[j]);
    im(out[j]) = -im(in[j]);
  }
}

void complexify_array(const real_prec *in, complex_prec *out, ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG i = 0; i < size; i++) {
    re(out[i]) = in[i];
    im(out[i]) = 0.;
  }
}

void real_part_array(const complex_prec *in, real_prec *out, ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG i = 0; i < size; i++)
    out[i] = re(in[i]);
}

void imaginary_part_array(const complex_prec *in, real_prec *out, ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG i = 0; i < size; i++)
    out[i] = im(in[i]);
}


bool contains_nan(const real_prec *in, ULONG size) {
  for (ULONG i = 0; i < size; ++i)
    if (std::isnan(in[i])) {
      std::cout << "found a NaN at i = " << i << "!" << std::endl;
      return (true);
    }
  return false;
}

bool contains_nan(const complex_prec *in, ULONG size) {
  for (ULONG i = 0; i < size; ++i)
    if (std::isnan(re(in[i])) || std::isnan(im(in[i]))) {
      std::cout << "found a NaN at i = " << i << "!" << std::endl;
      return (true);
    }
  return false;
}


// Operations in Fourier space //

void conjugate_array(const complex_prec *in, complex_prec *out, ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG i = 0; i < size; i++) {
    re(out[i]) = re(in[i]);
    im(out[i]) = -im(in[i]);
  }
}

void times_i_array(const complex_prec *in, complex_prec *out, ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG i = 0; i < size; i++) {
    real_prec temp = re(in[i]);  // for in-place usage (out == in)
    re(out[i]) = -im(in[i]);
    im(out[i]) = temp;
  }
}


// TODO: This should go to math_funcs or a new scalar_math file
real_prec absolute_squared(const complex_prec in) {
  return re(in)*re(in) + im(in)*im(in);
}

void absolute_squared_array(const complex_prec *in, real_prec *out, ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG i = 0; i < size; i++)
    out[i] = absolute_squared(in[i]);
}

void absolute_squared_array(const complex_prec *in, complex_prec *out, ULONG size) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
  for (ULONG i = 0; i < size; i++)
    re(out[i]) = absolute_squared(in[i]);
}


// split string (http://stackoverflow.com/a/236803/1199693)
std::vector<std::string> &split(const std::string &s, char delim,
                                std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}
