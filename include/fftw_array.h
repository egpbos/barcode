/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#pragma once

#include "fftw3.h"
//#include <iostream>
#include "define_opt.h"
#include <cstddef>  // std::ptrdiff_t, std::size_t
/* This template allows to avoid writing fftw_free */


/*
 * Trying to figure out a fix for the behavior of fftw_array as it is on 14 April 2017 (and before).
 * Specifically, when compiling with clang (possibly with gcc as well), I get a warning about the implicit conversion:
 *   /Users/pbos/code/barcode/EqSolvers.cc:582:21: warning: implicit conversion changes signedness: 'unsigned long' to 'long' [-Wsign-conversion]
 *     m2v[i]=LapPhivx[i]*LapPhivy[i]-LapPhivxy[i]*LapPhivxy[i]+LapPhivx[i]*LapPhivz[i]-LapPhivxz[i]*LapPhivxz[i]+LapPhivy[i]*LapPhivz[i]-LapPhivyz[i]*LapPhivyz[i];
 *            ~~~~~~~~ ^
 * This happens with all fftw_arrays when an unsigned integer is used as an index.
 *
 * It is probably caused by the built-in subscript operator: http://en.cppreference.com/w/cpp/language/operator_member_access
 * Since fftw_array has no operator[] (because it is neither a pointer nor an array, so operator[] cannot take fftw_array as operand),
 * the implicit conversion operator `operator T*()` (see fftw_array.h) kicks in as the only available option. This returns the `data`
 * member which is a pointer of type T*.
 *
 * The builtin subscript operator T& operator[](T*, std::ptrdiff_t) then converts the index to a std::ptrdiff_t, which though
 * implementation dependent is usually a signed integer of some size, since it needs to be able to do negative differences too.
 * This causes the signedness warning.
 *
 * A solution would be to overload the subscript operator for unsigned types.
 *
 * However, fftw_array is often passed to functions using its data member only. In this case, we cannot overload the operator from the
 * fftw_array class.
 *
 * Can we have the class template additionally create an overloaded global function T& operator[](T*, std::ptrdiff_t)? I.e. one outside
 * of the class fftw_array, but defined by the use of fftw_array<T> somewhere in the code...
 *
 * Do we have to overload the template function? :D
 *
 * 15 April:
 * It seems like there are no changed signedness warnings for the real_prec* arrays. Perhaps the compiler only warns about this
 * for user defined classes. Let's just overload the operator here and see what happens.
 */


template<typename T> class fftw_array
{
  // N.B.: fftw_array does not copy/assign properly! The data pointers will be copied, meaning that the new array will point to the same memory slots! This gives errors when destroying the objects when they go out of scope.
public:
  T *data;

  fftw_array(ULONG size)
  {
#ifdef SINGLE_PREC
    data = static_cast<T *>( fftwf_malloc(size*sizeof(T)) );
#endif
#ifdef DOUBLE_PREC
    data = static_cast<T *>( fftw_malloc(size*sizeof(T)) );
#endif
  }
  ~fftw_array()
  {
#ifdef SINGLE_PREC
    fftwf_free(data);
#endif
#ifdef DOUBLE_PREC
    fftw_free(data);
#endif
  }

  // "implicit conversion operator": works for functions, but not for templates. For the latter case write: <array_name>.data
  operator T*()
  {
    return data;
  }

  // overloaded subscript operators, see story at the top of the file
  T& operator[](std::size_t N) {
    // *(ptr+N) is equivalent to ptr[N], but since ptr[N] will again call the builtin subscript operator, we don't
    // want to use that, since it will again need casting to signed integer (std::ptrdiff_t) and thus will limit the
    // range of the index to half the total possible range.
    return *(data + N);
  }

};

//template<typename T> class fftw_array_debug
//{
//  // N.B.: fftw_array_debug does not copy/assign properly! The data pointers will be copied, meaning that the new array will point to the same memory slots! This gives errors when destroying the objects when they go out of scope.
//public:
//  T *data;
//
//  fftw_array_debug(ULONG size)
//  {
//    std::cout << "creating fftw_array_debug with address " << this << " and data address " << data << std::endl;
//
//#ifdef SINGLE_PREC
//    data = static_cast<T *>( fftwf_malloc(size*sizeof(T)) );
//#endif
//#ifdef DOUBLE_PREC
//    data = static_cast<T *>( fftw_malloc(size*sizeof(T)) );
//#endif
//  }
//  ~fftw_array_debug()
//  {
//    std::cout << "destroying fftw_array_debug with address " << this << std::endl;
//#ifdef SINGLE_PREC
//    fftwf_free(data);
//#endif
//#ifdef DOUBLE_PREC
//    fftw_free(data);
//#endif
//  }
//
//  // implicit conversion: works for functions, but not for templates. For the latter case write: <array_name>.data
//  operator T*()
//  {
//    return data;
//  }
//
//};
