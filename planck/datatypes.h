/*
 *  This file defines various platform-independent data types.
 *  If any of the requested types is not available, compilation aborts
 *  with an error (unfortunately a rather obscure one).
 *
 *  Copyright (C) 2004 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef PLANCK_DATATYPES_H
#define PLANCK_DATATYPES_H

#include <string>
#include "message_error.h"

// Template magic to select the proper data types. These templates
// should not be used outside this file.

template<typename T> inline const char *type2typename ()
  { return "unknown type"; }
template<> inline const char *type2typename<signed char> ()
  { return "signed char"; }
template<> inline const char *type2typename<unsigned char> ()
  { return "unsigned char"; }
template<> inline const char *type2typename<short> ()
  { return "short"; }
template<> inline const char *type2typename<unsigned short> ()
  { return "unsigned short"; }
template<> inline const char *type2typename<int> ()
  { return "int"; }
template<> inline const char *type2typename<unsigned int> ()
  { return "unsigned int"; }
template<> inline const char *type2typename<long> ()
  { return "long"; }
template<> inline const char *type2typename<unsigned long> ()
  { return "unsigned long"; }
template<> inline const char *type2typename<long long> ()
  { return "long long"; }
template<> inline const char *type2typename<unsigned long long> ()
  { return "unsigned long long"; }
template<> inline const char *type2typename<float> ()
  { return "float"; }
template<> inline const char *type2typename<double> ()
  { return "double"; }
template<> inline const char *type2typename<bool> ()
  { return "bool"; }
template<> inline const char *type2typename<std::string> ()
  { return "std::string"; }

#endif
