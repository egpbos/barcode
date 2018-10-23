/*! \file cxxutils.h
 *  Various convenience functions used by the Planck LevelS package.
 *
 *  Copyright (C) 2002, 2003, 2004, 2005, 2006, 2007 Max-Planck-Society
 *  \author Martin Reinecke \author Reinhard Hell
 */

#ifndef PLANCK_CXXUTILS_H
#define PLANCK_CXXUTILS_H

#include <algorithm>
#include <string>
#include <map>
#include <sstream>
#include <iomanip>
#include "message_error.h"


//! Throws a Message_error containing \a msg if \a testval is false.
inline void planck_assert (bool testval, const std::string &msg)
  {
  if (testval) return;
  throw Message_error ("Assertion failed: "+msg);
  }


/*! \defgroup stringutilsgroup String handling helper functions */
/*! \{ */

//! Returns the string \a orig without leading and trailing whitespace.
std::string trim (const std::string &orig);

//! Parses the file \a filename and returns the key/value pairs in \a dict.
void parse_file (const std::string &filename,
  std::map<std::string,std::string> &dict);

/*! \} */


//! Returns a string containing the text representation of \a x.
/*! Care is taken that no information is lost in the conversion. */

template<typename T> std::string dataToString (const T &x)
{
  std::ostringstream strstrm;
  strstrm << x;
  return trim(strstrm.str());
}


//! Reads a value of a given datatype from a string

template<typename T> void stringToData (const std::string &x, T &value)
{
  std::istringstream strstrm(x);
  strstrm >> value;
  if (!strstrm) throw Message_error("conversion error in stringToData");

  std::string rest;
  strstrm >> rest;
  if (rest.length()>0) throw Message_error("conversion error in stringToData (rest data left)");
}


#endif
