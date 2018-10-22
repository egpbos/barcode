/*
 *  This file contains the implementation of various convenience functions
 *  used by the Planck LevelS package.
 *
 *  Copyright (C) 2002, 2003, 2004, 2005, 2006, 2007 Max-Planck-Society
 *  Authors: Martin Reinecke, Reinhard Hell
 */

/*
 * TODO: (EGP) clean up unused things
 */

// if we are using g++, check for version 3.0 or higher
#ifdef __GNUC__
#if (__GNUC__<3)
#error your C++ compiler is too old. g++ version 3.0 or higher is required.
#endif
#endif

#include <fstream>
#include <iostream>
#include "cxxutils.h"


std::string trim (const std::string &orig) {
  std::string::size_type p1=orig.find_first_not_of(" \t");
  if (p1==std::string::npos) return "";
  std::string::size_type p2=orig.find_last_not_of(" \t");
  return orig.substr(p1,p2-p1+1);
}


void parse_file (const std::string &filename, std::map<std::string, std::string> &dict)
  {
  int lineno=0;
  dict.clear();
    std::ifstream inp(filename.c_str());
  planck_assert(inp.good(),"Could not open parameter file "+filename);
  while (inp)
    {
    std::string line;
    getline(inp, line);
    ++lineno;
    line=line.substr(0,line.find_first_of("#"));
    line=trim(line);
    if (line.size()>0) {
      std::string::size_type eqpos=line.find("=");
      if (eqpos!=std::string::npos) {
        std::string key=trim(line.substr(0,eqpos)),
               value=trim(line.substr(eqpos+1,std::string::npos));
        if (key=="")
          std::cerr << "Warning: empty key in " << filename << ", line "
               << lineno << std::endl;
        else
          {
          if (dict.find(key)!=dict.end())
            std::cerr << "Warning: key " << key << " multiply defined in "
                 << filename << ", line " << lineno << std::endl;
          dict[key]=value;
          }
        }
      else
        std::cerr << "Warning: unrecognized format in " << filename << ", line "
             << lineno << ":\n" << line << std::endl;
      }
    }
  }


template<> void stringToData (const std::string &x, std::string &value)
{ value = trim(x); }

template<> void stringToData (const std::string &x, bool &value)
{
  if ( x=="F" || x=="f" || x=="n" || x=="N" || x=="false" || x==".false."
       || x=="FALSE" || x==".FALSE.")
    value=false;
  else if (x=="T" || x=="t" || x=="y" || x=="Y" || x=="true" || x==".true."
           || x=="TRUE" || x==".TRUE.")
    value=true;
  else
  {
    std::string error = std::string("conversion error in stringToData<bool>(\"")+x+"\")";
    throw Message_error (error);
  }
}


template<> std::string dataToString (const bool &x)
{ return x ? "T" : "F"; }
template<> std::string dataToString (const std::string &x)
{ return trim(x); }
template<> std::string dataToString (const float &x)
{
  std::ostringstream strstrm;
  strstrm << std::setprecision(8) << x;
  return trim(strstrm.str());
}
template<> std::string dataToString (const double &x)
{
  std::ostringstream strstrm;
  strstrm << std::setprecision(16) << x;
  return trim(strstrm.str());
}

