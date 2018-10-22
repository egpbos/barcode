/*
 *  Class for parsing parameter files
 *
 *  Copyright (C) 2003, 2004, 2005 Max-Planck-Society
 *  Authors: Martin Reinecke, Reinhard Hell
 */

#ifndef PLANCK_PARAMFILE_H
#define PLANCK_PARAMFILE_H

#include <map>
#include <string>
#include <iostream>
#include "cxxutils.h"

class paramfile
  {
  private:
    typedef std::map<std::string,std::string> params_type;
    params_type params;
    bool verbose;

    std::string get_valstr(const std::string &key) const
      {
      params_type::const_iterator loc=params.find(key);
      if (loc!=params.end()) return loc->second;
      throw Message_error ("Error: Cannot find the key \"" + key + "\".");
      }

  public:
    paramfile (const std::string &filename, bool verbose_=true)
      : verbose(verbose_)
      { parse_file (filename, params); }

    paramfile (const params_type &par)
      : params (par), verbose(true)
      {}

    bool param_present(const std::string &key) const
      { return (params.find(key)!=params.end()); }

    template<typename T> T find (const std::string &key) const
      {
      T result;
      stringToData(get_valstr(key),result);
      //EGP if (verbose)
      //  std::cout << "Parser: " << key << " = " << dataToString(result)
            //      << std::endl;
      return result;
      }
    template<typename T> T find
      (const std::string &key, const T &deflt)
      {
      if (param_present(key)) return find<T>(key);
      if (verbose)
      //  std::cout << "Parser: " << key << " = " << dataToString(deflt)
          //        << " <default>" << std::endl;
      params[key]=dataToString(deflt);
      return deflt;
      }

  };

#endif
