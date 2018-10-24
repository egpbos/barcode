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
#include "cxxutils.h"

class paramfile {
  private:
    std::map<std::string,std::string> params;

    std::string get_valstr(const std::string &key) const {
      auto loc = params.find(key);
      if (loc!=params.end()) return loc->second;
      throw Message_error ("Error: Cannot find the key \"" + key + "\".");
    }

  public:
    paramfile (const std::string &filename) {
        parse_file (filename, params);
    }

    template<typename T> T find (const std::string &key) const {
      T result;
      stringToData(get_valstr(key), result);
      return result;
    }
};

#endif
