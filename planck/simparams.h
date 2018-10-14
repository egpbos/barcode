/*
 *  Class for storing parameter information for later use
 *
 *  Copyright (C) 2003 Max-Planck-Society
 *  Authors: Reinhard Hell, Martin Reinecke
 */

#ifndef PLANCK_SIMPARAMS_H
#define PLANCK_SIMPARAMS_H

#include <string>
#include <vector>
#include <iostream>
#include "cxxutils.h"
class fitshandle;

class simparams
  {
  private:
    class Param
      {
      public:
        std::string key, shortkey, value, comment;

        Param (const std::string &Key, const std::string &Shortkey,
               const std::string &Value, const std::string &Comment)
          : key(Key), shortkey(Shortkey), value(Value), comment(Comment) {}
      };

    std::vector<Param> paramMap;
    std::vector<std::string> source_files;
    std::vector<int> hdus;

  public:
    void add_comment (const std::string &comment)
      { paramMap.push_back(Param("","","",comment)); }
    template<typename T> void add(const std::string &key,
      const std::string &shortkey, const T &value, const std::string &comment)
      {
      paramMap.push_back(Param(key, shortkey, dataToString(value), comment));
      }
    template<typename T> void add(const std::string &key,
      const std::string &shortkey, const T &value)
      {
      paramMap.push_back(Param(key, shortkey, dataToString(value), ""));
      }

    void add_source_file (const std::string &filename, int hdu=2)
      {
      source_files.push_back(filename);
      hdus.push_back(hdu);
      }

    void add_keys (std::ostream &os) const;
    void add_keys (fitshandle &out) const;
  };

#endif
