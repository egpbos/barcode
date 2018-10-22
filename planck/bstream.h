#pragma once

#include <utility>

class bofstream: public std::ofstream
  {
  public:
    /*! */
    bofstream (const char *fname)
      : std::ofstream(fname,std::ios::binary) {}

    template<typename T> bofstream &put (const T *data, unsigned long num)
      {
      write (reinterpret_cast<const char *> (data), static_cast<std::streamsize>(num*sizeof(T)));
      return *this;
      }
  };

class bifstream: public std::ifstream
  {
  public:
    /*! */
    bifstream (const char *fname)
      : std::ifstream(fname,std::ios::binary) {}

    template<typename T> bifstream &get (T *data, unsigned long num) {
      read (reinterpret_cast<char *> (data), static_cast<std::streamsize>(num*sizeof(T)));
      return *this;
    }
  };
