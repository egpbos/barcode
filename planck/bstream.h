#pragma once

#include <utility>

typedef unsigned char byte;
const bool big_endian=false;
#define ULONG unsigned long

template<ULONG size> inline void byteswap_helper (byte *)
  {
  bool Error_unspecialized_template_used[-(size*size)];
  // compile time error would be nice here ...
  }
template<> inline void byteswap_helper<1> (byte *)
  {}
template<> inline void byteswap_helper<2> (byte *val)
  {
  using namespace std;
  swap (val[0],val[1]);
  }
template<> inline void byteswap_helper<4> (byte *val)
  {
  using namespace std;
  swap (val[0],val[3]); swap (val[1],val[2]);
  }
template<> inline void byteswap_helper<8> (byte *val)
  {
  using namespace std;
  swap (val[0],val[7]); swap (val[1],val[6]);
  swap (val[2],val[5]); swap (val[3],val[4]);
  }

template<typename T> inline void byteswap (T& val)
  {
  byteswap_helper<sizeof(T)> (reinterpret_cast<byte *> (&val));
  }

const bool file_is_lsb=big_endian, file_is_msb=!big_endian,
           file_is_natural=false;

class bofstream: public std::ofstream
  {
  private:
    bool doswap;

  public:
    /*! */
    bofstream (const char *fname, bool doswap_)
      : std::ofstream(fname,std::ios::binary), doswap(doswap_) {}

    template<typename T> bofstream &operator<< (const T &data)
      {
      if (doswap)
        {
        T tmp = data;
        byteswap (tmp);
        write (reinterpret_cast<const char *> (&tmp), sizeof(T));
        }
      else
        write (reinterpret_cast<const char *> (&data), sizeof(T));
      return *this;
      }
    template<typename T> bofstream &put (const T *data, ULONG num)
      {
      if (doswap)
        {
        for (ULONG m=0; m<num; ++m)
          {
	  T tmp=data[m];
          byteswap (tmp);
          write (reinterpret_cast<const char *> (&tmp), sizeof(T));
          }
        }
      else
        write (reinterpret_cast<const char *> (data), static_cast<std::streamsize>(num*sizeof(T)));
      return *this;
      }
  };

class bifstream: public std::ifstream
  {
  private:
    bool doswap;

  public:
    /*! */
    bifstream (const char *fname, bool doswap_)
      : std::ifstream(fname,std::ios::binary), doswap(doswap_) {}

    template<typename T> bifstream &operator>> (T &data) {
      read (reinterpret_cast<char *> (&data), sizeof(T));
      if (doswap)
        byteswap (data);
      return *this;
    }
    template<typename T> bifstream &get (T *data, ULONG num) {
      read (reinterpret_cast<char *> (data), static_cast<std::streamsize>(num*sizeof(T)));
      if (doswap)
        for (ULONG m=0; m<num; ++m)
          byteswap (data[m]);
      return *this;
    }
  };
