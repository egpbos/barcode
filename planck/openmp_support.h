/*
 *  Copyright (C) 2005, 2006, 2007 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_OPENMP_SUPPORT_H
#define PLANCK_OPENMP_SUPPORT_H

#ifdef _OPENMP
#include <omp.h>
#endif

inline bool openmp_enabled()
  {
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
  }

inline int openmp_max_threads ()
  {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
  }

inline int openmp_thread_num ()
  {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
  }

/*! Calculates the range of indices between \a glo and \a ghi which
    must be processed by this thread and returns it in \a lo and \a hi.

    The indices \a ghi and \a hi are "one past the last real index",
    in analogy to the STL iterators. */
inline void openmp_calc_share (int glo, int ghi, int &lo, int &hi)
  {
#ifdef _OPENMP
  int nwork = ghi-glo;
  int nproc = omp_get_num_threads();
  int me = omp_get_thread_num();
  int nbase = nwork/nproc;
  int additional = nwork%nproc;
  lo = glo+me*nbase + ((me<additional) ? me : additional);
  hi = lo+nbase+(me<additional);
#else
  lo=glo; hi=ghi;
#endif
  }

#endif
