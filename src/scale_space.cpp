/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <cmath> // M_PI
#include "scale_space.hpp"

////////////////////////////////////////////////////////////////
// k-vector calculators
////////////////////////////////////////////////////////////////

real_prec k_squared(unsigned int i, unsigned int j, unsigned int k, real_prec L1, real_prec L2, real_prec L3,
                    unsigned int N1, unsigned int N2, unsigned int N3)
{
  real_prec k2=0.;
  auto kfac=static_cast<real_prec>(2.*M_PI);

  real_prec kx=0.;
  real_prec ky=0.;
  real_prec kz=0.;

  if (i<=N1/2) kx = kfac/L1*static_cast<real_prec>(i);
  else kx= -kfac/L1*static_cast<real_prec>(N1-i);

  if (j<=N2/2) ky= kfac/L2*static_cast<real_prec>(j);
  else ky= -kfac/L2*static_cast<real_prec>(N2-j);

  if (k<=N3/2) kz = kfac/L3*static_cast<real_prec>(k);
  else kz = -kfac/L3*static_cast<real_prec>(N3-k);

  k2=kx*kx+ky*ky+kz*kz;

  return(k2);
}

// EGP: Note: calc_kx, y and z were exactly the same functions! Replaced with calc_ki
real_prec calc_ki(unsigned int i, real_prec Li, unsigned int Ni) {
  auto kfac=static_cast<real_prec>(2.*M_PI/Li);
  real_prec ki=0.;

  if (i<=Ni/2)
    ki = kfac*static_cast<real_prec>(i);
  else
    ki = -kfac*static_cast<real_prec>(Ni-i);

  return(ki);
}

real_prec calc_kx(unsigned int i, real_prec L1, unsigned int N1) {
  return(calc_ki(i, L1, N1));
}

real_prec calc_ky(unsigned int j, real_prec L2, unsigned int N2) {
  return(calc_ki(j, L2, N2));
}

real_prec calc_kz(unsigned int k, real_prec L3, unsigned int N3) {
  return(calc_ki(k, L3, N3));
}
