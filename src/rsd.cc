/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include "define_opt.h"
#include <cmath>
#include <iostream>
#include <stdexcept> // runtime_error
#include "pacman.hpp"

using namespace std;

void calc_pos_rsd(ULONG Npart, real_prec L3, real_prec xobs, real_prec yobs, real_prec zobs, const real_prec *x,
                  const real_prec *y,
                  const real_prec *z, const real_prec *vx, const real_prec *vy, const real_prec *vz, real_prec *xr, real_prec *yr,
                  real_prec *zr, real_prec ascale, real_prec Omega_M, real_prec Omega_L, bool planepar, bool periodic)
{
  // Units:
  // - xobs, yobs, zobs, x, y, z, xr, yr, zr in Mpc/h
  // - vx, vy, vz in km/s
  
  real_prec Omega_C = num_1 - Omega_M - Omega_L;        
  real_prec Hub = 100. * sqrt(Omega_M / ascale / ascale / ascale + Omega_L + Omega_C / ascale / ascale);           

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (ULONG i = 0; i < Npart; ++i)
  {
    // setting observation coordinates
    real_prec xpos = x[i] - xobs;
    real_prec ypos = y[i] - yobs;
    real_prec zpos = z[i] - zobs;

    real_prec rpos = sqrt(xpos*xpos+ypos*ypos+zpos*zpos);

    real_prec v_norm = 1. / Hub / ascale;

    if (!planepar) // planepar==false
    {
      real_prec ruxv = (xpos * vx[i] + ypos * vy[i] + zpos * vz[i]) / rpos * v_norm; // ru*v (v component in the direction of unit-r-vector)
      real_prec r_new = rpos + ruxv;
      xr[i] = (xpos / rpos) * r_new + xobs; // == r_new*sin(thetapos)*cos(phipos)+xobs, where phipos = atan2(ypos,xpos) and thetapos = acos(zpos/rpos)
      yr[i] = (ypos / rpos) * r_new + yobs; // == r_new*sin(thetapos)*sin(phipos)+yobs
      zr[i] = (zpos / rpos) * r_new + zobs; // == r_new*cos(thetapos)+zobs
    }
    else // planepar==true
    {
      real_prec ruxv = vz[i] * v_norm; // ru*v (v component in the direction of unit-r-vector)
      xr[i] = x[i];
      yr[i] = y[i];
      zr[i] = z[i] + ruxv;
    }

    if (periodic) {
      if (!planepar) {
        throw runtime_error("Periodic boundary conditions not yet implemented for non-plane-parallel RSDs!");
      } else {
        pacman_coordinate(&(zr[i]), L3);
      }
    }
  }
}

