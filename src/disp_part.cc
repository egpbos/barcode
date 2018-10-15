/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include "define_opt.h"

#include <iomanip>
#include <iostream>

//#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>

#include "fftw_array.h"

#include "math_funcs.h"

#include "convenience.h"

#include "BarcodeException.h"

/*

  F.S. Kitaura 2012
  E.G.P. Bos 2013-2017

*/

using namespace std;

bool planepar=false;// plane parallel approximation
bool periodic=true; // periodic boundary conditions


// ONLY USED IN Lag2Eul

void disp_part(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec L2, real_prec L3,
               real_prec d1, real_prec d2, real_prec d3, real_prec *posx, real_prec *posy, real_prec *posz,
               real_prec *psix, real_prec *psiy, real_prec *psiz, real_prec *dummyL, unsigned int facL, bool reggrid,
               gsl_rng *gBaseRand)
{
  if ((facL > 1) && (reggrid == false))
  {
    throw BarcodeException("facL > 1 and reggrid == false in disp_part; doesn't make sense, because it would then just put facL^3 particles in the same cell center positions.");
  }

  real_prec rx,ry,rz;

  ULONG nout=0;

  unsigned NL = facL*facL*facL;
  ULONG numpart = N1*N2*N3*NL;

  // Get grid center positions, or random positions if reggrid == false.
#ifdef MULTITHREAD_RNG
#pragma omp parallel for
#endif // MULTITHREAD_RNG
  for (unsigned i=0;i<N1;i++)
    for (unsigned j=0;j<N2;j++)
      for (unsigned k=0;k<N3;k++)
        for (unsigned ii=0; ii<NL; ii++)
        {
          ULONG jj = k + N3*(j+ N2*i);
          ULONG jjj = ii + NL*jj;

          if (reggrid==false)
          {
            rx=static_cast<real_prec>(gsl_rng_uniform(gBaseRand))*d1;
            ry=static_cast<real_prec>(gsl_rng_uniform(gBaseRand))*d2;
            rz=static_cast<real_prec>(gsl_rng_uniform(gBaseRand))*d3;
          }
          else
          {
            rx=static_cast<real_prec>(0.5)*d1;
            ry=static_cast<real_prec>(0.5)*d2;
            rz=static_cast<real_prec>(0.5)*d3;
          }

          real_prec posxi=d1*(real_prec(i))+rx;
          real_prec posyi=d2*(real_prec(j))+ry;
          real_prec poszi=d3*(real_prec(k))+rz;

          posx[jjj]=posxi;
          posy[jjj]=posyi;
          posz[jjj]=poszi;
        }

  // Interpolate the psi values to the particle positions (if reggrid == false)
  if (reggrid == false)
  {
    throw BarcodeException("This part of the code is not right, fix before use!");
    // below I explain why it's wrong
    fftw_array<real_prec> dri(numpart);
    interpolate_CIC(N1, N2, N3, L1, L2, L3, d1, d2, d3, posx, posy, posz, psix, numpart, dri);
    // WRONG: at this point we already add dri to posx:
    add_to_array(dri, posx, numpart);
    // WRONG: but in the next part we use this already updated posx to interpolate psiy!
    //        This interpolation should have been done before updating the pos arrays.
    interpolate_CIC(N1, N2, N3, L1, L2, L3, d1, d2, d3, posx, posy, posz, psiy, numpart, dri);
    add_to_array(dri, posy, numpart);
    interpolate_CIC(N1, N2, N3, L1, L2, L3, d1, d2, d3, posx, posy, posz, psiz, numpart, dri);
    add_to_array(dri, posz, numpart);
  }
  else
  {
  // If reggrid, just couple the psi values directly to the positions.
    add_to_array(psix, posx, numpart);
    add_to_array(psiy, posy, numpart);
    add_to_array(psiz, posz, numpart);
  }

  // Periodic boundary conditions
  if (periodic==true)
  {
#ifdef MULTITHREAD_RNG
#pragma omp parallel for
#endif // MULTITHREAD_RNG
    for (unsigned i=0;i<N1;i++)
      for (unsigned j=0;j<N2;j++)
        for (unsigned k=0;k<N3;k++)
          for (unsigned ii=0; ii<NL; ii++)
          {
            ULONG jj = k + N3*(j+ N2*i);
            ULONG jjj = ii + NL*jj;

            pacman_coordinate(&posx[jjj], L1);
            pacman_coordinate(&posy[jjj], L2);
            pacman_coordinate(&posz[jjj], L3);
          }
  }
  else
  {
  // Check for out of boundary particles (only if not periodic)
    bool bout=false;
#ifdef MULTITHREAD_RNG
#pragma omp parallel for reduction(+:nout)
#endif // MULTITHREAD_RNG
    for (unsigned i=0;i<N1;i++)
      for (unsigned j=0;j<N2;j++)
        for (unsigned k=0;k<N3;k++)
          for (unsigned ii=0; ii<NL; ii++)
          {
            ULONG jj = k + N3*(j+ N2*i);
            ULONG jjj = ii + NL*jj;

            if (posx[jjj]<0 || posy[jjj]<0 || posz[jjj]<0 || posx[jjj]>=L1 || posy[jjj]>=L1 || posz[jjj]>=L1)
            {
              nout++;
              bout=true;
            }

            dummyL[jjj]=1.0;
          }
    if (bout)
    {
      cout<<endl;
      cout<<"# of objects outside boundaries "<<nout<<endl;
    }
  }
}

