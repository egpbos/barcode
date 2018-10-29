/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include "define_opt.h"
#include "struct_main.h"

#include <math.h>
#include <iomanip>

#include <gsl/gsl_integration.h>

#include "fftw_array.h"

#include "math_funcs.h"

#include "convenience.h"
#include "BarcodeException.h"

// EGP: time.h for benchmark tests
#include <sys/time.h>

using namespace std;

void overdens(unsigned int N1, unsigned int N2, unsigned int N3, const real_prec *in, real_prec *out)
{
  ULONG N=N1*N2*N3;

  double nmeanD=0.;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:nmeanD)
#endif // MULTITHREAD
  for(ULONG i=0;i<N;i++)
    nmeanD+=static_cast<double>(in[i]);
  real_prec nmean=static_cast<real_prec>(nmeanD)/static_cast<real_prec>(N);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (ULONG i=0;i<N;i++)
    out[i]=in[i]/nmean-num_1;
}

void getDensity_NGP(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec L2, real_prec L3,
                    real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,
                    const real_prec *xp, const real_prec *yp, const real_prec *zp, const real_prec *Par_mass, ULONG N_OBJ, real_prec *delta)
{

  //EGP    real_prec xc, yc, zc;
  //EGP    real_prec dx, dy, dz, tx, tz, ty;
  unsigned i, j, k;//EGP, ii, jj, kk, n;


#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (i=0;i<N1; i++)
    for (j=0; j<N2; j++)
      for (k=0; k<N3; k++){
        delta[k+N3*(j+N2*i)] = 0.;  //-1 if we want to calculate overdensity
      }

  ULONG NLOSS=0;
#ifdef MULTITHREAD
#pragma omp parallel for private(i,j,k)
#endif // MULTITHREAD
  for (ULONG n=0; n<N_OBJ; n++){ // EGP: declared n here instead of above as int

    //check if particle is in selected Domain, else discard it
    if((xp[n]>=min1 && xp[n]<min1+L1) && (yp[n]>=min2 && yp[n]<min2+L2) && (zp[n]>=min3 && zp[n]<min3+L3))
    {

      i = static_cast<unsigned>(floor((xp[n]-min1)/d1)); // indices of the cell of the particle
      j = static_cast<unsigned>(floor((yp[n]-min2)/d2));
      k = static_cast<unsigned>(floor((zp[n]-min3)/d3));


      i = static_cast<unsigned>(fmod(real_prec(i),real_prec(N1)));
      j = static_cast<unsigned>(fmod(real_prec(j),real_prec(N2)));
      k = static_cast<unsigned>(fmod(real_prec(k),real_prec(N3)));


      real_prec mass=Par_mass[n];
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      delta[k+N3*(j+N2*i)]    +=mass;

    }
    //else NLOSS++;
  }
  if(NLOSS!=0) cout << " >>> mass assignment found "<<NLOSS<<" particles outside mesh boundary...."<<endl<<endl;
}

void getDensity_CIC(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3, const real_prec *xp, const real_prec *yp, const real_prec *zp,
                    const real_prec *Par_mass, ULONG N_OBJ, real_prec *delta, bool weightmass)
{
#define DELTA(i,j,k) delta[k[2]+N3*(j[1]+N2*i[0])]

  fillZero(delta, N1*N2*N3);

  ULONG NLOSS=0;
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (ULONG n=0; n<N_OBJ; n++){
    real_prec dx[3], tx[3];
    ULONG i[3], ii[3];

   //check if particle is in selected Domain, else discard it
    if((xp[n]>=min1 && xp[n]<min1+L1) && (yp[n]>=min2 && yp[n]<min2+L2) && (zp[n]>=min3 && zp[n]<min3+L3))
    {	
      getCICcells(N1, N2, N3, L1, L2, L3, d1, d2, d3, xp[n], yp[n], zp[n], i, ii);
      getCICweights(L1, L2, L3, d1, d2, d3, xp[n], yp[n], zp[n], i, dx, tx);

      // Take care this assumes periodic boundary conditions. This conserves Mass when using FFTs to deconvolve with CIC kernel, as the FFT assumes periodicity
      real_prec mass=num_1;
      if (weightmass==true)
        mass=Par_mass[n];

#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(i,i,i)    += mass*tx[0]*tx[1]*tx[2];
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(ii,i,i)   += mass*dx[0]*tx[1]*tx[2];
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(i,ii,i)   += mass*tx[0]*dx[1]*tx[2];
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(i,i,ii)   += mass*tx[0]*tx[1]*dx[2];
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(ii,ii,i)  += mass*dx[0]*dx[1]*tx[2];
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(ii,i,ii)  += mass*dx[0]*tx[1]*dx[2];
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(i,ii,ii)  += mass*tx[0]*dx[1]*dx[2];
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(ii,ii,ii) += mass*dx[0]*dx[1]*dx[2];
    }
    //else NLOSS++;	
  }
  if(NLOSS!=0) cout << " >>> ARGO: Mass assignment found "<<NLOSS<<" particles outside mesh boundary...."<<endl<<endl;

#undef DELTA
}


void getDensity_TSC(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec L2, real_prec L3,
                    real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,
                    const real_prec *xp, const real_prec *yp, const real_prec *zp, const real_prec *Par_mass, ULONG N_OBJ, real_prec *delta)
{

  real_prec xc, yc, zc;
  real_prec dx, dy, dz, hx0, hz0, hy0, hxp1,hyp1,hzp1, hxm1,hym1,hzm1;
  unsigned i, j, k, ii, jj, kk,iii,jjj,kkk;//EGP, n;

#define DELTA(i,j,k) delta[k+N3*(j+N2*i)]

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (i=0;i<N1; i++)
    for (j=0; j<N2; j++)
      for (k=0; k<N3; k++){
        DELTA(i,j,k) = 0.;  //-1 if we want to calculate overdensity
      }

  ULONG NLOSS=0;

#ifdef MULTITHREAD
#pragma omp parallel for private(xc, yc, zc, dx, dy, dz, hx0, hz0, hy0, hxp1,hyp1,hzp1, hxm1,hym1,hzm1, i, j, k, ii, jj, kk,iii,jjj,kkk)
#endif // MULTITHREAD
  for (ULONG n=0; n<N_OBJ; n++){// EGP: declared n here instead of above as int

    //check if particle is in selected Domain, else discard it
    if((xp[n]>=min1 && xp[n]<=min1+L1) && (yp[n]>=min2 && yp[n]<=min2+L2) && (zp[n]>=min3 && zp[n]<=min3+L3))
    {

      i = static_cast<unsigned>(floor((xp[n]-min1)/d1)); // indices of the cell of the particle
      j = static_cast<unsigned>(floor((yp[n]-min2)/d2));
      k = static_cast<unsigned>(floor((zp[n]-min3)/d3));


      i = static_cast<unsigned>(fmod(real_prec(i),real_prec(N1)));
      j = static_cast<unsigned>(fmod(real_prec(j),real_prec(N2)));
      k = static_cast<unsigned>(fmod(real_prec(k),real_prec(N3)));

      ii = static_cast<unsigned>(fmod(real_prec(i+1),real_prec(N1)));
      jj = static_cast<unsigned>(fmod(real_prec(j+1),real_prec(N2)));
      kk = static_cast<unsigned>(fmod(real_prec(k+1),real_prec(N3)));

      iii= static_cast<unsigned>(fmod(real_prec(i-1+N1),real_prec(N1)));
      jjj= static_cast<unsigned>(fmod(real_prec(j-1+N2),real_prec(N2)));
      kkk= static_cast<unsigned>(fmod(real_prec(k-1+N3),real_prec(N3)));


      //printf("%f, %f, %f\n", xp[n], yp[n], zp[n]);

      xc = static_cast<real_prec>(i+0.5); // centers of the cells
      yc = static_cast<real_prec>(j+0.5);
      zc = static_cast<real_prec>(k+0.5);

      dx = (xp[n]-min1)/d1 - xc; // distance of particle to center of the cell
      dy = (yp[n]-min2)/d2 - yc;
      dz = (zp[n]-min3)/d3 - zc;



      hx0=static_cast<real_prec>(0.75-dx*dx); // fraction of particle assigned
      hy0=static_cast<real_prec>(0.75-dy*dy); // fraction of particle assigned
      hz0=static_cast<real_prec>(0.75-dz*dz); // fraction of particle assigned

      hxp1=static_cast<real_prec>(0.5*(0.5+dx)*(0.5+dx)); // fraction of particle assigned
      hyp1=static_cast<real_prec>(0.5*(0.5+dy)*(0.5+dy)); // fraction of particle assigned
      hzp1=static_cast<real_prec>(0.5*(0.5+dz)*(0.5+dz)); // fraction of particle assigned

      hxm1=static_cast<real_prec>(0.5*(0.5-dx)*(0.5-dx)); // fraction of particle assigned
      hym1=static_cast<real_prec>(0.5*(0.5-dy)*(0.5-dy)); // fraction of particle assigned
      hzm1=static_cast<real_prec>(0.5*(0.5-dz)*(0.5-dz)); // fraction of particle assigned



      real_prec mass=Par_mass[n];
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(i,j,k)            += mass*hx0*hy0*hz0;
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(ii,jj,kk)         += mass*hxp1*hyp1*hzp1;
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(iii,jjj,kkk)         += mass*hxm1*hym1*hzm1;

#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(ii,j,k)           += mass*hxp1*hy0*hz0;
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(i,jj,k)           += mass*hx0*hyp1*hz0;
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(i,j,kk)           += mass*hx0*hy0*hzp1;

#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(ii,jj,k)          += mass*hxp1*hyp1*hz0;
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(ii,j,kk)          += mass*hxp1*hy0*hzp1;
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(i,jj,kk)          += mass*hx0*hyp1*hzp1;

#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(iii,jj,kk)           += mass*hxm1*hyp1*hzp1;
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(ii,jjj,kk)          += mass*hxp1*hym1*hzp1;
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(ii,jj,kkk)          += mass*hxp1*hyp1*hzm1;

#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(iii,jjj,kk)          += mass*hxm1*hym1*hzp1;
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(iii,jj,kkk)          += mass*hxm1*hyp1*hzm1;
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(ii,jjj,kkk)          += mass*hxp1*hym1*hzm1;

#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(i,jjj,kkk)           += mass*hx0*hym1*hzm1;
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(iii,j,kkk)          += mass*hxm1*hy0*hzm1;
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(iii,jjj,k)          += mass*hxm1*hym1*hz0;

#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(i,j,kkk)          += mass*hx0*hy0*hzm1;
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(i,jjj,k)          += mass*hx0*hym1*hz0;
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(iii,j,k)          += mass*hxm1*hy0*hz0;

#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(i,jj,kkk)           += mass*hx0*hyp1*hzm1;
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(ii,j,kkk)          += mass*hxp1*hy0*hzm1;

#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(iii,jj,k)           += mass*hxm1*hyp1*hz0;
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(ii,jjj,k)          += mass*hxp1*hym1*hz0;

#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(i,jjj,kk)           += mass*hx0*hym1*hzp1;
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
      DELTA(iii,j,kk)          += mass*hxm1*hy0*hzp1;
    }
  }
  if(NLOSS!=0) cout << " >>> mass assignment found "<<NLOSS<<" particles outside mesh boundary...."<<endl<<endl;
  // cout <<"TOTAL NUMBER OF OBJECTS = "<< N_OBJ <<endl;
}

real_prec SPH_kernel_3D(real_prec r, real_prec h)
{
  // N.B.: when using this for density estimation, you need to normalize the
  // result afterwards with V/N! See e.g. getDensity_SPH. Same goes for
  // grad_SPH_kernel_3D.
  
  // Monaghan kernel W_4
  real_prec result = 0.;
  real_prec q = r/h;
  if (q < 0.) {
    throw BarcodeException("In SPH_kernel_3D: r may not be smaller than zero!");
  }
  else if (q <= 1.)
    result = 1./M_PI/(h*h*h) * (1 - 3./2*q*q + 3./4*q*q*q);
  else if (q <= 2.)
    result = 1./M_PI/(h*h*h) * (1./4 * gsl_pow_3(2.-q));

  return(result);
}


//void getDensity_SPH(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,real_prec *xp,real_prec *yp,real_prec *zp, real_prec *Par_mass, ULONG N_OBJ, real_prec *delta, bool weightmass, real_prec kernel_h, bool old_cell_index) {
void getDensity_SPH(ULONG, ULONG, ULONG, real_prec, real_prec, real_prec, real_prec, real_prec, real_prec, real_prec, real_prec, real_prec, real_prec *, real_prec *,real_prec *, real_prec *, ULONG, real_prec *, bool, real_prec, bool) {
  std::cout << "This signature is deprecated, the old_cell_index parameter is no longer used. Use Mercurial revision 938 if you want to use the version with the old cell index for testing." << std::endl;
}

void getDensity_SPH(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,
                    const real_prec *xp,
                    const real_prec *yp,
                    const real_prec *zp,
                    const real_prec *Par_mass, ULONG N_OBJ, real_prec *delta, bool weightmass, real_prec kernel_h)
{
  // Note: current implementation assumes periodic boundary conditions!
  //
  // Note: all normalization stuff removed; it is unnecessary. This function now
  // just gives rho.
  //
  // Warning: don't use with kernel_h larger than N{1,2,3}/4, otherwise the kernel diameter will be larger than the box width in that direction!
  //          When using this function inside Barcode using the supplied INIT_PARAMS initialization function, this condition will always be enforced.
  // Warning 2: Make sure N{1,2,3} are smaller than LONG_MAX, since they need to be cast to signed integers.

  fillZero(delta, N1*N2*N3);

  // First determine the reach of the kernel, which determines how many cells
  // every particle needs to loop over to check for contributions. This is
  // based on the kernel radius ( == 2*kernel_h ).
  int reach1 = static_cast<int>(2*kernel_h/d1) + 1, reach2 = static_cast<int>(2*kernel_h/d2) + 1, reach3 = static_cast<int>(2*kernel_h/d3) + 1;
  // cout << "reach1: " << reach1 << endl;

  // For normalization later on:
  //real_prec mass_total = 0;

  ULONG NLOSS=0;
#ifdef MULTITHREAD
//#pragma omp parallel for reduction(+:mass_total)
#pragma omp parallel for
#endif // MULTITHREAD
  for (ULONG n=0; n<N_OBJ; n++)
  {
    // check if particle is in selected domain, else discard it
    if((xp[n] >= min1 && xp[n] < min1+L1) && (yp[n] >= min2 && yp[n] < min2+L2) && (zp[n] >= min3 && zp[n] < min3+L3))
    {
      real_prec mass = num_1;
      if (weightmass == true)
        mass = Par_mass[n];
      //mass_total += mass;

      // Determine central cell index where particle resides
      ULONG ix = static_cast<ULONG>(xp[n]/d1);
      ULONG iy = static_cast<ULONG>(yp[n]/d2);
      ULONG iz = static_cast<ULONG>(zp[n]/d3);
      // Central cell position:
      real_prec ccx = (static_cast<real_prec>(ix) + 0.5)*d1;
      real_prec ccy = (static_cast<real_prec>(iy) + 0.5)*d2;
      real_prec ccz = (static_cast<real_prec>(iz) + 0.5)*d3;

      // Loop over surrounding gridcells (including central cell itself) within
      // kernel radius.
      for(int i1 = -reach1; i1 <= reach1; ++i1)
        for(int i2 = -reach2; i2 <= reach2; ++i2)
          for(int i3 = -reach3; i3 <= reach3; ++i3)
          {
            // Cell position (relative to the central cell):
            real_prec cx = ccx + static_cast<real_prec>(i1)*d1;
            real_prec cy = ccy + static_cast<real_prec>(i2)*d2;
            real_prec cz = ccz + static_cast<real_prec>(i3)*d3;
            // Cell index, taking into account periodic boundary conditions:
            ULONG kx, ky, kz;
            kx = (static_cast<ULONG>(static_cast<long>(N1) + i1) + ix) % N1;  // checked in testcodes/signed_unsigned_periodic.cpp that signedness casts here doesn't cause bugs
            ky = (static_cast<ULONG>(static_cast<long>(N2) + i2) + iy) % N2;
            kz = (static_cast<ULONG>(static_cast<long>(N3) + i3) + iz) % N3;
            // The above casts are necessary to ensure implicit signedness casts don't cause trouble.
            // We assume here that kernel_h (and thus reach) <= N{1,2,3}/4 and N1 < LONG_MAX (the latter is true in Barcode, since Ni are usually uint, not ulong).
            ULONG index = kz + N3*(ky + N2*kx);

            real_prec diff_x = xp[n] - cx;
            real_prec diff_y = yp[n] - cy;
            real_prec diff_z = zp[n] - cz;
            real_prec r = sqrt(diff_x*diff_x + diff_y*diff_y + diff_z*diff_z);
            
            if (r/kernel_h <= 2.)
#ifdef MULTITHREAD
#pragma omp atomic
#endif // MULTITHREAD
              //delta[index] += normalize * SPH_kernel_3D(r, kernel_h) * mass;
              delta[index] += SPH_kernel_3D(r, kernel_h) * mass;
          }
    }
    else
      ++NLOSS;
  }
  if (NLOSS>0)
    cout << "Lost " << NLOSS << " particles in getDensity_SPH." << endl;

  // Normalization given that we want the integral over density to be V (left-
  // hand side) and the integral over all particle kernels to be N (right-hand
  // side).
  // FIXME: this really isn't necessary at all. It's only necessary if you
  // expect this function to give you rho/rho_c. If you want it to give you
  // simply rho (in units of mass/volume), then you don't need it!
  //real_prec normalize = L1*L2*L3/static_cast<real_prec>(N1*N2*N3); // Note: actually N1*N2*N3 was also wrong, should be N_OBJ!
  //real_prec normalize = L1*L2*L3/mass_total;

//#ifdef MULTITHREAD
//#pragma omp parallel for
//#endif // MULTITHREAD
  //for (ULONG n = 0; n < N_OBJ; ++n)
    //delta[n] *= normalize;

}


void cellbound(unsigned int N1, unsigned int N2, unsigned int N3, real_prec *v1, real_prec *v2, real_prec *v3)
{
  ULONG  N=N1*N2*N3;

  /* start: interpolate from cell center to cell boundaries */
  fftw_array<real_prec> vx(N),vy(N),vz(N);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (unsigned i=0;i<N1; i++)
    for (unsigned j=0; j<N2; j++)
      for (unsigned k=0; k<N3; k++)
      {
        ULONG l=k+N3*(j+N2*i);
        ULONG m=k-1+N3*(j-1+N2*(i-1));

        if (i>0 && j>0 && k>0)
        {
          vx[l]=static_cast<real_prec>(0.5)*(v1[m]+v1[l]);
          vy[l]=static_cast<real_prec>(0.5)*(v2[m]+v2[l]);
          vz[l]=static_cast<real_prec>(0.5)*(v3[m]+v3[l]);
        }

        /* periodic boundary conditions */

        long ii=k-1+N3*(j-1+N2*(N1-1));
        long jj=k-1+N3*(N2-1+N2*(i-1));
        long kk=N3-1+N3*(j-1+N2*(i-1));
        long ij=k-1+N3*(N2-1+N2*(N1-1));
        long ik=N3-1+N3*(j-1+N2*(N1-1));
        long jk=N3-1+N3*(N2-1+N2*(i-1));
        long ijk=N3-1+N3*(N2-1+N2*(N1-1));

        if (i==0 && j>0 && k>0)
        {
          vx[l]=static_cast<real_prec>(0.5)*(v1[ii]+v1[l]);
          vy[l]=static_cast<real_prec>(0.5)*(v2[ii]+v2[l]);
          vz[l]=static_cast<real_prec>(0.5)*(v3[ii]+v3[l]);
        }
        if (i==0 && j==0 && k>0)
        {
          vx[l]=static_cast<real_prec>(0.5)*(v1[ij]+v1[l]);
          vy[l]=static_cast<real_prec>(0.5)*(v2[ij]+v2[l]);
          vz[l]=static_cast<real_prec>(0.5)*(v3[ij]+v3[l]);
        }
        if (i==0 && j>0 && k==0)
        {
          vx[l]=static_cast<real_prec>(0.5)*(v1[ik]+v1[l]);
          vy[l]=static_cast<real_prec>(0.5)*(v2[ik]+v2[l]);
          vz[l]=static_cast<real_prec>(0.5)*(v3[ik]+v3[l]);
        }
        if (i==0 && j==0 && k==0)
        {
          vx[l]=static_cast<real_prec>(0.5)*(v1[ijk]+v1[l]);
          vy[l]=static_cast<real_prec>(0.5)*(v2[ijk]+v2[l]);
          vz[l]=static_cast<real_prec>(0.5)*(v3[ijk]+v3[l]);
        }
        if (i>0 && j==0 && k==0)
        {
          vx[l]=static_cast<real_prec>(0.5)*(v1[jk]+v1[l]);
          vy[l]=static_cast<real_prec>(0.5)*(v2[jk]+v2[l]);
          vz[l]=static_cast<real_prec>(0.5)*(v3[jk]+v3[l]);
        }
        if (i>0 && j==0 && k>0)
        {
          vx[l]=static_cast<real_prec>(0.5)*(v1[jj]+v1[l]);
          vy[l]=static_cast<real_prec>(0.5)*(v2[jj]+v2[l]);
          vz[l]=static_cast<real_prec>(0.5)*(v3[jj]+v3[l]);
        }
        if (i>0 && j>0 && k==0)
        {
          vx[l]=static_cast<real_prec>(0.5)*(v1[kk]+v1[l]);
          vy[l]=static_cast<real_prec>(0.5)*(v2[kk]+v2[l]);
          vz[l]=static_cast<real_prec>(0.5)*(v3[kk]+v3[l]);
        }
      }

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (ULONG l=0;l<N; l++) // EGP: long -> ULONG
  {
    v1[l]=vx[l];
    v2[l]=vy[l];
    v3[l]=vz[l];
  }
  /* end: interpolate from cell center to cell boundaries */
}

void cellboundcomp(unsigned int N1, unsigned int N2, unsigned int N3, real_prec *vi)
{
  ULONG  N=N1*N2*N3;

  /* start: interpolate from cell center to cell boundaries */
  fftw_array<real_prec> viout(N);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (unsigned i = 0; i < N1; i++)
    for (unsigned j = 0; j < N2; j++)
      for (unsigned k = 0; k < N3; k++)
      {
        ULONG l=k+N3*(j+N2*i);
        ULONG m=k-1+N3*(j-1+N2*(i-1));

        if (i>0 && j>0 && k>0)
        {
          viout[l]=static_cast<real_prec>(0.5)*(vi[m]+vi[l]);
        }

        /* periodic boundary conditions */

        ULONG ii=k-1+N3*(j-1+N2*(N1-1));
        ULONG jj=k-1+N3*(N2-1+N2*(i-1));
        ULONG kk=N3-1+N3*(j-1+N2*(i-1));
        ULONG ij=k-1+N3*(N2-1+N2*(N1-1));
        ULONG ik=N3-1+N3*(j-1+N2*(N1-1));
        ULONG jk=N3-1+N3*(N2-1+N2*(i-1));
        ULONG ijk=N3-1+N3*(N2-1+N2*(N1-1));

        if (i==0 && j>0 && k>0)
        {
          viout[l]=static_cast<real_prec>(0.5)*(vi[ii]+vi[l]);
        }
        if (i==0 && j==0 && k>0)
        {
          viout[l]=static_cast<real_prec>(0.5)*(vi[ij]+vi[l]);
        }
        if (i==0 && j>0 && k==0)
        {
          viout[l]=static_cast<real_prec>(0.5)*(vi[ik]+vi[l]);
        }
        if (i==0 && j==0 && k==0)
        {
          viout[l]=static_cast<real_prec>(0.5)*(vi[ijk]+vi[l]);
        }
        if (i>0 && j==0 && k==0)
        {
          viout[l]=static_cast<real_prec>(0.5)*(vi[jk]+vi[l]);
        }
        if (i>0 && j==0 && k>0)
        {
          viout[l]=static_cast<real_prec>(0.5)*(vi[jj]+vi[l]);
        }
        if (i>0 && j>0 && k==0)
        {
          viout[l]=static_cast<real_prec>(0.5)*(vi[kk]+vi[l]);
        }
      }

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (ULONG l = 0; l < N; l++)
  {
    vi[l]=viout[l];
  }
  /* end: interpolate from cell center to cell boundaries */
}

