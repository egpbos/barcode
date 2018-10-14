/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include "struct_main.h"
#include "fftw_array.h"

#include <math.h>
#include <iomanip>

#include "cosmo.h"
#include "IOfunctionsGen.h"
#include "math_funcs.h"

#include "convenience.h"
#include "BarcodeException.h"

// Everything in here is only used in Lag2Eul

using namespace std;

void PoissonSolver(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec L2, real_prec L3,
                   real_prec *delta, real_prec *Pot)
{
  // ULONG N=N1*N2*N3;
  unsigned int N3half = (N3/2+1);
  ULONG Nhalf = N1*N2*N3half;

  // fftw_array<complex_prec> deltaC(N), PotC(N), AUX(N);
  fftw_array<complex_prec> deltaC(Nhalf);

  // complexify_array(delta, deltaC, N);
  // FFT3d(N1,N2,N3, true, deltaC, AUX);
  fftR2C(N1, N2, N3, delta, deltaC);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (unsigned i=0; i<N1;i++)
    for (unsigned j=0; j<N2;j++)
      for (unsigned k = 0; k < N3half; ++k) {
        ULONG index=k + N3half * (j + N2*i);

        real_prec kmod2=k_squared(i,j,k,L1,L2,L3,N1,N2,N3);

        real_prec fackern=0.;
        if (kmod2>0.)
          fackern=static_cast<real_prec>(-1./kmod2);

        re(deltaC[index]) = re(deltaC[index]) * fackern;  // deltaC -> "PotC"
        im(deltaC[index]) = im(deltaC[index]) * fackern;
      }

  // FFT3d(N1,N2,N3, false, AUX, PotC);
  // real_part_array(PotC, Pot, N);
  fftC2R(N1, N2, N3, deltaC, Pot);
}


// Also used in webclass, but effectively Lag2Eul only.
void calc_LapPhiv(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, complex_prec *philv, real_prec *LapPhiv,
                  int index1, int index2)
{
  ULONG N=N1*N2*N3;
  // int N3half = (N3/2 + 1);
  // ULONG Nhalf = N1*N2*N3half;
  fftw_array<complex_prec> LapPhivl(N), LapPhivC(N);
  // fftw_array<complex_prec> LapPhivl(Nhalf);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD 
  for (unsigned i=0;i<N1;i++)
    for (unsigned j=0;j<N2;j++)
      for (unsigned k = 0; k < N3; ++k) {
        ULONG ii = k + N3 * (j + N2*i);
        real_prec k1 = 0., k2 = 0.;

        real_prec kx=calc_kx(i,L1,N1);
        real_prec ky=calc_ky(j,L1,N2);
        real_prec kz=calc_kz(k,L1,N3);

        switch (index1)
        {
          case 1:
            k1=kx;
            break;
          case 2:
            k1=ky;
            break;
          case 3:
            k1=kz;
            break;
          default:
            throw BarcodeException("In calc_LapPhiv : index1 must be either 1, 2 or 3!");
        }

        switch (index2)
        {
          case 1:
            k2=kx;
            break;
          case 2:
            k2=ky;
            break;
          case 3:
            k2=kz;
            break;
          default:
            throw BarcodeException("In calc_LapPhiv : index2 must be either 1, 2 or 3!");
        }

        re(LapPhivl[ii])=-k1*k2*re(philv[ii]); // EGP (20131213): note that this was wrong! philv[ihalf] was called,
        im(LapPhivl[ii])=-k1*k2*im(philv[ii]); // but should be ii, because philv was already a full N size array!
      }

  FFT3d(N1,N2,N3, false, LapPhivl, LapPhivC);
  real_part_array(LapPhivC, LapPhiv, N);
  // fftC2R(N1, N2, N3, LapPhivl, LapPhiv);
}


real_prec linearvel3d(int index, real_prec kx, real_prec ky, real_prec kz, real_prec phi)
{
  real_prec out, kl;

  switch(index)
  {
    case 1:
      {
        kl=kx;
        break;
      }
    case 2:
      {
        kl=ky;
        break;
      }
    case 3:
      {
        kl=kz;
        break;
      }
    default:
      {
        throw BarcodeException("In linearvel3d : index must be either 1, 2 or 3!");
      }
  }
  real_prec kmod2=kx*kx+ky*ky+kz*kz;

  real_prec fackern=0.0;
  if (kmod2>eps)
    fackern = kl/kmod2;

  out=fackern*phi;

  return out;
}


void theta2vel(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec L2, real_prec L3,
               real_prec scale, real_prec Omega_M,
               real_prec Omega_L, real_prec *delta, real_prec *vex, real_prec *vey, real_prec *vez, bool zeropad,
               bool norm, plan_pkg *R2Cplan, plan_pkg *C2Rplan)
{
  // The fastest way to call this function is to use:
  // 1. R2Cplan->R for delta
  // 2. C2Rplan->R for vez
  real_prec Lzp1=L1, Lzp2=L2, Lzp3=L3;
  unsigned Nzp1=N1, Nzp2=N2, Nzp3=N3;

  if (zeropad==true)
  {
    Lzp1=L1*static_cast<real_prec>(2.0);
    Lzp2=L2*static_cast<real_prec>(2.0);
    Lzp3=L3*static_cast<real_prec>(2.0);

    Nzp1=N1*2;
    Nzp2=N2*2;
    Nzp3=N3*2;
  }

  //real_prec H0=static_cast<real_prec>(100.*hconst *cgs_km/cgs_Mpc/cgs_sec);
  real_prec cpecvel=1.;

  if (norm==true)
    //cpecvel=c_pecvel(scale,Omega_M,Omega_L,H0);
    cpecvel = c_pecvel(scale, Omega_M, Omega_L, 1); // Note: only first order (Zel'dovich) f!

  unsigned N3half = (N3/2 + 1);
  ULONG Nhalf = N1 * N2 * N3half;
  // fftw_array<complex_prec> deltaC(Nhalf), velx(Nhalf), vely(Nhalf);
  fftw_array<complex_prec> velx(Nhalf), vely(Nhalf);

  // fftR2C(N1, N2, N3, delta, deltaC);
  fftR2Cplanned(delta, R2Cplan->C, R2Cplan);  // using R2Cplan->C as deltaC

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif // MULTITHREAD
  for (unsigned i=0;i<N1;i++) {
    real_prec kx=calc_kx(i,Lzp1,Nzp1);
    for (unsigned j=0;j<N2;j++) {
      real_prec ky=calc_ky(j,Lzp2,Nzp2);
      for (unsigned k = 0; k < N3half; ++k) {
        real_prec kz=calc_kz(k,Lzp3,Nzp3);
        real_prec ksq = kx*kx+ky*ky+kz*kz;

        if (ksq > eps) {  // eps is defined in define_opt.h as 1.e-14
          real_prec fac_kern = cpecvel / ksq;
          ULONG ix = k + N3half * (j + N2*i);
          // real_prec deltaC_ix_r = re(deltaC[ix]);
          // real_prec deltaC_ix_i = im(deltaC[ix]);
          real_prec deltaC_ix_r = re(R2Cplan->C[ix]);
          real_prec deltaC_ix_i = im(R2Cplan->C[ix]);

          real_prec fac_kern_x = fac_kern * kx;
          re(velx[ix]) = fac_kern_x * deltaC_ix_i;
          im(velx[ix]) = fac_kern_x * -deltaC_ix_r;

          real_prec fac_kern_y = fac_kern * ky;
          re(vely[ix]) = fac_kern_y * deltaC_ix_i;
          im(vely[ix]) = fac_kern_y * -deltaC_ix_r;

          // using deltaC as dummy instead of velz
          real_prec fac_kern_z = fac_kern * kz;
          // re(deltaC[ix]) = fac_kern_z * deltaC_ix_i;
          // im(deltaC[ix]) = fac_kern_z * -deltaC_ix_r;
          re(R2Cplan->C[ix]) = fac_kern_z * deltaC_ix_i;
          im(R2Cplan->C[ix]) = fac_kern_z * -deltaC_ix_r;
        } else {
          ULONG ix = k + N3half * (j + N2*i);
          re(velx[ix]) = 0;
          im(velx[ix]) = 0;
          re(vely[ix]) = 0;
          im(vely[ix]) = 0;
          // re(deltaC[ix]) = 0;
          // im(deltaC[ix]) = 0;
          re(R2Cplan->C[ix]) = 0;
          im(R2Cplan->C[ix]) = 0;
        }

        // From http://math.mit.edu/~stevenj/fft-deriv.pdf: Nyquist components
        // should be zero, in every dimension, otherwise the result is complex.
        // This function is not actually treated in that document, but I would
        // say it also counts as an odd-ordered derivative.
        if ((i==N1/2) || (j==N2/2) || (k==N3/2))
        {
          ULONG ix = k + N3half * (j + N2*i);
          re(velx[ix]) = 0.;
          im(velx[ix]) = 0.;
          re(vely[ix]) = 0.;
          im(vely[ix]) = 0.;
          // re(deltaC[ix]) = 0;
          // im(deltaC[ix]) = 0;
          re(R2Cplan->C[ix]) = 0;
          im(R2Cplan->C[ix]) = 0;
        }
      }
    }
  }

  // fftC2R(N1, N2, N3, velx, vex);
  // fftC2R(N1, N2, N3, vely, vey);
  // fftC2R(N1, N2, N3, deltaC, vez);
  // fftC2R(N1, N2, N3, R2Cplan->C, vez);
  fftC2Rplanned(R2Cplan->C, vez, C2Rplan);  // must go first, C2R destroys input!
  fftC2Rplanned(velx, vex, C2Rplan);
  fftC2Rplanned(vely, vey, C2Rplan);
}


void theta2velcomp(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec L2, real_prec L3,
                   real_prec scale, real_prec Omega_M, real_prec Omega_L, real_prec *delta, real_prec *vei,
                   bool zeropad, bool norm, int comp)
{
  real_prec Lzp1=L1, Lzp2=L2, Lzp3=L3;
  unsigned Nzp1=N1, Nzp2=N2, Nzp3=N3;

  if (zeropad==true)
  {
    Lzp1=L1*static_cast<real_prec>(2.0);
    Lzp2=L2*static_cast<real_prec>(2.0);
    Lzp3=L3*static_cast<real_prec>(2.0);

    Nzp1=N1*2;
    Nzp2=N2*2;
    Nzp3=N3*2;
  }

  //real_prec H0=static_cast<real_prec>(100.*hconst *cgs_km/cgs_Mpc/cgs_sec);
  real_prec cpecvel=1.;

  if (norm==true)
    //cpecvel=c_pecvel(scale,Omega_M,Omega_L,H0);
    cpecvel = c_pecvel(scale, Omega_M, Omega_L, 1); // Note: only first order (Zel'dovich) f!

  // ULONG N=N1*N2*N3;
  unsigned N3half = (N3/2 + 1);
  ULONG Nhalf = N1 * N2 * N3half;
  // fftw_array<complex_prec> dummyC(N), phi(N), veli(N);
  fftw_array<complex_prec> deltaC(Nhalf);

  // complexify_array(delta, dummyC, N);
  // FFT3d(N1,N2,N3, true, dummyC, phi);
  fftR2C(N1, N2, N3, delta, deltaC);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (unsigned i=0;i<N1;i++)
    for (unsigned j=0;j<N2;j++)
      for (unsigned k = 0; k < N3half; ++k) {
        ULONG index = k + N3half * (j + N2*i);

        real_prec kx=calc_kx(i,Lzp1,Nzp1);
        real_prec ky=calc_ky(j,Lzp2,Nzp2);
        real_prec kz=calc_kz(k,Lzp3,Nzp3);

        real_prec deltaC_ix_r = re(deltaC[index]);
        real_prec deltaC_ix_i = im(deltaC[index]);

        // deltaC is now used as dummy instead of veli
        switch (comp)
        {
          case 1:
            {
              re(deltaC[index])=cpecvel*linearvel3d(1, kx, ky, kz, deltaC_ix_i);
              im(deltaC[index])=cpecvel*linearvel3d(1, kx, ky, kz, -deltaC_ix_r);
            }
            break;
          case 2:
            {
              re(deltaC[index])=cpecvel*linearvel3d(2, kx, ky, kz, deltaC_ix_i);
              im(deltaC[index])=cpecvel*linearvel3d(2, kx, ky, kz, -deltaC_ix_r);

            }
            break;
          case 3:
            {
              re(deltaC[index])=cpecvel*linearvel3d(3, kx, ky, kz, deltaC_ix_i);
              im(deltaC[index])=cpecvel*linearvel3d(3, kx, ky, kz, -deltaC_ix_r);
            }
            break;
        }

        // From http://math.mit.edu/~stevenj/fft-deriv.pdf: Nyquist components
        // should be zero, in every dimension, otherwise the result is complex.
        // This function is not actually treated in that document, but I would
        // say it also counts as an odd-ordered derivative.
        if ((i==N1/2) || (j==N2/2) || (k==N3/2))
        {
          re(deltaC[index]) = 0.;
          im(deltaC[index]) = 0.;
        }
      }

  // FFT3d(N1,N2,N3, false, veli, dummyC);
  // real_part_array(dummyC, vei, N);
  fftC2R(N1, N2, N3, deltaC, vei);
}


/* start LPT packages*/

void calc_m2v_mem(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec *phiv, real_prec *m2v)
{
  ULONG N=N1*N2*N3;

  fftw_array<real_prec> LapPhivx(N), LapPhivy(N), LapPhivz(N);
  fftw_array<real_prec> LapPhivxy(N), LapPhivxz(N), LapPhivyz(N);

#ifdef GFFT
  fftw_array<complex_prec> philv(N);
  for(ULONG i=0;i<N;i++)
  {
    re(philv[i])=0.0;
    im(philv[i])=0.0;
  }

  FFT3dR2C (N1,N2,N3,phiv,philv);

  calc_LapPhiv(N1,N2,N3,L1,philv,LapPhivx,1,1);
  calc_LapPhiv(N1,N2,N3,L1,philv,LapPhivy,2,2);
  calc_LapPhiv(N1,N2,N3,L1,philv,LapPhivz,3,3);

  calc_LapPhiv(N1,N2,N3,L1,philv,LapPhivxy,1,2);
  calc_LapPhiv(N1,N2,N3,L1,philv,LapPhivxz,1,3);
  calc_LapPhiv(N1,N2,N3,L1,philv,LapPhivyz,2,3);
#endif  // GFFT

#ifdef GFINDIFF
  fftw_array<real_prec> dummy(N);
  gradfindif(N1, L1, phiv, dummy, 1);
  gradfindif(N1, L1, dummy, LapPhivx, 1);
  gradfindif(N1, L1, dummy, LapPhivxy, 2);
  gradfindif(N1, L1, dummy, LapPhivxz, 3);

  gradfindif(N1, L1, phiv, dummy, 2);
  gradfindif(N1, L1, dummy, LapPhivy, 2);
  gradfindif(N1, L1, dummy, LapPhivyz, 3);

  gradfindif(N1, L1, phiv, dummy, 3);
  gradfindif(N1, L1, dummy, LapPhivz, 3);
#endif  // GFINDIFF


#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for(ULONG i=0;i<N;i++)
  {
    m2v[i]=LapPhivx[i]*LapPhivy[i]-LapPhivxy[i]*LapPhivxy[i]+LapPhivx[i]*LapPhivz[i]-LapPhivxz[i]*LapPhivxz[i]+LapPhivy[i]*LapPhivz[i]-LapPhivyz[i]*LapPhivyz[i];
  }
}


