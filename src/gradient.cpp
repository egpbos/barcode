/*
 * Barcode
 * Copyright E.G.P. Bos
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <stdexcept>
#include <string>

#include "fftw_array.h"
#include "fftwrapper.h"
#include "scale_space.hpp"

#include "gradient.hpp"

////////////////////////////////////////////////////////////////
// Gradient and other 3D derivative functions
////////////////////////////////////////////////////////////////

void gradfft(unsigned N1, unsigned N2, unsigned N3, real_prec L1, real_prec L2, real_prec L3, real_prec *in,
             real_prec *out, unsigned int dim)
{
  // ULONG N=N1*N2*N3;
  unsigned N3half = N3/2 + 1;
  ULONG Nhalf = N1 * N2 * N3half;
  // fftw_array<complex_prec> AUX(N),AUX2(N);
  fftw_array<complex_prec> AUX(Nhalf);

  // complexify_array(in, AUX2, N);
  // FFT3d(N1,N2,N3, true, AUX2, AUX);
  fftR2C(N1, N2, N3, in, AUX);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (unsigned i=0;i<N1;i++)
    for (unsigned j=0;j<N2;j++)
      // for (unsigned k=0;k<N3;k++)
      for (unsigned k = 0; k < N3half; ++k) {
        real_prec kl;
        switch (dim)
        {
          case 1:
            kl=calc_kx(i,L1,N1);
            break;
          case 2:
            kl=calc_ky(j,L2,N2);
            break;
          case 3:
            kl=calc_kz(k,L3,N3);
            break;
          default:
            throw std::runtime_error("In gradfft: dim should be either 1, 2 or 3!");
            break;
        }

        ULONG ll = k + N3half * (j + N2*i);

        // EGP: switched minus sign between terms! Was wrong way round.
        // EGP: ... and switched it back; Fourier def. assumption was wrong!
        real_prec dummy = re(AUX[ll]);  // for in-place
        re(AUX[ll]) = -kl * im(AUX[ll]);
        im(AUX[ll]) =  kl * dummy;

        // From http://math.mit.edu/~stevenj/fft-deriv.pdf: Nyquist components
        // should be zero, in every dimension, otherwise the result is complex.
        if ((i==N1/2) || (j==N2/2) || (k==N3/2))
        {
          re(AUX[ll]) = 0.;
          im(AUX[ll]) = 0.;
        }
      }
  // FFT3d(N1,N2,N3, false, AUX2, AUX);
  // real_part_array(AUX, out, N);
  fftC2R(N1, N2, N3, AUX, out);
}


void gradfindif(unsigned N1, real_prec L1, const real_prec *in, real_prec *out, unsigned int dim)
{
  if (N1 > INT_MAX) {
    std::string msg = "Box dimensions must not be larger than INT_MAX in gradfindif!";
    throw std::logic_error(msg);
    // otherwise the signed integer indices below will overflow
  }
  auto fac=static_cast<real_prec>(N1/(2.*L1));

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for(unsigned x = 0; x < N1; x++)
    for(unsigned y = 0; y < N1; y++)
      for(unsigned z = 0; z < N1; z++)
      {
        unsigned xr, xrr;
        int xl, xll;
        unsigned yr, yrr;
        int yl, yll;
        unsigned zr, zrr;
        int zl, zll;

        xrr = xr = x;
        yrr = yr = y;
        zrr = zr = z;

        xll = xl = static_cast<int>(x);
        yll = yl = static_cast<int>(y);
        zll = zl = static_cast<int>(z);

        int *il, *ill;
        unsigned *ir, *irr, *ii;

        switch (dim) {
          case 1: {
            ii = &x; il = &xl; ill = &xll; ir = &xr; irr = &xrr;
            break;
          }
          case 2: {
            ii = &y; il = &yl; ill = &yll; ir = &yr; irr = &yrr;
            break;
          }
          case 3: {
            ii = &z; il = &zl; ill = &zll; ir = &zr; irr = &zrr;
            break;
          }
          default: {
            throw std::runtime_error("dim must be 1, 2 or 3 in gradfindif!");
          }
        }

        *ir = *ii + 1;
        *il = static_cast<int>(*ii) - 1;
        *irr = *ii + 2;
        *ill = static_cast<int>(*ii) - 2;
        if(*ir >= N1)
          *ir -= N1;
        if(*irr >= N1)
          *irr -= N1;
        if(*il < 0)
          *il += N1;
        if(*ill < 0)
          *ill += N1;

        ULONG ix = z+N1*(y+N1*x);
        ULONG l = static_cast<unsigned>(zl)+N1*(static_cast<unsigned>(yl)+N1*static_cast<unsigned>(xl));
        ULONG r = zr+N1*(yr+N1*xr);
        ULONG ll = static_cast<unsigned>(zll)+N1*(static_cast<unsigned>(yll)+N1*static_cast<unsigned>(xll));
        ULONG rr = zrr+N1*(yrr+N1*xrr);

        out[ix] = -static_cast<real_prec>(fac*((4.0 / 3)*(in[l] - in[r]) - (1.0 / 6) * (in[ll] - in[rr])));
      }
}

// EGP: added this simultaneous (a) one component fourier space gradient and (b) inverse laplacian calculator, in Fourier space:
void grad_inv_lap_FS(unsigned N1, unsigned N2, unsigned N3, real_prec L1, real_prec L2, real_prec L3, complex_prec *in,
                     complex_prec *out, unsigned int index, bool rfft)
{
  unsigned kz_max = N3;
  if (rfft)
    kz_max = N3/2 + 1;

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (unsigned i=0;i<N1;i++)
    for (unsigned j=0;j<N2;j++)
      for (unsigned k = 0; k < kz_max; ++k) {
        ULONG ii=k + kz_max * (j + N2*i);

        real_prec kx=calc_kx(i,L1,N1);
        real_prec ky=calc_ky(j,L2,N2);
        real_prec kz=calc_kz(k,L3,N3);

        real_prec kmod = kx*kx + ky*ky + kz*kz;
        real_prec fac_kmod = 0.;
        if (kmod > 0)
          fac_kmod = 1/kmod;

        real_prec ki_over_kmod=0.;
        switch (index)
        {
          case 1:
            ki_over_kmod=kx*fac_kmod;
            break;
          case 2:
            ki_over_kmod=ky*fac_kmod;
            break;
          case 3:
            ki_over_kmod=kz*fac_kmod;
            break;
          default:
            throw std::runtime_error("In grad_inv_lap_FS: index must be either 1, 2 or 3!");
            break;
        }
        real_prec dummy = re(in[ii]); // to make in-place possible (in == out)
        re(out[ii]) =  ki_over_kmod*im(in[ii]);
        im(out[ii]) = -ki_over_kmod*dummy;

        // From http://math.mit.edu/~stevenj/fft-deriv.pdf: Nyquist components
        // should be zero, in every dimension, otherwise the result is complex.
        // This function is not actually treated in that document, but I would
        // say it also counts as an odd-ordered derivative.
        if ((i==N1/2) || (j==N2/2) || (k==N3/2))
        {
          re(out[ii]) = 0.;
          im(out[ii]) = 0.;
        }
      }
}
