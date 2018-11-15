/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include "define_opt.h"
#include "fftw_array.h"
#include "fftwrapper.h" // fftC2R
#include "convenience.h" // copyArray
#include "random.hpp" // resolution independent random grid


////////////////////////////////////////////////////////////////
// Random functions
////////////////////////////////////////////////////////////////

real_prec GR_NUM(gsl_rng * SEED, real_prec sigma ,int GR_METHOD)
{
  real_prec val=0.;
  switch (GR_METHOD)
  {
    case 0:
    {
      val=static_cast<real_prec>(gsl_ran_gaussian(SEED, sigma));
    }
      break;

    case 1:
    {
      val=static_cast<real_prec>(gsl_ran_gaussian_ziggurat(SEED, sigma));
    }
      break;

    case 2:
    {
      val=static_cast<real_prec>(gsl_ran_gaussian_ratio_method(SEED, sigma));
    }
      break;

  }
  return(val);
}


void create_GARFIELD(unsigned N1, unsigned N2, unsigned N3, real_prec L1, real_prec L2, real_prec L3,real_prec *delta,
                     const real_prec * Power, gsl_rng * seed)
{
  /// FOLLOWING HUGO MARTEL
  //
  // New version, EGP: uses a "resolution independent random grid",
  // which makes sure that, when given the same seed, the random grid
  // that is generated for gridsize N is also generated for gridsize
  // 2*N, plus the higher k terms (smaller scale structure).
  // Even though there are no more random calls in the for loop, there
  // is the possibility that one array location is accessed simultaneously
  // by multiple threads, so no multithreading here! Unless you want to
  // put atomics around all the calls...


  ULONG N=N1*N2*N3;
  ULONG  Nhalf=N1*N2*(N3/2+1);

  real_prec Vol = L1*L2*L3;

  fftw_array<complex_prec> GRF_array(N),GRF_array2(Nhalf);

  // Build a resolution independent random grid
  std::vector< std::complex<real_prec> > random_grid = resolution_independent_random_grid_FS<real_prec>( N1, seed, false );
  auto *random_grid_array = reinterpret_cast<complex_prec *>(random_grid.data());
  copyArray(random_grid_array, GRF_array, N);

// Factor for conversion from power spectrum to discrete fourier transform
// amplitude (see Martel 2005):
  real_prec ps2dft_amp = 1.;
#ifdef FOURIER_DEF_1
  ps2dft_amp = 1./Vol;
#endif // FOURIER_DEF_1
#ifdef FOURIER_DEF_2
  ps2dft_amp = static_cast<real_prec>(N*N)/Vol;
#endif // FOURIER_DEF_2

  // TODO:
  // We can shorten the code considerably by just multiplying the whole array
  // by sigma first in a separate (multithreaded) for-loop (regardless of all
  // the special cases in the complicated loop below) and only after that
  // doing all the special cases below, but we can then leave out all the
  // *= sigma parts. Actually, I think what is done here is just two things:
  // 1. making the array hermitian, and 2. a couple of special cases that need
  // to be non-imaginary. The former we can do in a separate function (which is
  // already in the code, I believe) and the latter we can just do with a few
  // direct assignments, instead of having to go through all the if's at every
  // step of the loop.

  //EGP  real_prec kr=0.;
  real_prec sigma=0.;
//#ifdef MULTITHREAD
//#pragma omp parallel for
//#endif // MULTITHREAD
  for (unsigned i=0 ; i<=N1/2;i++)
    for (unsigned j=0 ; j<=N2/2;j++)
      for (unsigned k=0 ; k<=N3/2;k++)
      {
        sigma = sqrt(ps2dft_amp * Power[k+N3*(j+N2*i)]/num_2);

        if( (i>0 && i<N1/2) && (j>0 && j<N2/2) && (k>0 && k<N3/2))
        {
          unsigned ii=N1-i; unsigned jj=N2-j; unsigned kk=N3-k;

          re(GRF_array[k+N3*(j+N2*i)]) *= sigma;
          im(GRF_array[k+N3*(j+N2*i)]) *= sigma;

          re(GRF_array[kk+N3*(jj+N2*ii)])=re(GRF_array[k+N3*(j+N2*i)]);
          im(GRF_array[kk+N3*(jj+N2*ii)])=-im(GRF_array[k+N3*(j+N2*i)]);

          ///**********************************************************

          re(GRF_array[k+N3*(j+N2*ii)]) *= sigma;
          im(GRF_array[k+N3*(j+N2*ii)]) *= sigma;

          re(GRF_array[kk+N3*(jj+N2*i)])=re(GRF_array[k+N3*(j+N2*ii)]);
          im(GRF_array[kk+N3*(jj+N2*i)])=-im(GRF_array[k+N3*(j+N2*ii)]);

          ///**********************************************************

          re(GRF_array[k+N3*(jj+N2*i)]) *= sigma;
          im(GRF_array[k+N3*(jj+N2*i)]) *= sigma;

          re(GRF_array[kk+N3*(j+N2*ii)])=re(GRF_array[k+N3*(jj+N2*i)]);
          im(GRF_array[kk+N3*(j+N2*ii)])=-im(GRF_array[k+N3*(jj+N2*i)]);

          ///**********************************************************

          re(GRF_array[kk+N3*(j+N2*i)]) *= sigma;
          im(GRF_array[kk+N3*(j+N2*i)]) *= sigma;

          re(GRF_array[k+N3*(jj+N2*ii)])=re(GRF_array[kk+N3*(j+N2*i)]);
          im(GRF_array[k+N3*(jj+N2*ii)])=-im(GRF_array[kk+N3*(j+N2*i)]);

        }

        if( (i>0 && i<N1/2) && (j>0 && j<N2/2) && (k==N3/2))
        {
          unsigned ii=N1-i; unsigned jj=N2-j;

          re(GRF_array[N3/2+N3*(j+N2*i)]) *= sigma;
          im(GRF_array[N3/2+N3*(j+N2*i)]) *= sigma;

          re(GRF_array[N3/2+N3*(jj+N2*ii)])=re(GRF_array[N3/2+N3*(j+N2*i)]);
          im(GRF_array[N3/2+N3*(jj+N2*ii)])=-im(GRF_array[N3/2+N3*(j+N2*i)]);

          ///**********************************************************

          re(GRF_array[N3/2+N3*(j+N2*ii)]) *= sigma;
          im(GRF_array[N3/2+N3*(j+N2*ii)]) *= sigma;

          re(GRF_array[N3/2+N3*(jj+N2*i)])=re(GRF_array[N3/2+N3*(j+N2*ii)]);
          im(GRF_array[N3/2+N3*(jj+N2*i)])=-im(GRF_array[N3/2+N3*(j+N2*ii)]);

        }

        if( (i>0 && i<N1/2) && (j==N2/2) && (k>0 && k<N3/2))
        {
          unsigned ii=N1-i; unsigned kk=N3-k;

          re(GRF_array[k+N3*(N2/2+N2*i)]) *= sigma;
          im(GRF_array[k+N3*(N2/2+N2*i)]) *= sigma;

          re(GRF_array[kk+N3*(N2/2+N2*ii)])=re(GRF_array[k+N3*(N2/2+N2*i)]);
          im(GRF_array[kk+N3*(N2/2+N2*ii)])=-im(GRF_array[k+N3*(N2/2+N2*i)]);

          ///**********************************************************

          re(GRF_array[k+N3*(N2/2+N2*ii)]) *= sigma;
          im(GRF_array[k+N3*(N2/2+N2*ii)]) *= sigma;

          re(GRF_array[kk+N3*(N2/2+N2*i)])=re(GRF_array[k+N3*(N2/2+N2*ii)]);
          im(GRF_array[kk+N3*(N2/2+N2*i)])=-im(GRF_array[k+N3*(N2/2+N2*ii)]);

        }

        if( (i==N1/2) && (j>0 && j<N2/2) && (k>0 && k<N3/2))
        {
          unsigned jj=N2-j; unsigned kk=N3-k;

          re(GRF_array[k+N3*(j+N2*N1/2)]) *= sigma;
          im(GRF_array[k+N3*(j+N2*N1/2)]) *= sigma;

          re(GRF_array[kk+N3*(jj+N2*N1/2)])=re(GRF_array[k+N3*(j+N2*N1/2)]);
          im(GRF_array[kk+N3*(jj+N2*N1/2)])=-im(GRF_array[k+N3*(j+N2*N1/2)]);

          ///**********************************************************

          re(GRF_array[k+N3*(jj+N2*N1/2)]) *= sigma;
          im(GRF_array[k+N3*(jj+N2*N1/2)]) *= sigma;

          re(GRF_array[kk+N3*(j+N2*N1/2)])=re(GRF_array[k+N3*(jj+N2*N1/2)]);
          im(GRF_array[kk+N3*(j+N2*N1/2)])=-im(GRF_array[k+N3*(jj+N2*N1/2)]);

        }

        if( (i==N1/2) && (j==N2/2) && (k>0 && k<N3/2))
        {
          unsigned kk=N3-k;

          re(GRF_array[k+N3*(N2/2+N2*N1/2)]) *= sigma;
          im(GRF_array[k+N3*(N2/2+N2*N1/2)]) *= sigma;

          re(GRF_array[kk+N3*(N2/2+N2*N1/2)])=re(GRF_array[k+N3*(N2/2+N2*N1/2)]);
          im(GRF_array[kk+N3*(N2/2+N2*N1/2)])=-im(GRF_array[k+N3*(N2/2+N2*N1/2)]);

          ///**********************************************************
        }

        if( (i>0 && i<N1/2) && (j==N2/2) && (k==N3/2))
        {
          unsigned ii=N1-i;

          re(GRF_array[N3/2+N3*(N2/2+N2*i)]) *= sigma;
          im(GRF_array[N3/2+N3*(N2/2+N2*i)]) *= sigma;

          re(GRF_array[N3/2+N3*(N2/2+N2*ii)])=re(GRF_array[N3/2+N3*(N2/2+N2*i)]);
          im(GRF_array[N3/2+N3*(N2/2+N2*ii)])=-im(GRF_array[N3/2+N3*(N2/2+N2*i)]);

          ///**********************************************************
        }

        if( (i==N1/2) && (j>0 && j<N2/2) && (k==N3/2))
        {
          unsigned jj=N2-j;

          re(GRF_array[N3/2+N3*(j+N2*N1/2)]) *= sigma;
          im(GRF_array[N3/2+N3*(j+N2*N1/2)]) *= sigma;

          re(GRF_array[N3/2+N3*(jj+N2*N1/2)])=re(GRF_array[N3/2+N3*(j+N2*N1/2)]);
          im(GRF_array[N3/2+N3*(jj+N2*N1/2)])=-im(GRF_array[N3/2+N3*(j+N2*N1/2)]);

          ///**********************************************************
        }

        if( (i==N1/2) && (j==N2/2) && (k==N3/2))
        {
          re(GRF_array[N3/2+N3*(N2/2+N2*N1/2)]) *= sqrt(num_2)*sigma;
          im(GRF_array[N3/2+N3*(N2/2+N2*N1/2)])=0.;

          ///**********************************************************
        }

        if( (i>0 && i<N1/2) && (j>0 && j<N2/2) && k==0)
        {
          unsigned ii=N1-i; unsigned jj=N2-j;

          re(GRF_array[0+N3*(j+N2*i)]) *= sigma;
          im(GRF_array[0+N3*(j+N2*i)]) *= sigma;

          re(GRF_array[0+N3*(jj+N2*ii)])=re(GRF_array[k+N3*(j+N2*i)]);
          im(GRF_array[0+N3*(jj+N2*ii)])=-im(GRF_array[k+N3*(j+N2*i)]);

          ///**********************************************************

          re(GRF_array[0+N3*(j+N2*ii)]) *= sigma;
          im(GRF_array[0+N3*(j+N2*ii)]) *= sigma;

          re(GRF_array[0+N3*(jj+N2*i)])=re(GRF_array[0+N3*(j+N2*ii)]);
          im(GRF_array[0+N3*(jj+N2*i)])=-im(GRF_array[0+N3*(j+N2*ii)]);
        }

        if( (i>0 && i<N1/2) && (j==0) && (k>0 && k<N3/2))
        {
          unsigned ii=N1-i; unsigned kk=N3-k;

          re(GRF_array[k+N3*(0+N2*i)]) *= sigma;
          im(GRF_array[k+N3*(0+N2*i)]) *= sigma;

          re(GRF_array[kk+N3*(0+N2*ii)])=re(GRF_array[k+N3*(0+N2*i)]);
          im(GRF_array[kk+N3*(0+N2*ii)])=-im(GRF_array[k+N3*(0+N2*i)]);

          ///**********************************************************

          re(GRF_array[k+N3*(0+N2*ii)]) *= sigma;
          im(GRF_array[k+N3*(0+N2*ii)]) *= sigma;

          re(GRF_array[kk+N3*(0+N2*i)])=re(GRF_array[k+N3*(0+N2*ii)]);
          im(GRF_array[kk+N3*(0+N2*i)])=-im(GRF_array[k+N3*(0+N2*ii)]);
        }

        if( (i==0) && (j>0 && j<N2/2) && (k>0 && k<N3/2))
        {
          unsigned jj=N2-j; unsigned kk=N3-k;

          re(GRF_array[k+N3*(j+N2*0)]) *= sigma;
          im(GRF_array[k+N3*(j+N2*0)]) *= sigma;

          re(GRF_array[kk+N3*(jj+N2*0)])=re(GRF_array[k+N3*(j+N2*0)]);
          im(GRF_array[kk+N3*(jj+N2*0)])=-im(GRF_array[k+N3*(j+N2*0)]);

          ///**********************************************************

          re(GRF_array[k+N3*(jj+N2*0)]) *= sigma;
          im(GRF_array[k+N3*(jj+N2*0)]) *= sigma;

          re(GRF_array[kk+N3*(j+N2*0)])=re(GRF_array[k+N3*(jj+N2*0)]);
          im(GRF_array[kk+N3*(j+N2*0)])=-im(GRF_array[k+N3*(jj+N2*0)]);
        }

        if( (i>0 && i<N1/2) && (j==0) && (k==0))
        {
          unsigned ii=N1-i;

          re(GRF_array[0+N3*(0+N2*i)]) *= sigma;
          im(GRF_array[0+N3*(0+N2*i)]) *= sigma;

          re(GRF_array[0+N3*(0+N2*ii)])=re(GRF_array[0+N3*(0+N2*i)]);
          im(GRF_array[0+N3*(0+N2*ii)])=-im(GRF_array[0+N3*(0+N2*i)]);
          ///**********************************************************

        }

        if( (i==0) && (j==0) && (k>0 && k<N3/2))
        {
          unsigned kk=N3-k;

          re(GRF_array[k+N3*(0+N2*0)]) *= sigma;
          im(GRF_array[k+N3*(0+N2*0)]) *= sigma;

          re(GRF_array[kk+N3*(0+N2*0)])=re(GRF_array[k+N3*(0+N2*0)]);
          im(GRF_array[kk+N3*(0+N2*0)])=-im(GRF_array[k+N3*(0+N2*0)]);
          ///**********************************************************

        }

        if( (i==0) && (j>0 && j<N2/2) && (k==0))
        {
          unsigned jj=N2-j;

          re(GRF_array[0+N3*(j+N2*0)]) *= sigma;
          im(GRF_array[0+N3*(j+N2*0)]) *= sigma;

          re(GRF_array[0+N3*(jj+N2*0)])=re(GRF_array[0+N3*(j+N2*0)]);
          im(GRF_array[0+N3*(jj+N2*0)])=-im(GRF_array[0+N3*(j+N2*0)]);
          ///**********************************************************

        }

        if( (i==0) && (j==0) && (k==0))
        {
          re(GRF_array[0])=0.;
          im(GRF_array[0])=0.;
          ///**********************************************************

        }

        if( (i>0 && i<N1/2) && (j==N2/2) && (k==0))
        {
          unsigned ii=N1-i;

          re(GRF_array[0+N3*(N2/2+N2*i)]) *= sigma;
          im(GRF_array[0+N3*(N2/2+N2*i)]) *= sigma;

          re(GRF_array[0+N3*(N2/2+N2*ii)])=re(GRF_array[0+N3*(N2/2+N2*i)]);
          im(GRF_array[0+N3*(N2/2+N2*ii)])=-im(GRF_array[0+N3*(N2/2+N2*i)]);

          ///**********************************************************

        }

        if( (i>0 && i<N1/2) && (j==0) && (k==N3/2))
        {
          unsigned ii=N1-i;

          re(GRF_array[N3/2+N3*(0+N2*i)]) *= sigma;
          im(GRF_array[N3/2+N3*(0+N2*i)]) *= sigma;

          re(GRF_array[N3/2+N3*(0+N2*ii)])=re(GRF_array[N3/2+N3*(0+N2*i)]);
          im(GRF_array[N3/2+N3*(0+N2*ii)])=-im(GRF_array[N3/2+N3*(0+N2*i)]);

          ///**********************************************************

        }

        if( (i==N1/2) && (j>0 && j<N2/2) && (k==0))
        {
          unsigned jj=N2-j;

          re(GRF_array[0+N3*(j+N2*N1/2)]) *= sigma;
          im(GRF_array[0+N3*(j+N2*N1/2)]) *= sigma;

          re(GRF_array[0+N3*(jj+N2*N1/2)])=re(GRF_array[0+N3*(j+N2*N1/2)]);
          im(GRF_array[0+N3*(jj+N2*N1/2)])=-im(GRF_array[0+N3*(j+N2*N1/2)]);

          ///**********************************************************

        }


        if( (i==0) && (j>0 && j<N2/2) && (k==N3/2))
        {
          unsigned jj=N2-j;

          re(GRF_array[N3/2+N3*(j+N2*0)]) *= sigma;
          im(GRF_array[N3/2+N3*(j+N2*0)]) *= sigma;

          re(GRF_array[N3/2+N3*(jj+N2*0)])=re(GRF_array[N3/2+N3*(j+N2*0)]);
          im(GRF_array[N3/2+N3*(jj+N2*0)])=-im(GRF_array[N3/2+N3*(j+N2*0)]);

          ///**********************************************************
        }

        if( (i==N1/2) && (j==0) && (k>0 && k<N3/2))
        {
          unsigned kk=N3-k;

          re(GRF_array[k+N3*(0+N2*N1/2)]) *= sigma;
          im(GRF_array[k+N3*(0+N2*N1/2)]) *= sigma;

          re(GRF_array[kk+N3*(0+N2*N1/2)])=re(GRF_array[k+N3*(0+N2*N1/2)]);
          im(GRF_array[kk+N3*(0+N2*N1/2)])=-im(GRF_array[k+N3*(0+N2*N1/2)]);

          ///**********************************************************

        }

        if( (i==0) && (j==N2/2) && (k>0 && k<N3/2))
        {
          unsigned kk=N3-k;

          re(GRF_array[k+N3*(N2/2+N2*0)]) *= sigma;
          im(GRF_array[k+N3*(N2/2+N2*0)]) *= sigma;

          re(GRF_array[kk+N3*(N2/2+N2*0)])=re(GRF_array[k+N3*(N2/2+N2*0)]);
          im(GRF_array[kk+N3*(N2/2+N2*0)])=-im(GRF_array[k+N3*(N2/2+N2*0)]);

          ///**********************************************************

        }

        if( (i==0) && (j==0) && (k==N3/2))
        {
          re(GRF_array[N3/2+N3*(0+N2*0)]) *= sqrt(num_2)*sigma;
          im(GRF_array[N3/2+N3*(0+N2*0)])=0.;
          ///**********************************************************

        }

        if( (i==0) && (j==N2/2) && (k==0))
        {
          re(GRF_array[0+N3*(N2/2+N2*0)]) *= sqrt(num_2)*sigma;
          im(GRF_array[0+N3*(N2/2+N2*0)])=0.;
          ///**********************************************************

        }

        if( (i==N1/2) && (j==0) && (k==0))
        {
          re(GRF_array[0+N3*(0+N2*N1/2)]) *= sqrt(num_2)*sigma;
          im(GRF_array[0+N3*(0+N2*N1/2)])=0.;
          ///**********************************************************

        }

        if( (i==N1/2) && (j==N2/2) && (k==0))
        {
          re(GRF_array[0+N3*(N2/2+N2*N1/2)]) *= sqrt(num_2)*sigma;
          im(GRF_array[0+N3*(N2/2+N2*N1/2)])=0.;
          ///**********************************************************

        }

        if( (i==N1/2) && (j==0) && (k==N3/2))
        {
          re(GRF_array[N3/2+N3*(0+N2*N1/2)]) *= sqrt(num_2)*sigma;
          im(GRF_array[N3/2+N3*(0+N2*N1/2)])=0.;
          ///**********************************************************

        }

        if( (i==0) && (j==N2/2) && (k==N3/2))
        {
          re(GRF_array[N3/2+N3*(N2/2+N2*0)]) *= sqrt(num_2)*sigma;
          im(GRF_array[N3/2+N3*(N2/2+N2*0)])=0.;
          ///**********************************************************
        }

      }

  /*
     FFT3d (N1,N2,N3, to_Rspace, GRF_array, GRF_array);

     for (ULONG i=0 ; i<N;i++)
     delta[i] = re(GRF_array[i]);
     */

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (unsigned i=0 ; i<N1;i++)
    for (unsigned j=0 ; j<N2;j++)
      for (unsigned k=0 ; k<=N3/2;k++)
      {
        ULONG ihalf=k+(N3/2+1)*(j+N2*i);
        ULONG iind=k+N3*(j+N2*i);

        re(GRF_array2[ihalf])=re(GRF_array[iind]);
        im(GRF_array2[ihalf])=im(GRF_array[iind]);
      }

  fftC2R(N1,N2,N3,GRF_array2,delta);

}
