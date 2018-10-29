/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include "define_opt.h"
#include "struct_main.h"
#include "fftw_array.h"

#include <cmath>
#include <iomanip>
#include <cassert>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>

#include "cosmo.h"
#include "IOfunctionsGen.h"
#include "scale_space.hpp" // k_squared

#include "convenience.h"

using namespace std;

void readtab(struct DATA *data)
{
  real_prec L1=data->numerical->L1, L2=data->numerical->L2, L3=data->numerical->L3;
  unsigned N1=data->numerical->N1, N2=data->numerical->N2, N3=data->numerical->N3;
  //ULONG N=N1*N2*N3;

  string fnamePS=data->numerical->fnamePS;

  ULONG nbinPS=count_lines(fnamePS);

  fftw_array<float> kPS(nbinPS);
  fftw_array<float> powPS(nbinPS);

  cout<<"---> reading from ... "<<fnamePS<<endl;

  ifstream inStreamPS;
  inStreamPS.open(fnamePS.data());
  assert(inStreamPS.is_open());
  for(ULONG i=0;i<nbinPS;i++)
  {
    inStreamPS >>kPS[i]>>powPS[i];
    if (i<10)
      cout<<kPS[i]<<" "<<powPS[i]<<endl;
  }
  inStreamPS.close();

  fftw_array<gsl_real> kPSg(nbinPS);
  fftw_array<gsl_real> powPSg(nbinPS);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for(ULONG i=0;i<nbinPS;i++)
  {
    kPSg[i]=static_cast<gsl_real>(kPS[i]);
    powPSg[i]=static_cast<gsl_real>(powPS[i]);
  }

  //EGP  real_prec kmax=sqrt(k_squared(N1/2,N2/2,N3/2,L1,L2,L3,N1,N2,N3));
  //EGP  real_prec dk=static_cast<real_prec>(kmax/real_prec(nbinPS));

  gsl_interp *interpolation = gsl_interp_alloc (gsl_interp_linear,nbinPS);
  gsl_interp_init(interpolation,kPSg,powPSg,nbinPS);
  gsl_interp_accel * accelerator = gsl_interp_accel_alloc();

  //EGP  real_prec Vol=L1*L2*L3;
  //EGP  real_prec CONT_NORM=Vol;
  //EGP  real_prec Norm=1.0;

  //EGP  #ifdef          FOURIER_DEF_1
  //EGP   Norm=static_cast<real_prec>(1./CONT_NORM);
  //EGP #endif

  //EGP #ifdef          FOURIER_DEF_2
  //EGP   Norm=static_cast<real_prec>(real_prec(N)*real_prec(N)/CONT_NORM);
  //EGP #endif

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for(unsigned i=0;i<N1;i++)
    for(unsigned j=0;j<N2;j++)
      for(unsigned k=0;k<N3;k++)
      {
        real_prec k2=k_squared(i,j,k,L1,L2,L3,N1,N2,N3);
        gsl_real ktot=sqrt(k2);

        ULONG iid=k+N3*(j+N2*i);

        data->observational->Power[iid]=static_cast<real_prec>(gsl_interp_eval(interpolation,kPSg,powPSg,ktot,accelerator));//EGP *Norm);

        if(i==0 && j==0 && k==0)
          data->observational->Power[iid]=0.;
      }
  gsl_interp_accel_free (accelerator);
  gsl_interp_free(interpolation);
}

