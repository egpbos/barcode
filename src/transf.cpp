/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <iostream>

#include "transf.h"
#include "fftw_array.h"
//#include "math_funcs.h"
#include "scale_space.hpp"
#include "fftwrapper.h"
#include <cassert>
#include <fstream>
#include <cmath>
#include "IOfunctionsGen.h"

// barcoderunner only:
void transflpt(unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec L2, real_prec L3,
               int filtertype, std::string dir)
{
  // EGP:
  // Watch out, there is still a (probably wrong) FOURIER_DEF normalisation factor below here!!!

  bool zeld=false;
  bool secordlpt=false;
  bool denslpt=false;

  switch (filtertype)
  {
    case 1:
      zeld=true;
      break;
    case 2:
      secordlpt=true;
      break;
    case 3:
      denslpt=true;
      break;
  }


  ULONG N=N1*N2*N3;

  fftw_array<real_prec> out(N);

  int bmax=100;
  char buffer1[bmax];
  sprintf(buffer1,"init_spec.dat");
  std::string inputFileName= dir + std::string(buffer1);
  std::ifstream inStream;
  inStream.open(inputFileName.data());
  assert(inStream.is_open());
  ULONG Nkr=10000;
  auto dkr=static_cast<real_prec>(1.e-3);
  fftw_array<real_prec> kr(Nkr), pkr(Nkr);
  for(ULONG i=0;i<Nkr;i++) // EGP: int i -> ULONG i
  {
    inStream >> kr[i] >> pkr[i];
  }
  inStream.close();
  real_prec minkr=kr[0];

  real_prec kNL=0.;//0.25;//attention for z=0!! from fig 3 in http://arxiv.org/pdf/1203.5785.pdf

  {
    real_prec dvar=0.;
    ULONG ikk=0;
    while (dvar<1.)
    {
      dvar+=static_cast<real_prec>(4.*M_PI*dkr*kr[ikk]*kr[ikk]*(pkr[ikk]/(4.*M_PI)));
      ikk++;
    }
    if (ikk>0) {
      std::cout << "delta^2=1-> kNL=" << kr[ikk - 1] << std::endl;
    }

    kNL=kr[ikk-1];
  }


#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (unsigned i=0;i<N1;i++)
    for (unsigned j=0;j<N2;j++)
      for (unsigned k=0;k<N3;k++)
      {
        ULONG ii=k+N3*(j+N2*i);

        real_prec kmod=sqrt(k_squared(i,j,k,L1,L2,L3,N1,N2,N3));

        if (zeld==true)
        {
          auto fac=static_cast<real_prec>(0.085);
          out[ii]=exp(-fac*kmod*kmod/(kNL*kNL));
        }

        if (secordlpt==true)
        {
          real_prec kmodnorm1=kmod/kNL;
          real_prec kmodnorm2=kmodnorm1*kmodnorm1;
          real_prec kmodnorm3=kmodnorm2*kmodnorm1;
          real_prec kmodnorm4=kmodnorm3*kmodnorm1;

          out[ii]=static_cast<real_prec>(exp(0.6*kmodnorm1-1.7*kmodnorm2+0.623*kmodnorm3-0.078*kmodnorm4));
        }

        if (denslpt==true)
        {
          //real_prec Vol=L1*L2*L3;
          //real_prec Norm;
//#ifdef FOURIER_DEF_1
          //Norm=1./Vol;
//#endif
//#ifdef FOURIER_DEF_2
          //Norm=real_prec(N)*real_prec(N)/Vol;
//#endif
          real_prec dvar=0.;
          ULONG ikr=0;
          if (kmod*0.5>minkr)
          {
            ikr=static_cast<ULONG>(floor((kmod*0.5-minkr)/dkr));//attention for z=0!!
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:dvar)
#endif // MULTITHREAD
            for (ULONG ikk=0;ikk<ikr;ikk++)
              dvar+=static_cast<real_prec>(dkr*kr[ikk]*kr[ikk]*(pkr[ikk]/(4.*M_PI)));//normalisation??
          }

          auto dmod=static_cast<real_prec>(4.*M_PI*dvar);

          out[ii]=static_cast<real_prec>(exp(0.58*dmod));
        }
      }


  {
    fftw_array<complex_prec>  AUX(N);
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
    for(ULONG i=0;i<N;i++)
    {
      re(AUX[i])=out[i];
      im(AUX[i])=0.;
    }
    FFT3d (N1,N2,N3, to_Rspace, AUX, AUX);

    real_prec wtot=0.;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:wtot)
#endif // MULTITHREAD
    for(ULONG i=0;i<N;i++)
      wtot+=re(AUX[i]);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
    for(ULONG i=0;i<N;i++)
    {
      out[i]/=wtot;
    }
  }

  std::string fname;
  switch (filtertype)
  {
    case 1:
      fname= dir + std::string("auxtransfzeld");
      break;
    case 2:
      fname= dir + std::string("auxtransf2lpt");
      break;
    case 3:
      fname= dir + std::string("auxtransfdens");
      break;
  }

  dump_scalar(out, N1, N2, N3, fname);
}
