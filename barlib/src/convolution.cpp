/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <cmath>

#include "fftw_array.h"
#include "fftwrapper.h"
#include "scale_space.hpp"
#include "struct_main.h"
#include "IOfunctionsGen.h"

#include "convolution.hpp"

////////////////////////////////////////////////////////////////
// Convolution and kernel functions
////////////////////////////////////////////////////////////////

void convolve(real_prec L1, real_prec L2, real_prec L3, unsigned N1, unsigned N2, unsigned N3, const real_prec *in, real_prec *out,
              real_prec smol, bool zeropad, int filtertype)
{
  bool gauss=false;
  bool tophat=false;
  bool errfunc=false;

  switch (filtertype)
  {
    case 1:
      gauss=true;
      break;
    case 2:
      tophat=true;
      break;
    case 3:
      errfunc=true;
      break;
  }

  //EGP  ULONG N=N1*N2*N3;

  unsigned Nzp1=N1;
  unsigned Nzp2=N2;
  unsigned Nzp3=N3;

  real_prec Lzp1=L1;
  real_prec Lzp2=L2;
  real_prec Lzp3=L3;

  if (zeropad==true)
  {
    Lzp1=L1*num_2;
    Lzp2=L2*num_2;
    Lzp3=L3*num_2;

    Nzp1=N1*2;
    Nzp2=N2*2;
    Nzp3=N3*2;
  }

  ULONG Nzp=Nzp1*Nzp2*Nzp3;

  fftw_array<complex_prec>  AUX(Nzp), AUXb(Nzp);
  fftw_array<real_prec> wkernel(Nzp),dummy(Nzp);

  if (zeropad==true)
  {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
    for(ULONG i=0;i<Nzp;i++)
      dummy[i]=0.0;

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
    for(unsigned i=0;i<N1;i++)
      for(unsigned j=0;j<N2;j++)
        for(unsigned k=0;k<N3;k++)
        {
          unsigned k_zp=k+N3/2;
          unsigned j_zp=Nzp3*(j+N2/2);
          unsigned i_zp=Nzp3*Nzp2*(i+N1/2);

          unsigned k_eff=k;
          unsigned j_eff=N3*j;
          unsigned i_eff=N3*N2*i;

          dummy[i_zp+j_zp+k_zp]=in[i_eff+j_eff+k_eff];
        }
  }
  else
  {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
    for(ULONG i=0;i<Nzp;i++)
      dummy[i]=in[i];
  }

  FFT3dR2C (Nzp1,Nzp2,Nzp3,dummy,AUXb);

  //EGP  real_prec asmth=1.0;
  real_prec u;

#ifdef MULTITHREAD
#pragma omp parallel for private(u)
#endif // MULTITHREAD
  for (unsigned i=0;i<Nzp1;i++)
    for (unsigned j=0;j<Nzp2;j++)
      for (unsigned k=0;k<Nzp3;k++)
      {
        real_prec k2=k_squared(i,j,k,Lzp1,Lzp2,Lzp3,Nzp1,Nzp2,Nzp3);

        real_prec rS=smol;
        real_prec rS2=rS*rS;
        real_prec kcut=smol;//2.*M_PI/rS;
        auto sigma=static_cast<real_prec>(.3);

        if (tophat==true)
        {
          u = std::sqrt(k2);

          if (u>kcut)
            wkernel[k+(Nzp3)*(j+Nzp2*i)]=0.0;
          else
            wkernel[k+(Nzp3)*(j+Nzp2*i)]=1.0;
        }

        if (errfunc==true)
        {
          u = static_cast<real_prec>((std::sqrt(k2)-kcut)/(std::sqrt(2.)*sigma));
          auto fac = static_cast<real_prec>(std::erfc(u));
          wkernel[k+(Nzp3)*(j+Nzp2*i)]=fac;
        }

        if (gauss==true)
          wkernel[k+(Nzp3)*(j+Nzp2*i)]=static_cast<real_prec>(std::exp(-k2*rS2/2.));
      }


#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for(ULONG i=0;i<Nzp;i++)
  {
    re(AUX[i])=wkernel[i];
    im(AUX[i])=0.0;
  }

  {
    FFT3d (Nzp1,Nzp2,Nzp3, to_Rspace, AUX, AUX);

    real_prec wtot=0.;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:wtot)
#endif // MULTITHREAD
    for(ULONG i=0;i<Nzp;i++)
      wtot+=re(AUX[i]);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
    for(ULONG i=0;i<Nzp;i++)
    {
      wkernel[i]/=wtot;
    }
  }


#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for(ULONG i=0;i<Nzp;i++)
  {
    re(AUX[i])=wkernel[i];
    im(AUX[i])=0.0;
  }

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for(ULONG i=0;i<Nzp;i++)
  {
    re(AUXb[i])= re(AUXb[i])*re(AUX[i]);
    im(AUXb[i])= im(AUXb[i])*re(AUX[i]);
  }

  FFT3d (Nzp1,Nzp2,Nzp3, to_Rspace, AUXb, AUXb);
  if (zeropad==true)
  {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
    for(unsigned i=0;i<N1;i++)
      for(unsigned j=0;j<N2;j++)
        for(unsigned k=0;k<N3;k++)
        {
          unsigned k_zp=k+N3/2;
          unsigned j_zp=Nzp3*(j+N2/2);
          unsigned i_zp=Nzp3*Nzp2*(i+N1/2);

          unsigned k_eff=k;
          unsigned j_eff=N3*j;
          unsigned i_eff=N3*N2*i;

          out[i_eff+j_eff+k_eff]=re(AUXb[i_zp+j_zp+k_zp]);
        }
  }
  else
  {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
    for(ULONG i=0;i<Nzp;i++)
      out[i]=re(AUXb[i]);
  }
}


void kernelcomp(real_prec L1, real_prec L2, real_prec L3, unsigned N1, unsigned N2, unsigned N3, real_prec smol, int filtertype,
                struct DATA *data)
{
  bool gauss=false;
  bool errfunc=false;
  bool tophat=false;

  switch (filtertype)
  {
    case 1:
      gauss=true;
      break;
    case 2:
      tophat=true;
      break;
    case 3:
      errfunc=true;
      break;
  }

  ULONG N=N1*N2*N3;

  fftw_array<real_prec> out(N);
  fftw_array<complex_prec> AUX(N);

  //EGP  real_prec asmth=1.0;
  real_prec u;

#ifdef MULTITHREAD
#pragma omp parallel for private(u)
#endif // MULTITHREAD
  for (unsigned i=0;i<N1;i++)
    for (unsigned j=0;j<N2;j++)
      for (unsigned k=0;k<N3;k++)
      {
        ULONG ii=k+N3*(j+N2*i);

        real_prec k2=k_squared(i,j,k,L1,L2,L3,N1,N2,N3);

        real_prec rS=smol;
        real_prec rS2=rS*rS;
        real_prec kcut=smol;//2.*M_PI/rS;
        auto sigma=static_cast<real_prec>(.3);

        if (tophat==true)
        {
          u = sqrt(k2);

          if (u>kcut)
            re(AUX[ii])=0.0;
          else
            re(AUX[ii])=1.0;
        }

        if (errfunc==true)
        {
          u = static_cast<real_prec>((std::sqrt(k2)-kcut)/(std::sqrt(2.)*sigma));
          auto fac = static_cast<real_prec>(std::erfc(u));
          re(AUX[ii])=fac;
        }

        if (gauss==true)
          re(AUX[ii])=static_cast<real_prec>(std::exp(-k2*rS2/2.));


        im(AUX[ii])=0.0;
      }


#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for(ULONG i=0;i<N;i++)
  {
    out[i]=re(AUX[i]);
    im(AUX[i])=0.0;
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

  int bmax=100;
  char buffsl[bmax];
  sprintf(buffsl,"r%d",int(smol));
  std::string fname = data->numerical->dir + "auxkernel" + buffsl;
  dump_scalar(out, N1, N2, N3, fname);
}


void convcomp(unsigned N1, unsigned N2, unsigned N3, real_prec *in, real_prec *out, real_prec smol, const std::string & dir)
{
  //bool gauss=false;
  //bool errfunc=false;
  //bool tophat=false;

  //switch (filtertype)
  //{
  //case 1:
  //gauss=true;
  //break;
  //case 2:
  //tophat=true;
  //break;
  //case 3:
  //errfunc=true;
  //break;
  //}


  ULONG N=N1*N2*N3;

  fftw_array<complex_prec> AUX(N);

  FFT3dR2C (N1,N2,N3,in,AUX);

  {
    int bmax=100;
    char buffsl[bmax];
    sprintf(buffsl,"r%d",int(smol));
    //EGP    char * fileN;
    string fname= dir + string("auxkernel")+buffsl;
    get_scalar(fname,out,N1,N2,N3);
  }

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (unsigned i=0;i<N1;i++)
    for (unsigned j=0;j<N2;j++)
      for (unsigned k=0;k<N3;k++)
      {
        ULONG kk=k+(N3)*(j+N2*i);
        ULONG ii=k+(N3)*(j+N2*i);

        re(AUX[kk])= re(AUX[kk])*out[ii];
        im(AUX[kk])= im(AUX[kk])*out[ii];
      }

  FFT3dC2R (N1,N2,N3,AUX,out);
}


void convcompb(unsigned N1, unsigned N2, unsigned N3, real_prec *in, real_prec *out)
{
  ULONG N=N1*N2*N3;

  fftw_array<complex_prec> AUX(N);

  FFT3dR2C (N1,N2,N3,in,AUX);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (unsigned i=0;i<N1;i++)
    for (unsigned j=0;j<N2;j++)
      for (unsigned k=0;k<N3;k++)
      {
        ULONG kk=k+(N3)*(j+N2*i);
        ULONG ii=k+(N3)*(j+N2*i);

        re(AUX[kk])= re(AUX[kk])*out[ii];
        im(AUX[kk])= im(AUX[kk])*out[ii];
      }

  FFT3dC2R (N1,N2,N3,AUX,out);
}

