/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include <climits>  // INT_MAX
#include <exception>  // std::exception
#include <iostream>
#include <string>

#include "define_opt.h"
#include "fftw_array.h"

#include "math_funcs.h"

#include "convenience.h"

#include "fftwrapper.h"


// Take care: multidimensional C2R fftw transforms destroy the input array!
void fftC2R(unsigned int N1, unsigned int N2, unsigned int N3, complex_prec *in, real_prec *out)
{
  if (N1 > INT_MAX || N2 > INT_MAX || N3 > INT_MAX) {
    std::string msg = "FFT dimensions must not be larger than INT_MAX!";
    throw std::logic_error(msg);
  }
#ifdef SINGLE_PREC
  fftwf_plan fftp;
  fftp = fftwf_plan_dft_c2r_3d(static_cast<int>(N1), static_cast<int>(N2), static_cast<int>(N3), in,out,FFTW_OPTION);
  fftwf_execute(fftp);
#endif  
#ifdef DOUBLE_PREC
  fftw_plan fftp;
  fftp = fftw_plan_dft_c2r_3d(static_cast<int>(N1), static_cast<int>(N2), static_cast<int>(N3), in,out,FFTW_OPTION);
  fftw_execute(fftp);
#endif  

  ULONG  N=N1*N2*N3;
  real_prec fac = 1 / static_cast<real_prec>(N);
  multiply_factor_array(fac, out, out, N);

#ifdef SINGLE_PREC
  fftwf_destroy_plan(fftp);
#endif
#ifdef DOUBLE_PREC
  fftw_destroy_plan(fftp);
#endif
}


void fftR2C(unsigned int N1, unsigned int N2, unsigned int N3, real_prec *in, complex_prec *out)
{
  if (N1 > INT_MAX || N2 > INT_MAX || N3 > INT_MAX) {
    std::string msg = "FFT dimensions must not be larger than INT_MAX!";
    throw std::logic_error(msg);
  }
#ifdef SINGLE_PREC
  fftwf_plan fftp;
  fftp = fftwf_plan_dft_r2c_3d(static_cast<int>(N1), static_cast<int>(N2), static_cast<int>(N3), in,out,FFTW_OPTION);
  fftwf_execute(fftp);
#endif
#ifdef DOUBLE_PREC
  fftw_plan fftp;
  fftp = fftw_plan_dft_r2c_3d(static_cast<int>(N1), static_cast<int>(N2), static_cast<int>(N3), in,out,FFTW_OPTION);
  fftw_execute(fftp);
#endif

#ifdef SINGLE_PREC
  fftwf_destroy_plan(fftp);
#endif
#ifdef DOUBLE_PREC
  fftw_destroy_plan(fftp);
#endif
}





// Real FFT versions that take a plan as input. This greatly speeds up the FFT

// Take care: multidimensional C2R fftw transforms destroy the input array!
void fftC2Rplanned(complex_prec *in, real_prec *out, struct plan_pkg *plan) {
  if (in != plan->C) {
    copyArray(in, plan->C, plan->Nhalf);
  }

  #ifdef SINGLE_PREC
  fftwf_execute(plan->plan);
  #endif
  #ifdef DOUBLE_PREC
  fftw_execute(plan->plan);
  #endif

  real_prec fac = 1 / static_cast<real_prec>(plan->N);
  multiply_factor_array(fac, plan->R, out, plan->N);
}

void fftR2Cplanned(real_prec *in, complex_prec *out, struct plan_pkg *plan) {
  if (in != plan->R) {
    copyArray(in, plan->R, plan->N);
  }

  #ifdef SINGLE_PREC
  fftwf_execute(plan->plan);
  #endif
  #ifdef DOUBLE_PREC
  fftw_execute(plan->plan);
  #endif

  if (out != plan->C) {
    copyArray(plan->C, out, plan->Nhalf);
  }
}









void FFT3d(unsigned int N1, unsigned int N2, unsigned int N3, bool direction, complex_prec *in, complex_prec *out)
{
  if (N1 > INT_MAX || N2 > INT_MAX || N3 > INT_MAX) {
    std::string msg = "FFT dimensions must not be larger than INT_MAX!";
    throw std::logic_error(msg);
  }

  ULONG factor=N1*N2*N3;

#ifdef SINGLE_PREC
  fftwf_plan fftp;
  if (direction == true) // EGP: to Fourier space
  {
    fftp = fftwf_plan_dft_3d(N1,N2,N3,in,out,FORWARD,FFTW_OPTION);	
    fftwf_execute(fftp);
  }
  else // EGP: back to real space
  {
    fftp = fftwf_plan_dft_3d(N1,N2,N3,in,out,BACKWARD,FFTW_OPTION);	
    fftwf_execute(fftp);
    complexfactor_mult(factor,static_cast<real_prec>(1./real_prec(factor)),out,out);
  }
#endif

#ifdef DOUBLE_PREC
  fftw_plan fftp;
  if (direction == true) // EGP: to Fourier space
  {
    fftp = fftw_plan_dft_3d(static_cast<int>(N1), static_cast<int>(N2), static_cast<int>(N3),in,out,FORWARD,FFTW_OPTION);
    fftw_execute(fftp);
  }
  else // EGP: back to real space
  {
    fftp = fftw_plan_dft_3d(static_cast<int>(N1), static_cast<int>(N2), static_cast<int>(N3),in,out,BACKWARD,FFTW_OPTION);
    fftw_execute(fftp);
    complexfactor_mult(factor,1./real_prec(factor),out,out);
  }
#endif
#ifdef SINGLE_PREC
  fftwf_destroy_plan(fftp);
#endif
#ifdef DOUBLE_PREC
  fftw_destroy_plan(fftp);
#endif
}


void FFT3dR2C(unsigned int N1, unsigned int N2, unsigned int N3, real_prec *in, complex_prec *out)
{
  if (N1 > INT_MAX || N2 > INT_MAX || N3 > INT_MAX) {
    std::string msg = "FFT dimensions must not be larger than INT_MAX!";
    throw std::logic_error(msg);
  }

  ULONG N=N1*N2*N3;

  fftw_array<complex_prec> inC(N);

  for(ULONG i=0;i<N;i++){
    re(inC[i])=in[i];
    im(inC[i])=0.;
  }

#ifdef SINGLE_PREC
  fftwf_plan fftp;  
  fftp = fftwf_plan_dft_3d(N1,N2,N3,inC,out,FORWARD,FFTW_OPTION);	
  fftwf_execute(fftp);
#endif
#ifdef DOUBLE_PREC
  fftw_plan fftp;  
  fftp = fftw_plan_dft_3d(static_cast<int>(N1), static_cast<int>(N2), static_cast<int>(N3),inC,out,FORWARD,FFTW_OPTION);
  fftw_execute(fftp);
#endif

#ifdef SINGLE_PREC
  fftwf_destroy_plan(fftp);
#endif
#ifdef DOUBLE_PREC
  fftw_destroy_plan(fftp);
#endif
}

void FFT3dC2R(unsigned int N1, unsigned int N2, unsigned int N3, complex_prec *in, real_prec *out)
{
  if (N1 > INT_MAX || N2 > INT_MAX || N3 > INT_MAX) {
    std::string msg = "FFT dimensions must not be larger than INT_MAX!";
    throw std::logic_error(msg);
  }

  ULONG N=N1*N2*N3;

  fftw_array<complex_prec> inC(N);

  for(ULONG i=0;i<N;i++){
    re(inC[i])=re(in[i]);
    im(inC[i])=im(in[i]);
  }  

#ifdef SINGLE_PREC
  fftwf_plan fftp;  
  fftp = fftwf_plan_dft_3d(N1,N2,N3,inC,inC,BACKWARD,FFTW_OPTION);	
  fftwf_execute(fftp);
#endif
#ifdef DOUBLE_PREC
  fftw_plan fftp;  
  fftp = fftw_plan_dft_3d(static_cast<int>(N1), static_cast<int>(N2), static_cast<int>(N3),inC,inC,BACKWARD,FFTW_OPTION);
  fftw_execute(fftp);
#endif

  complexfactor_mult(N,static_cast<real_prec>(1./real_prec(N)),inC,inC);
#ifdef SINGLE_PREC
  fftwf_destroy_plan(fftp);
#endif
#ifdef DOUBLE_PREC
  fftw_destroy_plan(fftp);
#endif

  for(ULONG i=0;i<N;i++){
    out[i]=re(inC[i]);
  }  
}


//void FFTW_accumulate_wisdom (int N1,int N2,int N3)
//{
  //ULONG N=N1*N2*N3;

  //fftw_array<complex_prec> inC(N), out(N);
  //fill_one(inC, N);

  //// Options below: FFTW_MEASURE, _PATIENT and _EXHAUSTIVE

//#ifdef SINGLE_PREC
  //fftwf_plan fftp;  
  //fftp = fftwf_plan_dft_3d(N1,N2,N3,inC,out,FORWARD,FFTW_MEASURE);	
  //fftwf_execute(fftp);
//#endif
//#ifdef DOUBLE_PREC
  //fftw_plan fftp;  
  //fftp = fftw_plan_dft_3d(N1,N2,N3,inC,out,FORWARD,FFTW_MEASURE);	
  //fftw_execute(fftp);
//#endif

//#ifdef SINGLE_PREC
  //fftwf_destroy_plan(fftp);
//#endif
//#ifdef DOUBLE_PREC
  //fftw_destroy_plan(fftp);
//#endif
//}

// R2C
plan_pkg::plan_pkg(unsigned int _N1, unsigned int _N2, unsigned int _N3, real_prec *in, complex_prec *out) {
  std::cout << "plan_pkg R2C" << std::endl;

  if (_N1 > INT_MAX || _N2 > INT_MAX || _N3 > INT_MAX) {
    std::string msg = "FFT dimensions must not be larger than INT_MAX!";
    throw std::logic_error(msg);
  }

  N1 = _N1, N2 = _N2, N3 = _N3;
  N = (static_cast<ULONG>(N1) * N2) * N3;
  Nhalf = (static_cast<ULONG>(N1) * N2) * (N3/2 + 1);
  R = in;
  C = out;

#ifdef SINGLE_PREC
  plan = fftwf_plan_dft_r2c_3d(static_cast<int>(N1), static_cast<int>(N2), static_cast<int>(N3), in, out, FFTW_PATIENT);
#endif
#ifdef DOUBLE_PREC
  plan = fftw_plan_dft_r2c_3d(static_cast<int>(N1), static_cast<int>(N2), static_cast<int>(N3), in, out, FFTW_PATIENT);
#endif
}

// C2R
plan_pkg::plan_pkg(unsigned int _N1, unsigned int _N2, unsigned int _N3, complex_prec *in, real_prec *out) {
  std::cout << "plan_pkg C2R" << std::endl;

  if (_N1 > INT_MAX || _N2 > INT_MAX || _N3 > INT_MAX) {
    std::string msg = "FFT dimensions must not be larger than INT_MAX!";
    throw std::logic_error(msg);
  }

  N1 = _N1, N2 = _N2, N3 = _N3;
  N = static_cast<ULONG>(N1*N2*N3);
  Nhalf = static_cast<ULONG>(N1 * N2 * (N3/2 + 1));
  C = in;
  R = out;

#ifdef SINGLE_PREC
  plan = fftwf_plan_dft_c2r_3d(static_cast<int>(N1), static_cast<int>(N2), static_cast<int>(N3), in, out, FFTW_PATIENT);
#endif
#ifdef DOUBLE_PREC
  plan = fftw_plan_dft_c2r_3d(static_cast<int>(N1), static_cast<int>(N2), static_cast<int>(N3), in, out, FFTW_PATIENT);
#endif
}

plan_pkg::~plan_pkg() {
#ifdef SINGLE_PREC
  fftwf_destroy_plan(plan);
#endif
#ifdef DOUBLE_PREC
  fftw_destroy_plan(plan);
#endif
}
