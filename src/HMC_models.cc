/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include "struct_main.h"
#include "struct_hamil.h"

#include <cmath>
#include <iomanip>
#include <algorithm>  // std::max_element, std::find
#include <iterator>  // std::distance

#include <gsl/gsl_integration.h>

#include "fftw_array.h"

#include "math_funcs.h"
#include "Lag2Eul.h"
#include "cosmo.h"

#include "HMC_help.h"

#include "convenience.h"

using namespace std;

// TODO: split this file into different files for each model:
//       HMC/model/gaussian.cc, HMC/model/poisson.cc, etc.

// Note: in HMC_{momenta,mass}.cc the mass is tuned (if mass_fs)
// to the Gaussian likelihood model!


// Model: Poissonian likelihood //
void poissonian_likelihood_partial_f_delta_x_log_like(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *dummy)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;

  for(ULONG i = 0; i < n->N; i++)
  {
    auto dens = static_cast<real_prec>(1. + hd->biasP * deltaX[i]);// must be positive!
    real_prec Lambda = hd->window[i] * hd->rho_c * pow(dens, hd->biasE);

    if ((hd->window[i] > 0.0) && (dens > 0.0))
      //dummy[i] = (Lambda - hd->nobs[i]) * (hd->biasE * hd->biasP / dens);
      dummy[i] = (1 - hd->nobs[i]/Lambda) * hd->rho_c * hd->biasE * hd->biasP * pow(dens, hd->biasE - 1);
    else
      dummy[i] = 0.0;
  }
}

// N.B.: grad_f_delta_x_comp should no longer be used. Only necessary for calc_h != 2, which is bogus!
void poissonian_likelihood_grad_f_delta_x_comp(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out,
                                               unsigned int component)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  gradfindif(n->N1, n->L1, deltaX, out, component);
}

real_prec poissonian_likelihood_log_like(struct HAMIL_DATA *hd, real_prec *delta)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;

  // calculate deltaX and pos
  {
    unsigned facL=1;
    bool reggrid=true;
    gsl_rng *seed = nullptr; // empty: reggrid is true anyway
    real_prec kernel_scale = SPH_kernel_scale(hd);
    Lag2Eul(delta, hd->deltaX, hd->posx, hd->posy, hd->posz, n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1, n->d2,
            n->d3, n->min1, n->min2, n->min3, hd->D1, hd->D2, hd->ascale, hd->OM, hd->OL, hd->sfmodel, n->mk, n->kth, facL,
            reggrid, seed, "", kernel_scale, n->R2Cplan, n->C2Rplan);
  }

  // actual likelihood calculation
  real_prec out = 0.;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:out)
#endif // MULTITHREAD 
  for(ULONG i = 0; i < n->N; i++)
  {
    auto dens = static_cast<real_prec>(1. + hd->biasP * hd->deltaX[i]);// must be positive!
    real_prec Lambda = hd->window[i] * hd->rho_c * pow(dens, hd->biasE);

    if ((hd->window[i] > 0.) && (Lambda > 0.0))
      out += Lambda - hd->nobs[i] * log(Lambda);
  }

  return (out);
}

// END model Poissonian likelihood //


// Model: Gaussian random field (no evolution, just lagrangian field) //
// N.B.: grad_f_delta_x_comp should no longer be used. Only necessary for calc_h != 2, which is bogus!
void grf_likelihood_grad_f_delta_x_comp(struct HAMIL_DATA *, real_prec *, real_prec *, unsigned int)
{
}

void grf_likelihood_partial_f_delta_x_log_like(struct HAMIL_DATA *, real_prec *, real_prec *)
{
}

void grf_likelihood_grad_log_like(struct HAMIL_DATA *hd, real_prec *delta, real_prec *out)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;

#ifdef MULTITHREAD
#pragma omp parallel for 
#endif // MULTITHREAD 
  for(ULONG i=0;i<n->N;i++)
    if (hd->window[i]>0.)
      out[i] = (delta[i] - hd->nobs[i]) / (hd->noise[i] * hd->noise[i]);
    else
      out[i] = 0;
}

real_prec grf_likelihood_log_like(struct HAMIL_DATA *hd, real_prec *delta)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;

  real_prec out=0.;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:out)
#endif // MULTITHREAD 
  for(ULONG i=0;i<n->N;i++)
    if (hd->window[i]>0.)
      out += num_0_5 * gsl_pow_2((delta[i] - hd->nobs[i])/hd->noise[i]);

  return (out);
}

// END model Gaussian random field //

  

// Model: Gaussian likelihood //
void gaussian_likelihood_partial_f_delta_x_log_like_interpolate(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  // EGP: partial -log L/partial delta_x = (Lambda-nobs)/sigma**2
  
  fftw_array<real_prec> deltaXQ(n->N), nobsXQ(n->N), noiseXQ(n->N), windowXQ(n->N);

  interpolate_CIC(n->N1,n->N2,n->N3,n->L1,n->L2,n->L3,n->d1,n->d2,n->d3, hd->posx, hd->posy, hd->posz, deltaX, n->N, deltaXQ);
  interpolate_CIC(n->N1,n->N2,n->N3,n->L1,n->L2,n->L3,n->d1,n->d2,n->d3, hd->posx, hd->posy, hd->posz, hd->nobs, n->N, nobsXQ);
  interpolate_CIC(n->N1,n->N2,n->N3,n->L1,n->L2,n->L3,n->d1,n->d2,n->d3, hd->posx, hd->posy, hd->posz, hd->noise, n->N, noiseXQ);
  interpolate_CIC(n->N1,n->N2,n->N3,n->L1,n->L2,n->L3,n->d1,n->d2,n->d3, hd->posx, hd->posy, hd->posz, hd->window, n->N, windowXQ);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD 
  for(ULONG i=0;i<n->N;i++)
  {
    real_prec Lambda = windowXQ[i] * hd->rho_c * pow(num_1 + hd->biasP * deltaXQ[i], hd->biasE);
    if ((windowXQ[i]>0.) && (Lambda > 0.0))
    {
      real_prec resid = nobsXQ[i]-Lambda;
      out[i] = resid/(noiseXQ[i]*noiseXQ[i]);
    }
    else
      out[i] = 0.0;
  }
}
// N.B.: grad_f_delta_x_comp should no longer be used. Only necessary for calc_h != 2, which is bogus!
void gaussian_likelihood_grad_f_delta_x_comp_interpolate(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out,
                                                         unsigned int component)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  fftw_array<real_prec> deltaXQ(n->N);
  interpolate_CIC(n->N1,n->N2,n->N3,n->L1,n->L2,n->L3,n->d1,n->d2,n->d3, hd->posx, hd->posy, hd->posz, deltaX, n->N, deltaXQ);

  gradfft(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, deltaXQ, out, component);
  //gradfindif(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, deltaX, out, component);
}

// NOTE:
// partial_f_delta_x_log_like is also used in Gaussian RSD likelihood!
// Take care when introducing rank ordering, then it might become different.
void gaussian_likelihood_partial_f_delta_x_log_like(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  // EGP: partial -log L/partial delta_x = (Lambda-nobs)/sigma**2
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD 
  for(ULONG i = 0; i < n->N; i++)
  {
    real_prec Lambda = hd->window[i] * hd->rho_c * pow(num_1 + hd->biasP * deltaX[i], hd->biasE);
    if ((hd->window[i] > 0.) && (Lambda > 0.0))
    {
      real_prec resid = hd->nobs[i] - Lambda;
      out[i] = resid / (hd->noise[i] * hd->noise[i]);
    }
    else
      out[i]=0.0;
  }
}
// N.B.: grad_f_delta_x_comp should no longer be used. Only necessary for calc_h != 2, which is bogus!
void gaussian_likelihood_grad_f_delta_x_comp(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out,
                                             unsigned int component)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  //gradfindif(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, deltaX, out, component);
  gradfft(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, deltaX, out, component);
}
real_prec gaussian_likelihood_log_like(struct HAMIL_DATA *hd, real_prec *deltaQ)
{
  // -log L = 0.5*(Lambda-nobs)^2/sigma^2
  struct HAMIL_NUMERICAL *n = hd->numerical;
  fftw_array<real_prec> delta_growing(n->N);

  // single out growing mode (Peebles) -> delta_growing
  multiply_factor_array(n->deltaQ_factor, deltaQ, delta_growing, n->N);

  // calculate deltaX and pos
  {
    unsigned facL = 1;
    bool reggrid = true;
    gsl_rng *seed = nullptr; // empty: reggrid is true anyway
    real_prec kernel_scale = SPH_kernel_scale(hd);
    if (hd->rsd_model)
      Lag2Eul_rsd_zeldovich(delta_growing, hd->deltaX, hd->posx, hd->posy, hd->posz, n->N1, n->N2, n->N3, n->L1, n->L2,
                            n->L3, n->d1, n->d2, n->d3, n->min1, n->min2, n->min3, hd->D1, hd->ascale, hd->OM, hd->OL,
                            n->mk, facL,
                            reggrid, seed, kernel_scale, n->xobs, n->yobs, n->zobs, n->planepar, n->periodic,
                            n->R2Cplan, n->C2Rplan);
    else
      Lag2Eul(delta_growing, hd->deltaX, hd->posx, hd->posy, hd->posz, n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1,
              n->d2, n->d3, n->min1, n->min2, n->min3, hd->D1, hd->D2, hd->ascale, hd->OM, hd->OL, hd->sfmodel, n->mk, n->kth,
              facL,
              reggrid, seed, "", kernel_scale, n->R2Cplan, n->C2Rplan);
  }

  // actual likelihood calculation
  real_prec out=0.;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:out)
#endif // MULTITHREAD
  for(ULONG i=0;i<n->N;i++)
  {
    real_prec Lambda = hd->window[i] * hd->rho_c * pow(num_1+hd->biasP*hd->deltaX[i],hd->biasE);
    if ((hd->window[i]>0.) && (Lambda > 0.0))
      out+=num_0_5 * gsl_pow_2( (Lambda - hd->nobs[i]) / hd->noise[i] );
  }

  return (out);
}
// End model: Gaussian likelihood  //



// Model: log-normal likelihood //

// TEST: with lambda function (c++11) for Lambda 
//       This is also a useful structure for when/if we put these models in
//       full-fledged classes; one could then define the function for Lambda
//       as a class method and reuse it in the different functions.
//void lognormal_likelihood_partial_f_delta_x_log_like(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *dummy)
//{
  //struct HAMIL_NUMERICAL *n = hd->numerical;
  //// EGP: partial -log L/partial log(1+delta_x) = (Lambda-nobs)/sigma**2
  //auto Lambda = [&] (ULONG index) -> real_prec {log(hd->rho_c * pow(num_1 + hd->biasP * deltaX[index], hd->biasE));};
//#ifdef MULTITHREAD
//#pragma omp parallel for
//#endif // MULTITHREAD
  //for(ULONG i=0; i < n->N;i++)
    //if (hd->window[i]>0.)
      //dummy[i] = (hd->nobs[i] - Lambda(i))/(hd->noise[i] * hd->noise[i]);
    //else
      //dummy[i] = 0.0;
//}

void lognormal_likelihood_partial_f_delta_x_log_like(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *dummy)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  // EGP: partial -log L/partial log(1+delta_x) = (Lambda-nobs)/sigma**2
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for(ULONG i=0; i < n->N;i++)
  {
    real_prec Lambda = log(hd->rho_c * pow(num_1 + hd->biasP * deltaX[i], hd->biasE));
    if (hd->window[i]>0.)
      dummy[i]=(hd->nobs[i] - Lambda)/(hd->noise[i] * hd->noise[i]);
    else
      dummy[i]=0.0;
  }
}

real_prec lognormal_likelihood_f_delta_x_i_calc(real_prec rho_c, real_prec delta_min, real_prec deltaX_i)
{
  // keep above zero density
  if (deltaX_i < delta_min)
    deltaX_i = delta_min;

  return log(rho_c * (num_1 + deltaX_i));
}

real_prec lognormal_likelihood_f_delta_x_i(struct HAMIL_DATA *hd, real_prec deltaX_i)
{
  return lognormal_likelihood_f_delta_x_i_calc(hd->rho_c, hd->delta_min, deltaX_i);
}

void lognormal_likelihood_f_delta_x(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (ULONG i = 0; i < n->N; i++)
    out[i] = lognormal_likelihood_f_delta_x_i(hd, deltaX[i]);
}

// N.B.: grad_f_delta_x_comp should no longer be used. Only necessary for calc_h != 2, which is bogus!
void lognormal_likelihood_grad_f_delta_x_comp(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out,
                                              unsigned int component)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  fftw_array<real_prec> f_delta_x(n->N);

  lognormal_likelihood_f_delta_x(hd, deltaX, f_delta_x);

  gradfindif(n->N1, n->L1, f_delta_x, out, component);
}

real_prec lognormal_likelihood_log_like(struct HAMIL_DATA *hd, real_prec *delta)
{
  // -log L = 0.5*(Lambda-nobs)^2/sigma^2
  struct HAMIL_NUMERICAL *n = hd->numerical;

  // calculate deltaX and pos
  {
    unsigned facL=1;
    bool reggrid=true;
    gsl_rng *seed = nullptr; // empty: reggrid is true anyway
    real_prec kernel_scale = SPH_kernel_scale(hd);
    Lag2Eul(delta, hd->deltaX, hd->posx, hd->posy, hd->posz, n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1, n->d2,
            n->d3, n->min1, n->min2, n->min3, hd->D1, hd->D2, hd->ascale, hd->OM, hd->OL, hd->sfmodel, n->mk, n->kth, facL,
            reggrid, seed, "", kernel_scale, n->R2Cplan, n->C2Rplan);
  }

  // actual likelihood calculation
  real_prec out=0.;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:out)
#endif // MULTITHREAD 
  for(ULONG i=0;i<n->N;i++)
  {
    real_prec Lambda = lognormal_likelihood_f_delta_x_i(hd, hd->deltaX[i]);
    if (hd->window[i]>0.)
    {
      real_prec resid=Lambda - hd->nobs[i];
      out+=num_0_5*resid*resid/(hd->noise[i]*hd->noise[i]);
    }
  }

  return (out);
}
// End model: log-normal likelihood //


// Model: Gaussian likelihood (or any other that is dependent on delta_q only through delta_x) //
void likelihood_calc_h(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  fftw_array<real_prec> partLike(n->N), dummy(n->N);
  fftw_array<complex_prec> outC(n->Nhalf), dummyC(n->Nhalf);

  fillZero(outC, n->Nhalf);

  hd->partial_f_delta_x_log_like(hd, deltaX, partLike);

  for (unsigned i = 1; i <= 3; i++){
    // compute i-th component of x-gradient of f(deltaX) -> dummy
    hd->grad_f_delta_x_comp(hd, deltaX, dummy, i);
    // multiply with partLike -> dummy
    multiplyArrays(partLike, dummy, dummy, n->N);
    // transform to FS => g_i(k) -> dummyC
    fftR2C(n->N1, n->N2, n->N3, dummy, dummyC);
    // multiply by -ik_i/k^2
    bool rfft = true;
    grad_inv_lap_FS(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, dummyC, dummyC, i, rfft); 
    // sum the resulting h(k) term to outC
    add_to_array(dummyC, outC, n->Nhalf);
  }
  // transform h(k) to real space -> h(q)
  fftC2R(n->N1, n->N2, n->N3, outC, out);
}
void likelihood_calc_h_old(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  fftw_array<real_prec> partLike(n->N), dummy(n->N), dummy2(n->N);
  fftw_array<complex_prec> outC(n->Nhalf), dummyC(n->Nhalf);

  fillZero(outC, n->Nhalf);

  hd->partial_f_delta_x_log_like(hd, deltaX, partLike);

  for (unsigned i = 1; i <= 3; i++){
    // compute i-th component of x-gradient of f(deltaX) -> dummy
    hd->grad_f_delta_x_comp(hd, deltaX, dummy, i);
    //if (n->iGibbs == 1)
      //fill_one(dummy, n->N); // testing: when delta_1 == 0 (which it is with zero initial guess), then deltaX and its derivative will be zero, breaking the whole algorithm for this first step

    // multiply with partLike -> dummy
    multiplyArrays(partLike, dummy, dummy, n->N);
    // interpolate to x(q_i) => g_i(q)
    //{
     //switch (n->mk)
      //{
        //case 0:
          //throw BarcodeException("likelihood_calc_h: NGP interpolation not yet implemented!");
          //break;
        //case 1:
          interpolate_CIC(n->N1,n->N2,n->N3,n->L1,n->L2,n->L3,n->d1,n->d2,n->d3, hd->posx, hd->posy, hd->posz,dummy,n->N,dummy2);
          //break;
        //case 2:
          //throw BarcodeException("likelihood_calc_h: TSC interpolation not yet implemented!");
          //break;
      //}
    //}
    // transform to FS => g_i(k) -> AUX
    fftR2C(n->N1, n->N2, n->N3, dummy2, dummyC);
    // multiply by -ik_i/k^2
    bool rfft = true;
    grad_inv_lap_FS(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, dummyC, dummyC, i, rfft);
    // sum the resulting h(k) term to outC
    add_to_array(dummyC, outC, n->Nhalf);
  }
  // transform h(k) to real space -> h(q)
  fftC2R(n->N1, n->N2, n->N3, outC, out);
}


void likelihood_calc_h_simple(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  fftw_array<real_prec> partLike(n->N);

  hd->partial_f_delta_x_log_like(hd, deltaX, partLike);

  // interpolate to x(q_i) => g_i(q)
  //switch (n->mk)
  //{
    //case 0:
      //throw BarcodeException("likelihood_calc_h: NGP interpolation not yet implemented!");
      //break;
    //case 1:
      interpolate_CIC(n->N1,n->N2,n->N3,n->L1,n->L2,n->L3,n->d1,n->d2,n->d3, hd->posx, hd->posy, hd->posz, partLike, n->N, out);
      //break;
    //case 2:
      //throw BarcodeException("likelihood_calc_h: TSC interpolation not yet implemented!");
      //break;
  //}
}


// new SPH method stuff //
real_prec SPH_kernel_radius(struct HAMIL_DATA *hd) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  return(SPH_kernel_radius(n->particle_kernel, n->particle_kernel_h));
}

real_prec SPH_kernel_radius(int particle_kernel_type,
                            real_prec particle_kernel_h) {
  real_prec radius;
  switch (particle_kernel_type)
  {
    case 0: // SPH kernel
      radius = particle_kernel_h * 2;
      break;
    default:
      throw BarcodeException("\"Eat my shorts!\" -- B. J. Simpson");
  }
  return(radius);
}

// double inline __attribute__((fastcall)) sqrt14(double n) {
//   __asm__ ("fld qword ptr [esp+4];"
//            "fsqrt;"
//            "ret 8");
//         // _asm fld qword ptr [esp+4]
//         // _asm fsqrt
//         // _asm ret 8
// }

void grad_SPH_kernel_3D(real_prec x, real_prec y, real_prec z, real_prec h_inv, real_prec h_sq_inv, real_prec norm,
                        real_prec &out_x, real_prec &out_y, real_prec &out_z) {
  // N.B.: when using this for density estimation, you need to normalize the
  // result afterwards with V/N! See e.g. getDensity_SPH. Same goes for
  // SPH_kernel_3D.

  // norm must be 1/(PI * h^4)
  
  // derivative of Monaghan kernel W_4
  real_prec r_sq = x*x + y*y + z*z;
  real_prec partial;
  // if (r_sq > 4*h_sq)
  //   partial = 0.;
  // else if (r_sq > h_sq) {
  //   real_prec r = sqrt(r_sq);
  //   real_prec q = r * h_inv;
  //   real_prec qmin2 = q - 2;
  //   partial = -0.75*qmin2*qmin2 * norm / r;
  // } else {
  //   real_prec r = sqrt(r_sq);
  //   partial = (2.25*r*h_inv*h_inv - 3*h_inv) * norm;
  // }
  auto q_sq_i = static_cast<int>(r_sq * h_sq_inv);
  switch (q_sq_i) {
    case 0:
      {
        real_prec r = sqrt(r_sq);
        partial = (2.25*r*h_inv*h_inv - 3*h_inv) * norm;
      }
      break;
    case 1:
    case 2:
    case 3:
      {
        real_prec r = sqrt(r_sq);
        real_prec q = r * h_inv;
        real_prec qmin2 = q - 2;
        partial = -0.75*qmin2*qmin2 * norm / r;
      }
      break;
    default:
      partial = 0.;
      break;
  }

  out_x = partial * x;
  out_y = partial * y;
  out_z = partial * z;
}


// #include <ia32intrin.h>


void grad_SPH_kernel_3D_h_units(real_prec x_h, real_prec y_h, real_prec z_h,
                                real_prec norm,
                        real_prec &out_x, real_prec &out_y, real_prec &out_z) {
  // N.B.: when using this for density estimation, you need to normalize the
  // result afterwards with V/N! See e.g. getDensity_SPH. Same goes for
  // SPH_kernel_3D.

  // norm must be 1/(PI * h^4)
  
  // derivative of Monaghan kernel W_4
  // real_prec r_sq = x*x + y*y + z*z;
  real_prec q_sq = x_h*x_h + y_h*y_h + z_h*z_h;
  // __m128 d_h = _mm_setr_ps(static_cast<float>(x_h), static_cast<float>(y_h),
                          // static_cast<float>(z_h), 0.0);
  // real_prec q_sq = static_cast<real_prec>(_mm_cvtss_f32(_mm_dp_ps(d_h, d_h, 0x71)));

  real_prec partial;
  if (q_sq > 4)
    partial = 0.;
  else if (q_sq > 1) {
    // real_prec r = sqrt(r_sq);
    // real_prec q = r * h_inv;
    // real_prec qmin2 = q - 2;
    // partial = -0.75*qmin2*qmin2 * norm / r;
    real_prec q = sqrt(q_sq);
    real_prec qmin2 = q - 2;
    partial = -0.75*qmin2*qmin2 * norm / q;
  } else {
    // real_prec r = sqrt(r_sq);
    // partial = (2.25*r*h_inv*h_inv - 3*h_inv) * norm;
    real_prec q = sqrt(q_sq);
    partial = (2.25*q - 3) * norm;
  }
  // int q_sq_i = static_cast<int>(q_sq);
  // switch (q_sq_i) {
  //   case 0:
  //     {
  //       real_prec q = sqrt(q_sq);
  //       partial = (2.25*q - 3) * norm;
  //     }
  //     break;
  //   case 1:
  //   case 2:
  //   case 3:
  //     {
  //       real_prec q = sqrt(q_sq);
  //       // real_prec q = r * h_inv;
  //       real_prec qmin2 = q - 2;
  //       // partial = -0.75*qmin2*qmin2 * norm / r;
  //       partial = -0.75*qmin2*qmin2 * norm / q;
  //     }
  //     break;
  //   default:
  //     partial = 0.;
  //     break;
  // }

  out_x = partial * x_h;
  out_y = partial * y_h;
  out_z = partial * z_h;
}


int SPH_kernel_3D_cells_count(struct HAMIL_DATA *hd) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  return(SPH_kernel_3D_cells_count(n->particle_kernel, n->particle_kernel_h,
                                   n->d1, n->d2, n->d3));
}

int SPH_kernel_3D_cells_count(int particle_kernel_type,
                              real_prec particle_kernel_h,
                              real_prec d1, real_prec d2, real_prec d3) {

  // real_prec kernel_reach = SPH_kernel_radius(hd);
  real_prec kernel_reach = SPH_kernel_radius(particle_kernel_type,
                                             particle_kernel_h);

  int reach1 = static_cast<int>(kernel_reach/d1) + 1;
  int reach2 = static_cast<int>(kernel_reach/d2) + 1;
  int reach3 = static_cast<int>(kernel_reach/d3) + 1;
  real_prec kernel_reach_sq = kernel_reach * kernel_reach;

  int out = 0;

  for (int i1 = -reach1; i1 <= reach1; ++i1)
    for (int i2 = -reach2; i2 <= reach2; ++i2)
      for (int i3 = -reach3; i3 <= reach3; ++i3) {
        real_prec dx = (abs(static_cast<real_prec>(i1)) - 0.5)*d1;
        real_prec dy = (abs(static_cast<real_prec>(i2)) - 0.5)*d2;
        real_prec dz = (abs(static_cast<real_prec>(i3)) - 0.5)*d3;
        real_prec r_sq = dx*dx + dy*dy + dz*dz;  // squared for efficiency

        if (r_sq <= kernel_reach_sq) {
          ++out;
        }
      }

  return(out);
}


// std::vector<int> SPH_kernel_3D_cells(struct HAMIL_DATA *hd)
void SPH_kernel_3D_cells(struct HAMIL_DATA *hd, vector<int> &out_i,
                         vector<int> &out_j, vector<int> &out_k) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  SPH_kernel_3D_cells(n->particle_kernel, n->particle_kernel_h,
                      n->d1, n->d2, n->d3, out_i, out_j, out_k);
}

void SPH_kernel_3D_cells(int particle_kernel_type, real_prec particle_kernel_h,
                         real_prec d1, real_prec d2, real_prec d3,
                         vector<int> &out_i,
                         vector<int> &out_j, vector<int> &out_k) {
  // First determine the reach of the kernel, which determines how many cells
  // every particle needs to loop over to check for contributions. This is
  // based on the kernel radius (e.g. 2*kernel_h for SPH with splines).
  real_prec kernel_reach = SPH_kernel_radius(particle_kernel_type,
                                             particle_kernel_h);
  int reach1 = static_cast<int>(kernel_reach/d1) + 1;
  int reach2 = static_cast<int>(kernel_reach/d2) + 1;
  int reach3 = static_cast<int>(kernel_reach/d3) + 1;
  // Determine sphere to loop over (+ 0.5dx, because we need to accomodate all
  // possible positions in a cell; note that this is implemented as -0.5dx
  // instead of just increasing the kernel reach by 0.5dx, because in principle
  // dx can be different in each direction, so that can only be incorporated if
  // done in each direction independently).
  real_prec kernel_reach_sq = kernel_reach * kernel_reach;
  // std::vector<int> kernel_cells;
  // int ix_cell = 0;

  for(int i1 = -reach1; i1 <= reach1; ++i1)
    for(int i2 = -reach2; i2 <= reach2; ++i2)
      for(int i3 = -reach3; i3 <= reach3; ++i3)
      {
        real_prec dx = (abs(static_cast<real_prec>(i1)) - 0.5)*d1; // abs, otherwise the "- 0.5" doesn't work for negative indices
        real_prec dy = (abs(static_cast<real_prec>(i2)) - 0.5)*d2;
        real_prec dz = (abs(static_cast<real_prec>(i3)) - 0.5)*d3;
        real_prec r_sq = dx*dx + dy*dy + dz*dz; // squared for efficiency

        if (r_sq <= kernel_reach_sq)
        {
          out_i.push_back(i1);
          out_j.push_back(i2);
          out_k.push_back(i3);
          // ++ix_cell;
        }
      }

  // return(kernel_cells);
}



// twee mogelijkheden:
// 1. 3 vectors met i&j, k_begin en k_last (niet k_end, inclusive!)
// 2. 2 2D vectors met k_begin en k_end (vierkant) waar indices i en j zijn
// optie 1:
void SPH_kernel_3D_cells_hull_1(const vector<int> &i, const vector<int> &j, const vector<int> &k,
                                vector< pair<int, int> > &ij_out,
                                vector<int> &k_begin, vector<int> &k_last) {
  auto N_cells = i.size();

  pair<int, int> ij0(i[0], j[0]);
  ij_out.push_back(ij0);
  k_begin.push_back(k[0]);
  k_last.push_back(k[0]);

  for (unsigned int ix = 1; ix < N_cells; ++ix) {
    // cout << ix << endl;
    // cout << i[ix] << " " << j[ix] << endl;
    pair<int, int> ij(i[ix], j[ix]);
    auto ij_duplicate = find(ij_out.begin(), ij_out.end(), ij);
    if (ij_duplicate == ij_out.end()) {
      // cout << "nieuw" << endl << endl;
      ij_out.push_back(ij);
      k_begin.push_back(k[ix]);
      k_last.push_back(k[ix]);
    } else {
      // cout << "duplicaat" << endl << endl;
      ULONG ix_duplicate = static_cast<ULONG>(distance(ij_out.begin(), ij_duplicate));
      if (k_begin[ix_duplicate] > k[ix])
        k_begin[ix_duplicate] = k[ix];
      if (k_last[ix_duplicate] < k[ix])
        k_last[ix_duplicate] = k[ix];
    }
  }
}

// optie 2:
// void SPH_kernel_3D_cells_hull_2(vector<int> &i, vector<int> &j, vector<int> &k,
//                                 vector<int> &k_begin, vector<int> &k_last) {
  // auto N_cells = i.size();
  // for (int ix = 0; ix < N_cells; ++ix) {}
// }


// TODO: check changed signedness warnings in below code, but make sure the changes don't make it slower, since this is an arduously optimized piece of code! Original version of the function put in comments below.

// TODO: make sure somehow (using new type?) that ix/y/z are never larger than INT_MAX! Otherwise casts below can go wrong.
// TODO: same goes for index_xy_part below, which should not be larger than LONG_MAX.
/*
 * N2pad and N3pad are ULONG because we need that to promote to ULONG in the calculation of index_xy_part, but in fact
 * they are not expected to be higher than uint_max! I.e. conceptually they are unsigned int, not ulong.
 */
void _likelihood_calc_V_SPH_kernel_loop_h_units(ULONG N2pad, ULONG N3pad, real_prec d1_h, real_prec d2_h, real_prec d3_h,
                                                vector<pair<int, int> > &ij, vector<int> &k_begin, vector<int> &k_last,
                                                unsigned ix, unsigned iy, unsigned iz, real_prec dpcx_h, real_prec dpcy_h,
                                                real_prec dpcz_h, const real_prec *part_like_padded, unsigned padding,
                                                real_prec grad_SPH_kernel_norm, real_prec &_out_x_j, real_prec &_out_y_j,
                                                real_prec &_out_z_j) {
  real_prec out_x_j = 0., out_y_j = 0., out_z_j = 0.; // to avoid having to add to out_#[j]'s inside the loop, which is expensive, because these arrays are not at all contiguous => factor 0.82 timesaving
  unsigned ix_pad = ix + padding;
  unsigned iy_pad = iy + padding;
  unsigned iz_pad = iz + padding;
  // Note: don't use unsigned int as index below! Slower than signed. See:
  // http://stackoverflow.com/a/2044021/1199693.
  // Can't seem to be able to remove warning (unsigned vs signed comparison).
  // Tried converting the ij.size() to long and int, but is again slower...
  for (vector< pair<int, int> >::size_type ij_ix = 0;
       ij_ix < ij.size(); ++ij_ix) {
    int i1 = ij[ij_ix].first;
    int i2 = ij[ij_ix].second;
    auto kx = static_cast<unsigned>(static_cast<int>(ix_pad) + i1);
    auto ky = static_cast<unsigned>(static_cast<int>(iy_pad) + i2);
    // Cell position (relative to the central cell):
    real_prec diff_x_h = dpcx_h - static_cast<real_prec>(i1)*d1_h;
    real_prec diff_y_h = dpcy_h - static_cast<real_prec>(i2)*d2_h;
    // magic
    ULONG index_xy_part = N3pad*(ky + N2pad*kx);  // N2/3pad will promote kx/y to ULONG
    int kz_begin = k_begin[ij_ix];
    int kz_last = k_last[ij_ix];
    auto index_begin = static_cast<ULONG>(kz_begin + static_cast<long>(index_xy_part + iz_pad));
    ULONG index_end = index_begin + static_cast<ULONG>(kz_last - kz_begin);
    real_prec diff_z_h = dpcz_h - static_cast<real_prec>(kz_begin) * d3_h;
    // the actual loop
//    for (int i3 = kz_begin; i3 <= kz_last; ++i3) {
    for (ULONG index = index_begin; index <= index_end; ++index) {
      real_prec common_part = part_like_padded[index];

      real_prec grad_kernel_x, grad_kernel_y, grad_kernel_z;
      grad_SPH_kernel_3D_h_units(diff_x_h, diff_y_h, diff_z_h, grad_SPH_kernel_norm,
                                 grad_kernel_x, grad_kernel_y, grad_kernel_z);

      out_x_j += common_part * grad_kernel_x;
      out_y_j += common_part * grad_kernel_y;
      out_z_j += common_part * grad_kernel_z;

      // for next iteration:
//      ++index;
      diff_z_h -= d3_h;
    }
  }
  _out_x_j = out_x_j;
  _out_y_j = out_y_j;
  _out_z_j = out_z_j;
}


// NOTE:
// Below the original version of the above function. The original version was optimized by hand.
// However, the compiler complained about implicit conversions, which are now removed from the above
// version. TODO: check at some point which version is faster!
//void _likelihood_calc_V_SPH_kernel_loop_h_units(int N2pad, int N3pad, real_prec d1_h, real_prec d2_h, real_prec d3_h,
//                                               vector<pair<int, int> > &ij, vector<int> &k_begin, vector<int> &k_last,
//                                               int ix, int iy, int iz, real_prec dpcx_h, real_prec dpcy_h,
//                                               real_prec dpcz_h, real_prec *part_like_padded, int padding,
//                                               real_prec grad_SPH_kernel_norm, real_prec &_out_x_j, real_prec &_out_y_j,
//                                               real_prec &_out_z_j) {
//  real_prec out_x_j = 0., out_y_j = 0., out_z_j = 0.; // to avoid having to add to out_#[j]'s inside the loop, which is expensive, because these arrays are not at all contiguous => factor 0.82 timesaving
//  int ix_pad = ix + padding;
//  int iy_pad = iy + padding;
//  int iz_pad = iz + padding;
//  // Note: don't use unsigned int as index below! Slower than signed. See:
//  // http://stackoverflow.com/a/2044021/1199693.
//  // Can't seem to be able to remove warning (unsigned vs signed comparison).
//  // Tried converting the ij.size() to long and int, but is again slower...
//  for (vector< pair<int, int> >::size_type ij_ix = 0;
//       ij_ix < ij.size(); ++ij_ix) {
//    int i1 = ij[ij_ix].first;
//    int i2 = ij[ij_ix].second;
//    int kx = ix_pad + i1;
//    int ky = iy_pad + i2;
//    // Cell position (relative to the central cell):
//    real_prec diff_x_h = dpcx_h - static_cast<real_prec>(i1)*d1_h;
//    real_prec diff_y_h = dpcy_h - static_cast<real_prec>(i2)*d2_h;
//    // magic
//    ULONG index_xy_part = N3pad*(ky + N2pad*kx);
//    int kz_begin = k_begin[ij_ix];
//    int kz_last = k_last[ij_ix];
//    ULONG index = kz_begin + index_xy_part + iz_pad;
//    real_prec diff_z_h = dpcz_h - static_cast<real_prec>(kz_begin) * d3_h;
//    // the actual loop
//    for (int i3 = kz_begin; i3 <= kz_last; ++i3) {
//      real_prec common_part = part_like_padded[index];
//
//      real_prec grad_kernel_x, grad_kernel_y, grad_kernel_z;
//      // grad_SPH_kernel_3D(diff_x_h, diff_y_h, diff_z_h, particle_kernel_h,
//      //                    h_sq, h_inv, h_sq_inv, grad_SPH_kernel_norm,
//      //                    grad_kernel_x, grad_kernel_y, grad_kernel_z);
//      grad_SPH_kernel_3D_h_units(diff_x_h, diff_y_h, diff_z_h, grad_SPH_kernel_norm,
//                                 grad_kernel_x, grad_kernel_y, grad_kernel_z);
//
//      out_x_j += common_part * grad_kernel_x;
//      out_y_j += common_part * grad_kernel_y;
//      out_z_j += common_part * grad_kernel_z;
//
//      // for next iteration:
//      ++index;
//      diff_z_h -= d3_h;
//    }
//  }
//  _out_x_j = out_x_j;
//  _out_y_j = out_y_j;
//  _out_z_j = out_z_j;
//}

void pad_array_pacman(real_prec *input, unsigned int N1_in, real_prec *out,
                      unsigned int padding) {
  unsigned N1_out = N1_in + 2*padding;
  for (unsigned io = 0; io < N1_out; ++io)
    for (unsigned jo = 0; jo < N1_out; ++jo)
      for (unsigned ko = 0; ko < N1_out; ++ko) {
        ULONG ix_out = ko + N1_out*(jo + static_cast<ULONG>(N1_out)*io);
        unsigned ii = static_cast<unsigned>(static_cast<int>(io + N1_in) - static_cast<int>(padding)) % N1_in;
        unsigned ji = static_cast<unsigned>(static_cast<int>(jo + N1_in) - static_cast<int>(padding)) % N1_in;
        unsigned ki = static_cast<unsigned>(static_cast<int>(ko + N1_in) - static_cast<int>(padding)) % N1_in;
        ULONG ix_in = ki + N1_in*(ji + static_cast<ULONG>(N1_in)*ii);
        out[ix_out] = input[ix_in];
      }
}



// TODO: COMPARE THE OUTPUT OF THIS FUNCTION FROM AFTER 6 JUNE 2017 TO VERSION FROM BEFORE!!!
// TODO: COMPARE THE OUTPUT OF THIS FUNCTION FROM AFTER 6 JUNE 2017 TO VERSION FROM BEFORE!!!
// TODO: COMPARE THE OUTPUT OF THIS FUNCTION FROM AFTER 6 JUNE 2017 TO VERSION FROM BEFORE!!!
// TODO: COMPARE THE OUTPUT OF THIS FUNCTION FROM AFTER 6 JUNE 2017 TO VERSION FROM BEFORE!!!
// TODO: COMPARE THE OUTPUT OF THIS FUNCTION FROM AFTER 6 JUNE 2017 TO VERSION FROM BEFORE!!!
// TODO: COMPARE THE OUTPUT OF THIS FUNCTION FROM AFTER 6 JUNE 2017 TO VERSION FROM BEFORE!!!
// TODO: COMPARE THE OUTPUT OF THIS FUNCTION FROM AFTER 6 JUNE 2017 TO VERSION FROM BEFORE!!!
// TODO: also pad_array_pacman

void likelihood_calc_V_SPH(struct HAMIL_DATA *hd, real_prec *part_like, real_prec *posx, real_prec *posy, real_prec *posz, real_prec *out_x, real_prec *out_y, real_prec *out_z)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;

  // int N_cells = SPH_kernel_3D_cells_count(hd);
  // // const int kernel_cells_i[N_cells], kernel_cells_j[N_cells], kernel_cells_k[N_cells];
  // vector<int> kernel_cells_i, kernel_cells_j, kernel_cells_k;
  // SPH_kernel_3D_cells(hd, kernel_cells_i, kernel_cells_j, kernel_cells_k);
  const vector<int> kernel_cells_i = hd->kernel_cells_i;
  const vector<int> kernel_cells_j = hd->kernel_cells_j;
  const vector<int> kernel_cells_k = hd->kernel_cells_k;
  // const int N_cells = hd->N_cells;

  std::vector< std::pair<int, int> > ij;
  std::vector<int> k_begin, k_last;
  SPH_kernel_3D_cells_hull_1(kernel_cells_i, kernel_cells_j, kernel_cells_k,
                             ij, k_begin, k_last);

  unsigned padding = static_cast<unsigned>(*max_element(kernel_cells_i.begin(), kernel_cells_i.end()));
  fftw_array<real_prec> part_like_padded((n->N1+2*padding)*(n->N2+2*padding)*
                                         (n->N3+2*padding));
  pad_array_pacman(part_like, n->N1, part_like_padded, padding);

  // Normalization given that we want the integral over density to be V (left-
  // hand side) and the integral over all particle kernels to be N (right-hand
  // side).
  // Basically, this is the mass of all particles given that we want the
  // density to be rho_c.
  real_prec normalize = hd->rho_c * n->L1*n->L2*n->L3/static_cast<real_prec>(n->N1*n->N2*n->N3);

  real_prec f1;
  bool rsd_model = hd->rsd_model;
  bool planepar = n->planepar;
  if (rsd_model)
    f1 = fgrow(hd->ascale, hd->OM, hd->OL, 1); // first order growth factor (Zel'dovich only!)

  real_prec h = n->particle_kernel_h;
  real_prec h_sq = h*h;
  real_prec h_inv = 1. / h;
  // real_prec h_sq_inv = 1. / h_sq;
  real_prec grad_SPH_kernel_norm = 1. / (M_PI*h_sq*h_sq);
  unsigned N2 = n->N2, N3 = n->N3;
  unsigned N3pad = N3 + 2*padding;
  unsigned N2pad = N2 + 2*padding;
  real_prec d1 = n->d1, d2 = n->d2, d3 = n->d3;
  real_prec d1_h = d1 * h_inv, d2_h = d2 * h_inv, d3_h = d3 * h_inv;

  // Initialize this thing outside the for-loop for performance
  // std::vector<real_prec> grad_kernel(3);
  // Firstprivate below initializes grad_kernel in each thread to a copy
  // of the value outside the loop. Private only initializes an empty
  // vector, which then does not have the correct length.
  #ifdef MULTITHREAD
  // #pragma omp parallel for firstprivate(grad_kernel)
  #pragma omp parallel for
  #endif // MULTITHREAD 
  for (ULONG j = 0; j < n->N; j++)
  {
    // Load particle position
    real_prec px = posx[j], py = posy[j], pz=posz[j];
    // Determine central cell index where particle resides
    // TODO: this calculation will go wrong when posx,y,z are not already put
    // in periodic box coordinates; might be negative then.
    // Addition 6 Jun 2017: TODO: guarantee that positions are in periodic box coords by making special type for that.
    auto ix = static_cast<int>(px/d1);
    auto iy = static_cast<int>(py/d2);
    auto iz = static_cast<int>(pz/d3);
    // Central cell position:
    // real_prec ccx = (static_cast<real_prec>(ix) + 0.5)*d1;
    // real_prec ccy = (static_cast<real_prec>(iy) + 0.5)*d2;
    // real_prec ccz = (static_cast<real_prec>(iz) + 0.5)*d3;
    real_prec ccx_h = (static_cast<real_prec>(ix) + 0.5)*d1_h;
    real_prec ccy_h = (static_cast<real_prec>(iy) + 0.5)*d2_h;
    real_prec ccz_h = (static_cast<real_prec>(iz) + 0.5)*d3_h;
    // Optimization: precalculate diff. pos and cell
    // real_prec dpcx = px - ccx;
    // real_prec dpcy = py - ccy;
    // real_prec dpcz = pz - ccz;
    real_prec dpcx_h = px * h_inv - ccx_h;
    real_prec dpcy_h = py * h_inv - ccy_h;
    real_prec dpcz_h = pz * h_inv - ccz_h;

    real_prec out_x_j, out_y_j, out_z_j;

    // 6 Jun 2017: don't convert ix/y/z to unsigned above, since conversion between uint and double (for ccx_h) is
    //             slower than between int and double!
    _likelihood_calc_V_SPH_kernel_loop_h_units(N2pad, N3pad, d1_h, d2_h, d3_h, ij, k_begin, k_last,
                                               static_cast<unsigned>(ix), static_cast<unsigned>(iy),
                                               static_cast<unsigned>(iz), dpcx_h, dpcy_h, dpcz_h, part_like_padded,
                                               padding, grad_SPH_kernel_norm, out_x_j, out_y_j, out_z_j);

    out_x[j] = normalize * out_x_j;
    out_y[j] = normalize * out_y_j;
    out_z[j] = normalize * out_z_j;

    if (rsd_model) {
      if (!planepar) {
        throw BarcodeException("Non-plane-parallel RSD model is not yet implemented in calc_V! Use planepar = true.");
      } else {
        out_z[j] += f1 * out_z[j];
      }
    }
  }
}




// new Fourier based version with TSC interpolation
void likelihood_calc_V_SPH_fourier_TSC(struct HAMIL_DATA *hd, real_prec *part_like, real_prec *out_x, real_prec *out_y,
                                       real_prec *out_z) {
  // Works fastest if part_like == n->R2Cplan->R.
  // NOTE: n->R2Cplan->R is used in calc_h_SPH for posx as well, so part_like
  //       will be destroyed if it is also n->R2Cplan->R!

  // 0. initialization
  struct HAMIL_NUMERICAL *n = hd->numerical;
  real_prec h = n->particle_kernel_h;

  // kernel normalization
  real_prec norm_kernel = 24./(h*h*h);

  // Normalization given that we want the integral over density to be V (left-
  // hand side) and the integral over all particle kernels to be N (right-hand
  // side).
  // Basically, this is the mass of all particles given that we want the
  // density to be rho_c.
  real_prec norm_density = hd->rho_c * n->L1*n->L2*n->L3/static_cast<real_prec>(n->N1*n->N2*n->N3);

  real_prec norm = norm_kernel * norm_density;

  // 1. calculate convolution of part_like and kernel
  // 1.a. calculate fourier transform of part_like
  // fftw_array<complex_prec> part_like_F(n->Nhalf), conv_x_F(n->Nhalf),
  //                          conv_y_F(n->Nhalf);
  // fftw_array<complex_prec> conv_x_F(n->Nhalf), conv_y_F(n->Nhalf);
  fftw_array<complex_prec> conv_y_F(n->Nhalf);  // use n->C2Rplan->C as conv_x_F

  // fftR2C(n->N1, n->N2, n->N3, part_like, part_like_F);
  fftR2Cplanned(part_like, n->R2Cplan->C, n->R2Cplan);  // use n->R2Cplan->C as part_like_F

  real_prec L1 = n->L1, L2 = n->L2, L3 = n->L3;
  unsigned N1 = n->N1, N2 = n->N2, N3 = n->N3;
  unsigned N3half = N3/2 + 1;

  #ifdef MULTITHREAD
  #pragma omp parallel for //schedule(static,chunk) //collapse(3)
  #endif // MULTITHREAD
  for (unsigned i = 0; i < N1; ++i) {
    real_prec kx = calc_kx(i, L1, N1);
    for (unsigned j = 0; j < N2; ++j) {
      real_prec ky = calc_ky(j, L2, N2);
      for (unsigned k = 0; k < N3half; ++k) {
        real_prec kz = calc_kz(k, L3, N3);
        real_prec k_sq = kx*kx + ky*ky + kz*kz;
        // 1.b. calculate fourier transform of SPH kernel (it's real, so don't
        //      need complex number)
        real_prec SPH_kernel_F;
        if (k_sq == 0.) {
          SPH_kernel_F = 1./(h*h*h);
        } else {
          real_prec kk = sqrt(k_sq);
          real_prec ksink = kk * sin(kk);

          SPH_kernel_F = norm * (3 + cos(2*kk) - ksink + cos(kk) * (ksink - 4))
                         / (k_sq*k_sq*k_sq);
        }

        // 1.c. multiply fourier transforms of part_like and SPH kernel and take
        //      derivative (multiply by ik)
        // Also multiply with h to correct for wrong coordinate in derivative
        // (x(q_i) instead of q=(x(q_i)-x)/h which is the Fourier dual of k)
        ULONG ix = k + N3half*(j + static_cast<ULONG>(N2)*i);
        // x
        // re(conv_x_F[ix]) = h * kx * -im(part_like_F[ix]) * SPH_kernel_F;
        // im(conv_x_F[ix]) = h * kx *  re(part_like_F[ix]) * SPH_kernel_F;
        re(n->C2Rplan->C[ix]) = h * kx * -im(n->R2Cplan->C[ix]) * SPH_kernel_F;
        im(n->C2Rplan->C[ix]) = h * kx *  re(n->R2Cplan->C[ix]) * SPH_kernel_F;
        // y
        // re(conv_y_F[ix]) = h * ky * -im(part_like_F[ix]) * SPH_kernel_F;
        // im(conv_y_F[ix]) = h * ky *  re(part_like_F[ix]) * SPH_kernel_F;
        re(conv_y_F[ix]) = h * ky * -im(n->R2Cplan->C[ix]) * SPH_kernel_F;
        im(conv_y_F[ix]) = h * ky *  re(n->R2Cplan->C[ix]) * SPH_kernel_F;
        // z -> part_like_C itself, to save memory
        // real_prec dummy = re(part_like_F[ix]);
        // re(part_like_F[ix]) = h * kz * -im(part_like_F[ix]) * SPH_kernel_F;
        real_prec dummy = re(n->R2Cplan->C[ix]);
        re(n->R2Cplan->C[ix]) = h * kz * -im(n->R2Cplan->C[ix]) * SPH_kernel_F;
        im(n->R2Cplan->C[ix]) = h * kz * dummy                  * SPH_kernel_F;
      }
    }
  }

  ULONG N_part = n->N;

  // 1.d. transform back to real space => convolution done
  // fftC2R(n->N1, n->N2, n->N3, conv_x_F, out_z);  // out_z is dummy
  fftC2Rplanned(n->C2Rplan->C, n->C2Rplan->R, n->C2Rplan);  // use n->C2Rplan->R as dummy
  // 2. interpolate convolution of part_like and kernel to particle positions
  // interpolate_TSC(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1, n->d2, n->d3,
  //                 hd->posx, hd->posy, hd->posz, out_z, N_part, out_x);
  interpolate_TSC(n->N1, n->N2, n->N3, n->d1, n->d2, n->d3, hd->posx, hd->posy, hd->posz, n->C2Rplan->R, N_part, out_x);

  // fftC2R(n->N1, n->N2, n->N3, conv_y_F, out_z);  // again out_z dummy
  fftC2Rplanned(conv_y_F, n->C2Rplan->R, n->C2Rplan);  // again n->C2Rplan->R as dummy
  // interpolate_TSC(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1, n->d2, n->d3,
  //                 hd->posx, hd->posy, hd->posz, out_z, N_part, out_y);
  interpolate_TSC(n->N1, n->N2, n->N3, n->d1, n->d2, n->d3, hd->posx, hd->posy, hd->posz, n->C2Rplan->R, N_part, out_y);

  // fftw_array<real_prec> conv_x(n->N), conv_y(n->N);  // multi version
  // fftC2R(n->N1, n->N2, n->N3, conv_x_F, conv_x);     // multi version
  // fftC2R(n->N1, n->N2, n->N3, conv_y_F, conv_y);     // multi version
  // here part_like_F is still dummy and we also use part_like as dummy
  // fftC2R(n->N1, n->N2, n->N3, part_like_F, part_like);  // ALSO for multi!
  // fftC2R(n->N1, n->N2, n->N3, n->R2Cplan->C, part_like);  // ALSO for multi!
  fftC2Rplanned(n->R2Cplan->C, n->C2Rplan->R, n->C2Rplan);  // again n->C2Rplan->R as dummy  // ALSO for multi!
  // interpolate_TSC(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1, n->d2, n->d3,
  //                 hd->posx, hd->posy, hd->posz, part_like, N_part, out_z);
  interpolate_TSC(n->N1, n->N2, n->N3, n->d1, n->d2, n->d3, hd->posx, hd->posy, hd->posz, n->C2Rplan->R, N_part, out_z);

  // real_prec *input_fields[3]  = {conv_x, conv_y, part_like}; // multi version
  // real_prec *output_fields[3] = {out_x, out_y, out_z};       // multi version
  // interpolate_TSC_multi(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1,
  //                       n->d2, n->d3, hd->posx, hd->posy, hd->posz, N_part,
  //                       input_fields, output_fields, 3);     // multi version

  // RSD stuff
  if (hd->rsd_model) {
    if (!n->planepar) {
      throw BarcodeException("Non-plane-parallel RSD model is not yet "
                             "implemented in calc_V! Use planepar = true.");
    } else {
      // first order growth factor (Zel'dovich only!)
      real_prec f1 = fgrow(hd->ascale, hd->OM, hd->OL, 1);
      #ifdef MULTITHREAD
      #pragma omp parallel for
      #endif  // MULTITHREAD
      for (unsigned int ix = 0; ix < n->N; ++ix) {
        out_z[ix] += f1 * out_z[ix];
      }
    }
  }

}




void likelihood_calc_h_SPH(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out) {
  // fastest if out == n->C2Rplan->R
  struct HAMIL_NUMERICAL *n = hd->numerical;

  if (!(n->mk == 3)) {
    throw BarcodeException("Must use SPH mass kernel (masskernel = 3) when "
                           "using likelihood_calc_h_SPH (calc_h = 2 or 3)!");
  }

  bool rfft = true;

  // fftw_array<real_prec> part_like(n->N), V_x(n->N), V_y(n->N), V_z(n->N);
  // fftw_array<real_prec> V_x(n->N), V_y(n->N), V_z(n->N);
  fftw_array<real_prec> V_y(n->N), V_z(n->N);
  // hd->partial_f_delta_x_log_like(hd, deltaX, part_like);
  hd->partial_f_delta_x_log_like(hd, deltaX, n->R2Cplan->R);

  // The first fft below this switch works fastest if V_x == n->R2Cplan->R.
  // This is allowed in likelihood_calc_V_SPH_fourier_TSC, but it overwrites
  // part_like, because we also use n->R2Cplan->R for that!
  switch (n->calc_h) {
    case 2:
      // likelihood_calc_V_SPH(hd, part_like, hd->posx, hd->posy, hd->posz,
      //                       V_x, V_y, V_z);
      likelihood_calc_V_SPH(hd, n->R2Cplan->R, hd->posx, hd->posy, hd->posz,
                            n->R2Cplan->R, V_y, V_z);
      break;
    case 3:
      // likelihood_calc_V_SPH_fourier_TSC(hd, part_like, hd->posx, hd->posy,
      //                                   hd->posz, V_x, V_y, V_z);
      likelihood_calc_V_SPH_fourier_TSC(hd, n->R2Cplan->R, n->R2Cplan->R, V_y, V_z);
      break;
  }

  // fftw_array<complex_prec> outC(n->Nhalf), dummyC(n->Nhalf);

  // use outC as dummy in first run, saves one addition and one fillZero
  // fftR2C(n->N1, n->N2, n->N3, V_x, outC);
  fftR2Cplanned(n->R2Cplan->R, n->C2Rplan->C, n->R2Cplan);  // use n->C2Rplan->C as outC

  // grad_inv_lap_FS(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, outC, outC, 1, rfft); // multiply by -ik_i/k^2
  grad_inv_lap_FS(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, n->C2Rplan->C, n->C2Rplan->C, 1, rfft); // multiply by -ik_i/k^2

  // fftR2C(n->N1, n->N2, n->N3, V_y, dummyC);
  fftR2Cplanned(V_y, n->R2Cplan->C, n->R2Cplan);  // use n->R2Cplan->C as dummyC
  // grad_inv_lap_FS(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, dummyC, dummyC, 2, rfft); // multiply by -ik_i/k^2
  grad_inv_lap_FS(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, n->R2Cplan->C, n->R2Cplan->C, 2, rfft); // multiply by -ik_i/k^2
  // add_to_array(dummyC, outC, n->Nhalf); // sum the resulting h(k) term to outC
  add_to_array(n->R2Cplan->C, n->C2Rplan->C, n->Nhalf); // sum the resulting h(k) term to outC

  // fftR2C(n->N1, n->N2, n->N3, V_z, dummyC);
  fftR2Cplanned(V_z, n->R2Cplan->C, n->R2Cplan);  // use n->R2Cplan->C as dummyC
  // grad_inv_lap_FS(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, dummyC, dummyC, 3, rfft); // multiply by -ik_i/k^2
  grad_inv_lap_FS(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, n->R2Cplan->C, n->R2Cplan->C, 3, rfft); // multiply by -ik_i/k^2
  // add_to_array(dummyC, outC, n->Nhalf); // sum the resulting h(k) term to outC
  add_to_array(n->R2Cplan->C, n->C2Rplan->C, n->Nhalf); // sum the resulting h(k) term to outC

  // transform h(k) to real space -> h(q)
  // fftC2R(n->N1, n->N2, n->N3, outC, out);
  fftC2Rplanned(n->C2Rplan->C, out, n->C2Rplan);  // use n->R2Cplan->C as dummyC
}
// end SPH method stuff //


// Model: Gaussian likelihood (or any other that is dependent on delta_q only through delta_x) //
void likelihood_grad_log_like(struct HAMIL_DATA *hd, real_prec *delta, real_prec *out)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  // fftw_array<real_prec> dummy(n->N);  // use n->C2Rplan->R as dummy
  
  // single out growing mode (Peebles) -> n->C2Rplan->R (dummy)
  if (n->deltaQ_factor != 1.) {
    multiply_factor_array(n->deltaQ_factor, delta, n->C2Rplan->R, n->N);
  } else {
    copyArray(delta, n->C2Rplan->R, n->N);
  }

  // calculate deltaX and pos
  {
    unsigned facL=1;
    bool reggrid=true;
    gsl_rng *seed = nullptr; // empty: reggrid is true anyway
    real_prec kernel_scale = SPH_kernel_scale(hd);
    if (hd->rsd_model)
      Lag2Eul_rsd_zeldovich(n->C2Rplan->R, hd->deltaX, hd->posx, hd->posy, hd->posz, n->N1, n->N2, n->N3, n->L1, n->L2,
                            n->L3, n->d1, n->d2, n->d3, n->min1, n->min2, n->min3, hd->D1, hd->ascale, hd->OM, hd->OL,
                            n->mk, facL,
                            reggrid, seed, kernel_scale, n->xobs, n->yobs, n->zobs, n->planepar, n->periodic,
                            n->R2Cplan, n->C2Rplan);
    else
      Lag2Eul(n->C2Rplan->R, hd->deltaX, hd->posx, hd->posy, hd->posz, n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1,
              n->d2, n->d3, n->min1, n->min2, n->min3, hd->D1, hd->D2, hd->ascale, hd->OM, hd->OL, hd->sfmodel, n->mk, n->kth,
              facL,
              reggrid, seed, "", kernel_scale, n->R2Cplan, n->C2Rplan);
  }

  // h -> n->C2Rplan->R (dummy)
  switch (n->calc_h)
  {
    case 0:
      likelihood_calc_h(hd, hd->deltaX, n->C2Rplan->R);
      break;
    case 1:
      hd->partial_f_delta_x_log_like(hd, hd->deltaX, n->C2Rplan->R);
      break;
    case 2:
    case 3:
      likelihood_calc_h_SPH(hd, hd->deltaX, n->C2Rplan->R);
      break;
  }

  // now for some normalization terms
  real_prec norm = 1.;

  // ************************* WARNING **************************
  // No longer using heuristic factor, it was caused by a mismatch between mass
  // kernel used for density estimation and the one used for calc_h! As long as
  // the two match, the correspondence between fin. diff. and calculated h is
  // nearly perfect!
  // ************************* WARNING **************************

  // heuristically determined correction factor
  // EGP: I have no idea where this came from, but it seems necessary. See tests
  //      in ipynb "unittest - calc_V_fourier"
  //      NOTE: I'm using n->d1 now; this is of course only correct if all d's
  //            are the same!
  // real_prec heuristic_correction = 1.;
  // switch (n->calc_h)
  // {
  //   case 2:
  //     heuristic_correction = 0.52/0.33;
  //     break;
  //   case 3:
  //     if (n->d1 != n->d2 || n->d2 != n->d3) {
  //       throw BarcodeException("likelihood_grad_log_like: d1, d2 and d3 are not"
  //                              " equal, so heuristic_correction is unknown! "
  //                              "Aborting.");
  //     } else {
  //       heuristic_correction = 0.52/0.31/pow(n->d1, 1.04);
  //     }
  //     break;
  // }
  // norm *= heuristic_correction;

  // ************************* WARNING ************************** (see above)

  // EGP: for Zel'dovich, -h is the entire term. For other models more has to be done.
  // FIXME: hier moet nog een factor D1 bij als het signaal alleen delta^(1) is
  real_prec zeldovich_norm = -1.;  // EGP: Zel'dovich: -gradLogLike = -h
  norm *= zeldovich_norm;

  // TODO: check dit theoretisch!
  norm *= n->deltaQ_factor;

  if (n->correct_delta) {
    norm *= hd->D1;
  }

  multiply_factor_array(norm, n->C2Rplan->R, out, n->N);
}
void likelihood_grad_log_like_old(struct HAMIL_DATA *hd, real_prec *delta, real_prec *out)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  fftw_array<real_prec> dummy_h(n->N);
  
  // calculate deltaX and pos
  {
    unsigned facL=1;
    bool reggrid=true;
    gsl_rng *seed = nullptr; // empty: reggrid is true anyway
    real_prec kernel_scale = SPH_kernel_scale(hd);
    Lag2Eul(delta, hd->deltaX, hd->posx, hd->posy, hd->posz, n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1, n->d2,
            n->d3, n->min1, n->min2, n->min3, hd->D1, hd->D2, hd->ascale, hd->OM, hd->OL, hd->sfmodel, n->mk, n->kth, facL,
            reggrid, seed, "", kernel_scale, n->R2Cplan, n->C2Rplan);
  }

  likelihood_calc_h_old(hd, hd->deltaX, dummy_h);

  // EGP: for Zel'dovich, -h is the entire term. For other models more has to be done.
  multiply_factor_array(-1., dummy_h, out, n->N); // EGP: Zel'dovich: -gradLogLike = -h
}


//void prior_gaussian_grad_log_prior(struct HAMIL_DATA *hd, real_prec *signal, real_prec *out)
//{
  //// dPsiP/ds = S^-1 * s
  //convolveInvCorrFuncWithSignal(hd, signal, out, hd->signal_PS);
//}
void prior_gaussian_grad_log_prior(struct HAMIL_DATA *hd, real_prec *signal, real_prec *out)
{
  //struct HAMIL_NUMERICAL *n = hd->numerical;
  //fftw_array<real_prec> delta_growing(n->N);

  //// single out growing mode (Peebles) -> delta_growing
  //multiply_factor_array(n->deltaQ_factor, signal, delta_growing, n->N);
  // dPsiP/ds = S^-1 * s
  //convolveInvCorrFuncWithSignal(hd, delta_growing, out, hd->signal_PS);
  convolveInvCorrFuncWithSignal(hd, signal, out, hd->signal_PS);
}

//real_prec prior_gaussian_log_prior(struct HAMIL_DATA *hd, real_prec *signal)
//{
  //struct HAMIL_NUMERICAL *n = hd->numerical;
  //fftw_array<real_prec> dummy(n->N);

  //convolveInvCorrFuncWithSignal(hd, signal, dummy, hd->signal_PS);

  //// Psi_prior = 1/2 [s_k] * IFT[ 1/P*FT[s_k] ]

  //real_prec psi_prior=0.;
//#ifdef MULTITHREAD
//#pragma omp parallel for reduction(+:psi_prior)
//#endif // MULTITHREAD 
  //for(ULONG i=0;i<n->N;i++)
    //psi_prior += num_0_5 * signal[i] * dummy[i];

  //return(psi_prior);
//}
real_prec prior_gaussian_log_prior(struct HAMIL_DATA *hd, real_prec *signal)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  fftw_array<real_prec> dummy(n->N);
  //fftw_array<real_prec> dummy(n->N), delta_growing(n->N);

  //// single out growing mode (Peebles) -> delta_growing
  //multiply_factor_array(n->deltaQ_factor, signal, delta_growing, n->N);

  //convolveInvCorrFuncWithSignal(hd, delta_growing, dummy, hd->signal_PS);
  convolveInvCorrFuncWithSignal(hd, signal, dummy, hd->signal_PS);

  // Psi_prior = 1/2 [s_k] * IFT[ 1/P*FT[s_k] ]

  real_prec psi_prior=0.;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:psi_prior)
#endif // MULTITHREAD 
  for(ULONG i=0;i<n->N;i++)
    psi_prior += num_0_5 * signal[i] * dummy[i];
    //psi_prior += num_0_5 * delta_growing[i] * dummy[i];

  return(psi_prior);
}
