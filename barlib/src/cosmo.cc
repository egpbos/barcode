/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include "define_opt.h"

#include <cmath>
#include <fstream>

#include <gsl/gsl_integration.h>


real_prec E_Hubble_z(real_prec z, real_prec Omega_M, real_prec Omega_L)
{
  real_prec z_fac = 1+z;
  real_prec Omega_K = num_1 - Omega_M - Omega_L;
  real_prec E_Hubble = sqrt(Omega_M*z_fac*z_fac*z_fac + Omega_K*z_fac*z_fac + Omega_L);
  return(E_Hubble);
}

real_prec E_Hubble_a(real_prec a, real_prec Omega_M, real_prec Omega_L)
{
  real_prec Omega_K = num_1 - Omega_M - Omega_L;
  real_prec E_Hubble = sqrt(Omega_M/(a*a*a) + Omega_K/(a*a) + Omega_L);
  return(E_Hubble);
}

// Integration stuff
struct my_f_params {gsl_real OmM; gsl_real OmK; gsl_real OmL;};

gsl_real
e_of_z (gsl_real z, void * p) {
   /*struct my_f_params * params 
     = (struct my_f_params *)p;*/
  auto * params 
     = static_cast<struct my_f_params *>( p ); // EGP: static_cast

   gsl_real OmM = (params->OmM);
   gsl_real OmK = (params->OmK);
   gsl_real OmL = (params->OmL);

   gsl_real e = sqrt(OmM*(1.+z)*(1.+z)*(1.+z)+OmK*(1.+z)*(1.+z)+OmL);

   return 1./e;
}

gsl_real
e_of_z_cube (real_prec z, void * p) {
   /*struct my_f_params * params 
     = (struct my_f_params *)p;*/
  auto * params 
     = static_cast<struct my_f_params *>( p ); // EGP: static_cast

   gsl_real OmM = (params->OmM);
   gsl_real OmK = (params->OmK);
   gsl_real OmL = (params->OmL);

   gsl_real e = sqrt(OmM*(1.+z)*(1.+z)*(1.+z)+OmK*(1.+z)*(1.+z)+OmL);

   return (1.+z)*(1.+z)*(1.+z)/(e*e*e);
}

gsl_real
growth_var (gsl_real z, void * p) {
   /*struct my_f_params * params 
     = (struct my_f_params *)p;*/
  auto * params 
     = static_cast<struct my_f_params *>( p ); // EGP: static_cast

   gsl_real OmM = (params->OmM);
   gsl_real OmK = (params->OmK);
   gsl_real OmL = (params->OmL);

   gsl_real e = sqrt(OmM*(1.+z)*(1.+z)*(1.+z)+OmK*(1.+z)*(1.+z)+OmL);

   return (1.+z)/(e*e*e);
} 
// END integration stuff

real_prec calc_dcom(real_prec scale_factor, real_prec Omega_M, real_prec Omega_L, real_prec hconst)
{
  real_prec Omega_C=static_cast<real_prec>(1.)-Omega_M-Omega_L;


  // Getting Cosmology
  auto H0=static_cast<real_prec>(100.*hconst*cgs_km/cgs_Mpc/cgs_sec);	
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  auto redshift=static_cast<gsl_real>(1./scale_factor-1.);

  gsl_real result, error, epsabs, epsrel;
  size_t neval;
  
  epsabs=1e-4;
  epsrel=1e-8;

  auto omega_m=static_cast<gsl_real>(Omega_M);
  auto omega_k=static_cast<gsl_real>(Omega_C);
  auto omega_q=static_cast<gsl_real>(Omega_L);

  gsl_function F;
  struct my_f_params params = {omega_m,omega_k,omega_q};

  F.function = &e_of_z;
  F.params = &params;

  //gsl_integration_qagiu (&F, redshift, 0, 1e-4, 1000,
  //		 w, &result, &error);
  
  gsl_integration_qng (&F, 0.,redshift, epsabs, epsrel, 
  		       &result, &error, &neval);
  
  gsl_integration_workspace_free (w);
  return static_cast<real_prec>(result*cgs_clight/H0);
}
///____________________________________________________________________

/* --- function D_growth growth factor --- */
real_prec D_growth(real_prec scale_factor, real_prec Omega_M, real_prec Omega_L, real_prec hconst)
{
  // Getting Cosmology
  real_prec Omega_C=static_cast<real_prec>(1.)-Omega_M-Omega_L;

  auto H0=static_cast<real_prec>(hconst*(cgs_km/cgs_Mpc/cgs_sec));
  auto H=static_cast<real_prec>(H0*sqrt(Omega_M/scale_factor/scale_factor/scale_factor+Omega_L+Omega_C/scale_factor/scale_factor)); // Beware only LCDM
  
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  
  gsl_real result, result2, error, epsabs, epsrel;
  //EGP size_t neval;

  auto redshift=static_cast<gsl_real>(1./scale_factor-1.);

  epsabs=1e-4;
  epsrel=1e-8;

  auto omega_m=static_cast<gsl_real>(Omega_M);
  auto omega_k=static_cast<gsl_real>(Omega_C);
  auto omega_q=static_cast<gsl_real>(Omega_L);

  gsl_function F;
  struct my_f_params params = {omega_m,omega_k,omega_q};

  F.function = &growth_var;
  F.params = &params;

  gsl_integration_qagiu(&F,redshift, epsabs, epsrel, 1000, w, &result, &error); // EGP: neval -> 1000
  gsl_integration_qagiu(&F, 0., epsabs, epsrel, 1000, w, &result2, &error); // EGP: neval -> 1000
  
  gsl_integration_workspace_free (w);
  
  
  /*
  real_prec Om=Omega_M;
  real_prec Ov=data->cosmology->omega_q;
  real_prec Ot=Om+Ov;

  real_prec x=a*(1.-Ot)/Ot;

  //Bouchet et al 1995
  //real_prec D=1.+3./x+3.*sqrt((1.+x)/x/x/x)*log(sqrt(1.+x)-sqrt(x));
  
  //Carroll,Press and Turner 1992
  real_prec D=5./2.*a*Om/((pow(Om,4./7.)-Ov+(1.+0.5*Om)*(1.+1./70.*Ov)));
  */

      
  //return D;
  //printf("D: %f\n", static_cast<real_prec>(H/H0*result/result2)); // EGP: testing
  return static_cast<real_prec>(H/H0*result/result2);
  //return H*result/cgs_km;
}
/* --- end of function D_growth --- */

//real_prec fgrow(real_prec a, real_prec Omega_M, real_prec Omega_L, real_prec H0, int term)
real_prec fgrow(real_prec a, real_prec Omega_M, real_prec Omega_L, int term)
{
   // we use LCDM and the approximation given in Lahav 1991 MNRAS,252,128

  // E = H/H0
  real_prec E = E_Hubble_a(a, Omega_M, Omega_L); // Beware only LCDM

  // Om=Om_M/(H/H0)^2/a^3
  real_prec Omega=Omega_M/((E*E)*(a*a*a));

  // f=dlnD_dlna
  real_prec f=0.0;

  switch (term)
    {	      
    case 1:
      /*
       // EGP: There is in fact a general solution in Lahav+91, but it is complicated.
       //
       if(num_1 - Omega_M - Omega_L == 0.) // for flat universes we can use an improved fit
        {
          f=pow(Omega,0.6)+1./70.*(1.-0.5*Omega*(1.+Omega));
        }
      */
      f=static_cast<real_prec>(pow(Omega,5./9.));
      break;
    case 2:
      f=static_cast<real_prec>(2.*pow(Omega,6./11.));
      break;
    case 3:
      f=static_cast<real_prec>(3.*pow(Omega,13./24.));
      break;
    }

  return(f);
}

//real_prec c_pecvel(real_prec scale_factor, real_prec Omega_M, real_prec Omega_L, real_prec H0)
real_prec c_pecvel(real_prec scale_factor, real_prec Omega_M, real_prec Omega_L, int term)
{
  // EGP: put everything in units of Mpc/h etc., so that we don't need H0 (or h) to be given.
  // Note: this requires distance units to be in Mpc/h!
  auto H0 = static_cast<real_prec>(100.); // km/s/Mpc*h; i.e. in "Hubble-units"

  real_prec f = fgrow(scale_factor, Omega_M, Omega_L, term);
  real_prec E = E_Hubble_a(scale_factor, Omega_M, Omega_L);

  //real_prec out=f*H;
  // EGP: correction, needs to be times scale factor:
  //real_prec out = f * H * scale_factor;
  real_prec out = f * H0 * E * scale_factor;

  return(out);
}

