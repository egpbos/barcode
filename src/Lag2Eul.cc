/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <cmath> // sqrt

#include "define_opt.h"
#include "struct_main.h"
#include "fftw_array.h"

#include "cosmo.h"
#include "EqSolvers.h"
#include "massFunctions.h"
#include "convolution.hpp" // convcomp, convcompb

#include "disp_part.h"

#include "convenience.h"

#include "rsd.h"

using namespace std;

// A few ideas for ordering input parameters of functions in such a way that
// they aren't totally hidden away (like now with the data/hamil_data structs),
// but aren't overly explicit (and long) either.

struct GridSize {
  ULONG N;
  ULONG Nhalf;
  int x;
  int y;
  int z;
};

// for total simulation volume, cell sizes, origin
struct RealVec3 {
  real_prec x;
  real_prec y;
  real_prec z;
};

struct Cosmology {
  real_prec Omega_M;
  real_prec Omega_L;
  real_prec hconst;
  real_prec D1;
  real_prec D2;
  real_prec scale;
};

// Really though, the in and outs should be DensityField objects that have
// gridsize, boxsize, cellsize and origin builtin. That would really clean up
// the mess of parameters in functions. Oh well.
void Lag2Eul_zeldovich(const real_prec *in, real_prec *out, real_prec *posx,
                       real_prec *posy, real_prec *posz,
                       const GridSize &gridsize, const RealVec3 &boxsize,
                       const RealVec3 &cellsize, const RealVec3 &origin,
                       const Cosmology &cosmo, int masskernel,
                       real_prec kth, int facL, bool reggrid,
                       gsl_rng * gBaseRand, real_prec kernel_scale_factor,
                       plan_pkg *R2Cplan, plan_pkg *C2Rplan);


void Lag2Eul_zeldovich(real_prec *in, real_prec *out, real_prec *posx, real_prec *posy, real_prec *posz,
                       unsigned int N1, unsigned int N2,
                       unsigned int N3, real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2,
                       real_prec d3,
                       real_prec min1, real_prec min2, real_prec min3, real_prec D1, real_prec scale, real_prec Omega_M,
                       real_prec Omega_L, int masskernel, unsigned int facL, bool reggrid, gsl_rng *gBaseRand,
                       real_prec kernel_scale_factor, plan_pkg *R2Cplan, plan_pkg *C2Rplan)
{
  ULONG N=N1*N2*N3;

  ULONG NL1=N1*facL, NL2=N2*facL, NL3=N3*facL;
  ULONG NLL=NL1*NL2*NL3;

  fftw_array<real_prec>  psix(N), psiy(N), psiz(N);
  fftw_array<real_prec>  vex(N), vey(N), vez(N);

  // compute components of Psi
  // delta^(1) == input field (in)
  // phi = -D1*delta^(1) -> out (used as dummy array)
  multiply_factor_array(-D1, in, out, N);

  // Psi components
  bool zeropad = false;
  bool norm = false;  // true -> velocities, false -> displacements
  theta2vel(N1, N2, N3, L1, L2, L3, scale, Omega_M, Omega_L, out, psix, psiy, psiz, zeropad, norm, R2Cplan, C2Rplan);
  // EGP: cellbound interpolates to cell boundaries (corners), but we want
  // the values in the cell centers!
  //cellboundcomp(N1,N2,N3,L1,L2,L3,psix);
  //cellboundcomp(N1,N2,N3,L1,L2,L3,psiy);
  //cellboundcomp(N1,N2,N3,L1,L2,L3,psiz);

  // compute the resulting particle positions & velocities and density.
  {
    fftw_array<real_prec> dummyL(NLL);
    fill_one(dummyL, N);
    // FIXME
    // Eigenlijk moet deze alleen 1 zijn als er een deeltje in zit,
    // wat in disp_part wordt bepaald, maar alleen als de boel niet
    // periodic is, dus for now (nu alles periodic is) kunnen we dit
    // wel ff zo laten.
    // Het stelt namelijk de particle mass voor en die is nul als het
    // deeltje er niet meer bij zit.

    disp_part(N1,N2,N3,L1,L2,L3,d1,d2,d3,posx,posy,posz,psix,psiy,psiz,dummyL,facL,reggrid,gBaseRand);

    switch (masskernel)
    {
      case 0:
        getDensity_NGP(N1,N2,N3,L1,L2,L3,d1,d2,d3,min1,min2,min3,posx,posy,posz,dummyL,NLL,out);
        break;
      case 1:
        getDensity_CIC(N1,N2,N3,L1,L2,L3,d1,d2,d3,min1,min2,min3,posx,posy,posz,dummyL,NLL,out, true);
        break;
      case 2:
        getDensity_TSC(N1,N2,N3,L1,L2,L3,d1,d2,d3,min1,min2,min3,posx,posy,posz,dummyL,NLL,out);
        break;
      case 3:
        getDensity_SPH(N1, N2, N3, L1, L2, L3, d1, d2, d3, min1, min2, min3, posx, posy, posz, dummyL, NLL, out, true, kernel_scale_factor);
        break;
    }
  }

  overdens(N1,N2,N3,out,out);
}


// --------------------------------------------


void Lag2Eul_non_zeldovich(real_prec *in, real_prec *dummy, real_prec *posx, real_prec *posy, real_prec *posz,
                           unsigned int N1,
                           unsigned int N2, unsigned int N3, real_prec L1, real_prec L2, real_prec L3, real_prec d1,
                           real_prec d2,
                           real_prec d3, real_prec min1, real_prec min2, real_prec min3, real_prec D1, real_prec D2,
                           real_prec scale, real_prec Omega_M, real_prec Omega_L, int masskernel, real_prec kth,
                           unsigned int facL, bool reggrid, gsl_rng *gBaseRand, const std::string& dir,
                           real_prec kernel_scale_factor)
{
  real_prec kthsc=kth;
#ifdef TRANSF
  int ftype = 1;
#endif // TRANSF
#ifdef TRANSFSC
  int ftype = 1;
#endif // TRANSFSC

  ULONG N=N1*N2*N3;

  ULONG NL1=N1*facL, NL2=N2*facL, NL3=N3*facL;
  ULONG NLL=NL1*NL2*NL3;


  fftw_array<real_prec>  psix(N), psiy(N), psiz(N);
  fftw_array<real_prec>  vex(N), vey(N), vez(N);

  {
    fftw_array<real_prec> dummy2(N), dummy3(N);

    copyArray(in, dummy, N);

    // delta(1)-> Phi(1)
    PoissonSolver(N1,N2,N3,L1,L2,L3,dummy,dummy2);
    // Phi(1)  -> delta(2)
    calc_m2v_mem(N1, N2, N3, L1, dummy2, dummy3);

#ifdef TRANSF
    {
      {
        string fname= dir + string("auxtransfzeld");
        get_scalar(fname,dummy2,N1,N2,N3);
      }
      convcompb(L1,L2,L3,d1,d2,d3,N1,N2,N3,dummy,dummy2,ftype);
      copyArray(dummy2, dummy, N);
    }
    {
      {
        string fname= dir + string("auxtransf2lpt");
        get_scalar(fname,dummy2,N1,N2,N3);
      }
      convcompb(L1,L2,L3,d1,d2,d3,N1,N2,N3,dummy3,dummy2,ftype);
      copyArray(dummy2, dummy3, N);
    }
#endif /* TRANSF */

    //div Psi^2LPT
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
    for(ULONG i=0;i<N;i++)
      dummy2[i]=D1*dummy[i]-D2*dummy3[i];

    //K o div Psi^2LPT
    convcomp(N1, N2, N3, dummy2, dummy2, kth, dir);

    /* Psi^tot=K o Psi^2LPT + Psi^SC - K o Psi^SC */

    fftw_array<real_prec> dummy4(N);
    //div Psi^SC
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
    for(ULONG i=0;i<N;i++)
      {
        real_prec psilin=-D1*dummy[i];
        real_prec psisc=0.;

        if (1.+2./3.*psilin > 0.)
          psisc=static_cast<real_prec>(3.*(std::sqrt(1.+2./3.*psilin)-1.));
        else
          psisc=-3.;

        psisc*=static_cast<real_prec>(-1.);

        dummy4[i]=psisc;
      }


#ifdef TRANSFSC
    fftw_array<real_prec> dummy2(N);
    {
      string fname= dir + string("auxtransfzeld");
      get_scalar(fname,dummy,N1,N2,N3);
    }
    convcompb(L1,L2,L3,d1,d2,d3,N1,N2,N3,dummy4,dummy,ftype);
    copyArray(dummy, dummy4, N);
#endif

    /* Psi^tot x-component */
    //K o Psi^2LPT_x
    theta2velcomp(N1, N2, N3, L1, L2, L3, scale, Omega_M, Omega_L, dummy2, dummy3, false, false, 1);
   //Psi^SC_x
    theta2velcomp(N1, N2, N3, L1, L2, L3, scale, Omega_M, Omega_L, dummy4, dummy, false, false, 1);
    //K o Psi^2LPT_x + Psi^SC_x
    add_to_array(dummy, dummy3, N);
    //K o Psi^SC_x
    convcomp(N1, N2, N3, dummy, dummy, kthsc, dir);
    //K o Psi^2LPT_x + Psi^SC_x - K o Psi^SC_x
    subtract_arrays(dummy3, dummy, psix, N);
    cellboundcomp(N1, N2, N3, psix);

    /* Psi^tot y-component */

    //K o Psi^2LPT_y
    theta2velcomp(N1, N2, N3, L1, L2, L3, scale, Omega_M, Omega_L, dummy2, dummy3, false, false, 2);
    //Psi^SC_y
    theta2velcomp(N1, N2, N3, L1, L2, L3, scale, Omega_M, Omega_L, dummy4, dummy, false, false, 2);
    //K o Psi^2LPT_y + Psi^SC_y
    add_to_array(dummy, dummy3, N);
    //K o Psi^SC_y
    convcomp(N1, N2, N3, dummy, dummy, kthsc, dir);
    //K o Psi^2LPT_y + Psi^SC_y - K o Psi^SC_y
    subtract_arrays(dummy3, dummy, psiy, N);
    cellboundcomp(N1, N2, N3, psiy);

    /* Psi^tot z-component */
    //K o Psi^2LPT_z
    theta2velcomp(N1, N2, N3, L1, L2, L3, scale, Omega_M, Omega_L, dummy2, dummy3, false, false, 3);
    //Psi^SC_z
    theta2velcomp(N1, N2, N3, L1, L2, L3, scale, Omega_M, Omega_L, dummy4, dummy, false, false, 3);
    //K o Psi^2LPT_z + Psi^SC_z
    add_to_array(dummy, dummy3, N);
    //K o Psi^SC_z
    convcomp(N1, N2, N3, dummy, dummy, kthsc, dir);
    //K o Psi^2LPT_z + Psi^SC_z - K o Psi^SC_z
    subtract_arrays(dummy3, dummy, psiz, N);
    cellboundcomp(N1, N2, N3, psiz);
  }

  /* end Psi^tot */

  fftw_array<real_prec> dummyL(NLL);
  fill_one(dummyL, N);

  disp_part(N1,N2,N3,L1,L2,L3,d1,d2,d3,posx,posy,posz,psix,psiy,psiz,dummyL,facL,reggrid,gBaseRand);

  {
  switch (masskernel)
    {
    case 0:
      {
        getDensity_NGP(N1,N2,N3,L1,L2,L3,d1,d2,d3,min1,min2,min3,posx,posy,posz,dummyL,NLL,dummy);
      }
      break;
    case 1:
      {
        getDensity_CIC(N1,N2,N3,L1,L2,L3,d1,d2,d3,min1,min2,min3,posx,posy,posz,dummyL,NLL,dummy, true);
      }
      break;
    case 2:
      {
        getDensity_TSC(N1,N2,N3,L1,L2,L3,d1,d2,d3,min1,min2,min3,posx,posy,posz,dummyL,NLL,dummy);
      }
      break;
    case 3:
      getDensity_SPH(N1, N2, N3, L1, L2, L3, d1, d2, d3, min1, min2, min3, posx, posy, posz, dummyL, NLL, dummy, true, kernel_scale_factor);
      break;
    }
  }
  {
    fftw_array<real_prec> dummy2(N);
    overdens(N1,N2,N3,dummy,dummy2);
    copyArray(dummy2, dummy, N);
  }
}


// --------------------------------------------


void Lag2Eul(real_prec *in, real_prec *out, real_prec *posx, real_prec *posy, real_prec *posz, unsigned int N1,
             unsigned int N2, unsigned int N3, real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2,
             real_prec d3, real_prec min1, real_prec min2, real_prec min3, real_prec D1, real_prec D2, real_prec scale,
             real_prec Omega_M, real_prec Omega_L, int sfmodel, int masskernel, real_prec kth, unsigned int facL,
             bool reggrid, gsl_rng *gBaseRand, const std::string& dir, real_prec kernel_scale_factor, plan_pkg *R2Cplan,
             plan_pkg *C2Rplan)
{
  if (sfmodel == 1) {
    Lag2Eul_zeldovich(in, out, posx, posy, posz, N1, N2, N3, L1, L2, L3, d1, d2, d3, min1, min2, min3, D1, scale, Omega_M,
                      Omega_L, masskernel, facL,
                      reggrid, gBaseRand, kernel_scale_factor, R2Cplan, C2Rplan);
  } else {
    Lag2Eul_non_zeldovich(in, out, posx, posy, posz, N1, N2, N3, L1, L2, L3, d1, d2, d3, min1,min2, min3, D1, D2, scale, Omega_M, Omega_L, masskernel, kth, facL, reggrid, gBaseRand, dir, kernel_scale_factor);
  }
}


// --------------------------------------------


void Lag2Eul_rsd_zeldovich(real_prec *in, real_prec *out, real_prec *posx, real_prec *posy, real_prec *posz,
                           unsigned int N1, unsigned int N2,
                           unsigned int N3, real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2,
                           real_prec d3,
                           real_prec min1, real_prec min2, real_prec min3, real_prec D1, real_prec scale,
                           real_prec Omega_M,
                           real_prec Omega_L, int masskernel, unsigned int facL, bool reggrid, gsl_rng *gBaseRand,
                           real_prec kernel_scale_factor, real_prec xobs, real_prec yobs, real_prec zobs, bool planepar,
                           bool periodic, plan_pkg *R2Cplan, plan_pkg *C2Rplan)
{
  ULONG N=N1*N2*N3;

  ULONG NL1=N1*facL, NL2=N2*facL, NL3=N3*facL;
  ULONG NLL=NL1*NL2*NL3;

  fftw_array<real_prec>  psix(N), psiy(N), psiz(N);
  fftw_array<real_prec>  vex(N), vey(N), vez(N);

  // compute components of Psi
  // delta^(1) == input field (in)
  // phi = -D1*delta^(1) -> out (used as dummy array)
  // multiply_factor_array(-D1, in, out, N);
  // R2Cplan->R as dummy helps speed up theta2vel!
  multiply_factor_array(-D1, in, R2Cplan->R, N);

  // Psi components
  bool zeropad = false;
  bool norm = false;  // true -> velocities, false -> displacements
  theta2vel(N1, N2, N3, L1, L2, L3, scale, Omega_M, Omega_L, R2Cplan->R, psix, psiy, psiz, zeropad, norm, R2Cplan,
            C2Rplan);

  // EGP: cellbound interpolates to cell boundaries (corners), but we want
  // the values in the cell centers!
  //cellboundcomp(N1,N2,N3,L1,L2,L3,psix);
  //cellboundcomp(N1,N2,N3,L1,L2,L3,psiy);
  //cellboundcomp(N1,N2,N3,L1,L2,L3,psiz);

  // EGP: calculate velocities (disp_part doesn't do that!)
  // norm = true; // this makes sure it multiplies the psi's with the velocity factor
  // theta2vel(N1,N2,N3,L1,L2,L3,scale,Omega_M,Omega_L,hconst, out, vex, vey, vez, redshift, zeropad, norm);
  real_prec cpecvel = c_pecvel(scale, Omega_M, Omega_L, 1);
  multiply_factor_array(cpecvel, psix, vex, N);
  multiply_factor_array(cpecvel, psiy, vey, N);
  multiply_factor_array(cpecvel, psiz, vez, N);

  // compute the resulting particle positions and density.
  {
    fftw_array<real_prec> dummyL(NLL);
    fill_one(dummyL, N);
    // FIXME
    // Eigenlijk moet deze alleen 1 zijn als er een deeltje in zit,
    // wat in disp_part wordt bepaald, maar alleen als de boel niet
    // periodic is, dus for now (nu alles periodic is) kunnen we dit
    // wel ff zo laten.
    // Het stelt namelijk de particle mass voor en die is nul als het
    // deeltje er niet meer bij zit.

    disp_part(N1,N2,N3,L1,L2,L3,d1,d2,d3,posx,posy,posz,psix,psiy,psiz,dummyL,facL,reggrid,gBaseRand);
    //quick_dump_scalar(posz, N1, "posz_pre_rsd", 0, false);
    //quick_dump_scalar(vez, N1, "vez_pre_rsd", 0, false);

    // Here we do the RSD transformation
    calc_pos_rsd(NLL, L3, xobs, yobs, zobs, posx, posy, posz, vex, vey, vez, posx, posy, posz, scale, Omega_M, Omega_L,
                 planepar, periodic);

    //quick_dump_scalar(posz, N1, "posz_post_rsd", 0, false);
    //quick_dump_scalar(vez, N1, "vez_post_rsd", 0, false);

    switch (masskernel)
    {
      case 0:
        getDensity_NGP(N1,N2,N3,L1,L2,L3,d1,d2,d3,min1,min2,min3,posx,posy,posz,dummyL,NLL,out);
        break;
      case 1:
        getDensity_CIC(N1,N2,N3,L1,L2,L3,d1,d2,d3,min1,min2,min3,posx,posy,posz,dummyL,NLL,out, true);
        break;
      case 2:
        getDensity_TSC(N1,N2,N3,L1,L2,L3,d1,d2,d3,min1,min2,min3,posx,posy,posz,dummyL,NLL,out);
        break;
      case 3:
        getDensity_SPH(N1, N2, N3, L1, L2, L3, d1, d2, d3, min1, min2, min3, posx, posy, posz, dummyL, NLL, out, true, kernel_scale_factor);
        break;
    }
  }

  overdens(N1,N2,N3,out,out);
}


