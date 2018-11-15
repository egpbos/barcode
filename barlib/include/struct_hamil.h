/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#pragma once

// EGP: cleared out most includes, just left the two here:
#include <gsl/gsl_rng.h>  // For the gsl_rng type // EGP: added this one
#include "define_opt.h"  // For the ULONG and real_prec types
// #include "curses_funcs.h"
// EGP: toch maar weer terug gezet, gingen dingen mis
/* #include <math.h> */
/* #include <fstream> */
/* #include <iomanip> */
// #include <iostream>
/* #include <string>          */
/* #include <string.h> */
/* #include <cassert> */
/* #include <cfloat> */
/* #include <stdlib.h> */
/* #include <stdio.h> */
/* #include <fftw3.h> */
/* #include <omp.h> */
/* #include <sstream> */
/* #include <netinet/in.h> */

/* #include <gsl/gsl_errno.h> */
/* #include <gsl/gsl_spline.h> */
/* #include <gsl/gsl_randist.h> */
/* #include <gsl/gsl_integration.h> */
/* #include <gsl/gsl_sf.h> */

/* #include "fftw_array.h" */

/* #include "define_opt.h" */

#include <vector>
#include "HMC_models.h"
#include "struct_main.h"
#include <sstream>

using namespace std;

/* --- structures --- */

struct HAMIL_NUMERICAL {
  unsigned N1;  /* Number of cells per axis */
  unsigned N2;  /* Number of cells per axis */
  unsigned N3;  /* Number of cells per axis */
  ULONG N;
  ULONG Nhalf;

  int mk;

  real_prec L1;
  real_prec L2;
  real_prec L3;

  real_prec vol;

  real_prec d1;
  real_prec d2;
  real_prec d3;

  real_prec min1;
  real_prec min2;
  real_prec min3;

  // Redshift space distortions
  real_prec xobs;
  real_prec yobs;
  real_prec zobs;
  bool planepar;
  bool periodic;

  ULONG N_bin;

  ULONG iGibbs;
  ULONG rejections;
  ULONG total_steps_lim;

  // EGP: some extra parameters for statistics on performance:
  bool accepted;
  real_prec epsilon;
  ULONG Neps;
  real_prec dH, dK, dE, dprior, dlikeli, psi_prior, psi_likeli, psi_prior_i,
  psi_prior_f, psi_likeli_i, psi_likeli_f, H_kin_i, H_kin_f;

  int method;			/* which method */
  int filter;			/* which filter method */
  ULONG itmax;			/* maximal number of iterations */
  int i_reset;			/* number to reset */
  int INV_SUCCESS;		/* indicate inversion success */

  real_prec epsi;  /* precision */

  real_prec kth;

  // EGP: pseudo-timestep parameters
  real_prec N_eps_fac;
  real_prec eps_fac;
  real_prec eps_fac_target;
  real_prec eps_fac_power;
  ULONG s_eps_total;

  // mass type
  bool mass_fs;
  bool mass_rs;
  int mass_type;
  ULONG massnum_init;
  ULONG massnum_burn;

  // EGP: testing:
  real_prec mass_factor;

  /// TESTING ///
  real_prec grad_psi_prior_factor;
  real_prec grad_psi_likeli_factor;
  bool grad_psi_prior_conjugate;
  bool grad_psi_likeli_conjugate;
  bool grad_psi_prior_times_i;
  bool grad_psi_likeli_times_i;
  bool div_dH_by_N;
  int calc_h;
  real_prec deltaQ_factor;

  bool correct_delta;

  int particle_kernel;
  real_prec particle_kernel_h;

  // FFTW
  real_prec *in_r2c;
  real_prec *out_c2r;
  complex_prec *in_c2r;
  complex_prec *out_r2c;
  struct plan_pkg *R2Cplan;
  struct plan_pkg *C2Rplan;
};

struct HAMIL_DATA {
  struct HAMIL_NUMERICAL *numerical;  // numerical parameters

  // set arrays
  real_prec *mass_f;
  real_prec *mass_r;
  // real_prec *massi;
  real_prec *gradpsi;
  real_prec *x;
  // real_prec *growth;

  real_prec *posx;
  real_prec *posy;
  real_prec *posz;

  real_prec *signal_PS;
  real_prec *window;
  real_prec *func;
  real_prec *noise;
  real_prec *nobs;
  real_prec *deltaX;
  real_prec *corrf;

  // function pointers
  void (*partial_f_delta_x_log_like)(struct HAMIL_DATA*, real_prec*,
                                     real_prec*);
  real_prec (*log_like)(struct HAMIL_DATA*, real_prec*);
  void (*grad_f_delta_x_comp)(struct HAMIL_DATA *, real_prec *, real_prec *, unsigned int);

  real_prec (*log_prior)(struct HAMIL_DATA*, real_prec*);
  void (*grad_log_prior)(struct HAMIL_DATA*, real_prec*, real_prec*);

  // model flags
  int sfmodel;
  bool rsd_model;

  // set scalars
  real_prec Nmean_Gal;
  real_prec rho_c;
  real_prec code_norm;
  real_prec mu;

  real_prec sigma_fac;
  real_prec sigma_min;
  real_prec delta_min;

  real_prec biasP;
  real_prec biasE;

  real_prec z;
  real_prec ascale;

  real_prec D1;
  real_prec D2;

  real_prec OM;
  real_prec OL;
  real_prec h;

  // external inputs
  real_prec *M;

  // curses
  // struct CURSES_STRUCT *curses;

  // stuff for SPH_kernel_3D_cells
  const int N_cells;
  const std::vector<int> kernel_cells_i, kernel_cells_j, kernel_cells_k;


  // default constructor; initialize everything yourself!
  HAMIL_DATA() : N_cells(0), kernel_cells_i(0, 0), kernel_cells_j(0, 0),
  kernel_cells_k(0, 0) {
    numerical = new HAMIL_NUMERICAL;
  }


  // secondary constructor; used by call_hamil
  HAMIL_DATA(struct DATA *data, real_prec *A, real_prec *B_f, real_prec *B_r,
             real_prec *C, real_prec *D, real_prec *E, vector<int> cells_i,
             std::vector<int> cells_j, std::vector<int> cells_k,
             int N_cells_i) :
  N_cells(N_cells_i),
  kernel_cells_i(cells_i.begin(), cells_i.end()),
  kernel_cells_j(cells_j.begin(), cells_j.end()),
  kernel_cells_k(cells_k.begin(), cells_k.end()) {
    numerical = new HAMIL_NUMERICAL;

    numerical->N1 = data->numerical->N1;
    numerical->N2 = data->numerical->N2;
    numerical->N3 = data->numerical->N3;
    numerical->N = data->numerical->N;
    numerical->Nhalf = data->numerical->Nhalf;
    numerical->L1 = data->numerical->L1;
    numerical->L2 = data->numerical->L2;
    numerical->L3 = data->numerical->L3;
    numerical->vol = data->numerical->vol;
    numerical->d1 = data->numerical->d1;
    numerical->d2 = data->numerical->d2;
    numerical->d3 = data->numerical->d3;
    numerical->min1 = data->numerical->xllc;
    numerical->min2 = data->numerical->yllc;
    numerical->min3 = data->numerical->zllc;
    // Redshift space distortions
    numerical->xobs = data->numerical->xobs;
    numerical->yobs = data->numerical->yobs;
    numerical->zobs = data->numerical->zobs;
    numerical->planepar = data->numerical->planepar;
    numerical->periodic = data->numerical->periodic;
    
    numerical->N_bin = data->numerical->N_bin;
    numerical->mk = data->numerical->mk;

    numerical->kth = data->numerical->slength;

    numerical->iGibbs = data->numerical->iGibbs;

    numerical->rejections = 0; // EGP: added for rejection rate
    numerical->total_steps_lim = data->numerical->total_steps_lim;

    numerical->itmax = 2000;      
    numerical->INV_SUCCESS = 0;/* indicate inversion success initialize with 0*/
    
    numerical->N_eps_fac = data->numerical->N_eps_fac;
    numerical->eps_fac = data->numerical->eps_fac;
    numerical->eps_fac_target = data->numerical->eps_fac_target;
    numerical->eps_fac_power = data->numerical->eps_fac_power;
    numerical->s_eps_total = data->numerical->s_eps_total;

    numerical->mass_type = data->numerical->mass_type;
    switch (numerical->mass_type) {
      case 0:  // 0: R: all ones (essentially: no mass term)
      numerical->mass_rs = true;
      numerical->mass_fs = false;
      break;
      case 1:  // 1: FS: inverse power spectrum
      numerical->mass_rs = false;
      numerical->mass_fs = true;
      break;
      case 2:  // 2: FS+FS: inverse power spectrum + likelihood force spectrum
      numerical->mass_rs = false;
      numerical->mass_fs = true;
      break;
      case 3:  // 3: FS+FS: inverse power spectrum + mean likelihood force (Wang+13)
      numerical->mass_rs = false;
      numerical->mass_fs = true;
      break;
      case 4:  // 4: FS: power spectrum
      numerical->mass_rs = false;
      numerical->mass_fs = true;
      break;
      case 5:  // 5: FS+R: inverse power spectrum + 1st order likelihood force expansion (Jasche+13)
      numerical->mass_rs = true;
      numerical->mass_fs = true;
      break;
      case 6:  // 1st order likelihood force expansion (Jasche+13) (R)
      numerical->mass_rs = true;
      numerical->mass_fs = false;
      break;
      case 60:  // type 0 until burn-in, type 6 afterwards
      numerical->mass_rs = true;
      numerical->mass_fs = false;
      break;
      default:
      stringstream message;
      message << "mass_type " << numerical->mass_type << " is not a valid value!";
      throw runtime_error(message.str());
    }

    numerical->massnum_init = data->numerical->massnum_init;
    numerical->massnum_burn = data->numerical->massnum_burn;

    gradpsi = A;
    mass_f = B_f;
    mass_r = B_r;
    posx = C;
    posy = D;
    posz = E;  

    z = data->cosmology->z;
    ascale = data->cosmology->ascale;

    D1 = data->cosmology->D1;
    D2 = data->cosmology->D2;

    OM = data->cosmology->omega_m;
    OL = data->cosmology->omega_q;
    h = data->cosmology->h;

    biasP = data->observational->biasP;
    biasE = data->observational->biasE;

    x = data->observational->signal;
    nobs = data->observational->nobs;
    deltaX = data->observational->signalX;
    noise = data->observational->noise_sf;
    corrf = data->observational->corrf;
    
    Nmean_Gal = data->observational->Nmean_Gal;
    rho_c = data->observational->rho_c;
    code_norm = data->numerical->code_norm;
    mu = data->observational->mu;

    sigma_fac = data->observational->sigma_fac; 
    sigma_min = data->observational->sigma_min; 
    delta_min = data->observational->delta_min; 

    signal_PS = data->observational->Power;
    
    window = data->observational->window;

    // function pointers
    log_like = data->observational->log_like;
    partial_f_delta_x_log_like = data->observational->partial_f_delta_x_log_like;
    grad_f_delta_x_comp = data->observational->grad_f_delta_x_comp;
    log_prior = data->observational->log_prior;
    grad_log_prior = data->observational->grad_log_prior;

    // model flags
    sfmodel = data->numerical->sfmodel;
    rsd_model = data->numerical->rsd_model;

    // Testing:
    numerical->mass_factor = data->numerical->mass_factor;
    numerical->grad_psi_prior_factor = data->numerical->grad_psi_prior_factor;
    numerical->grad_psi_likeli_factor = data->numerical->grad_psi_likeli_factor;
    numerical->grad_psi_prior_conjugate = data->numerical->grad_psi_prior_conjugate;
    numerical->grad_psi_likeli_conjugate = data->numerical->grad_psi_likeli_conjugate;
    numerical->grad_psi_prior_times_i = data->numerical->grad_psi_prior_times_i;
    numerical->grad_psi_likeli_times_i = data->numerical->grad_psi_likeli_times_i;
    numerical->div_dH_by_N = data->numerical->div_dH_by_N;
    numerical->calc_h = data->numerical->calc_h;
    numerical->deltaQ_factor = data->numerical->deltaQ_factor;

    numerical->particle_kernel = data->numerical->particle_kernel;
    numerical->particle_kernel_h = data->numerical->particle_kernel_h;

    // FFTW
    numerical->in_r2c = data->numerical->in_r2c;
    numerical->out_c2r = data->numerical->out_c2r;
    numerical->in_c2r = data->numerical->in_c2r;
    numerical->out_r2c = data->numerical->out_r2c;
    numerical->R2Cplan = data->numerical->R2Cplan;
    numerical->C2Rplan = data->numerical->C2Rplan;

    // curses
    // curses = data->curses;

    numerical->correct_delta = data->numerical->correct_delta;
  }    

  ~HAMIL_DATA() {
    delete numerical;
  }
};

