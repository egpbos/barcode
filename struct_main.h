/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#pragma once

#include <string>         
#include <fstream>
#include <vector>

#include "define_opt.h"
#include "curses_funcs.h"
#include "fftwrapper.h"  // struct plan_pkg

using namespace std;

/* --- structures --- */
struct COSMOLOGY
{
  real_prec omega_r;			/* radiation density */
  real_prec omega_k;			/* curvature density */
  real_prec omega_m;			/* matter density */
  real_prec omega_b;			/* baryon density */
  real_prec omega_q;			/* dark energy density */
  real_prec w, wprime;		        /* de eos parameter w and first time derivative wprime=dw/da */
  real_prec sigma8;			/* mass fluctuation on 8Mpc/h scales */
  real_prec rsmooth;			/* smoothing radius: rsmooth = 8 yields sigma8 */
  real_prec n_s;			/* index of primordial dm power spectrum */
  real_prec h;			        /* Hubble parameter */
  real_prec beta;			/* lens galaxy distribution parameter 1 */
  real_prec A_spec;			/* Normalization of the power spectrum */  
  real_prec z;  
  real_prec ascale;
  real_prec mu;

  real_prec D1;
  real_prec D2;
};

struct OBSERVATIONAL
{	
  real_prec *Power;
  real_prec *signal;
  real_prec *signalX;
  real_prec *window;
  real_prec *noise_sf;
  real_prec *nobs;
  real_prec *corrf;

  real_prec biasP;
  real_prec biasE;
  real_prec Nmean_Gal;
  real_prec rho_c;
  real_prec mu;

  real_prec sigma_fac;
  real_prec sigma_min;
  real_prec delta_min;

  // EGP function pointers:
  void (*partial_f_delta_x_log_like)(struct HAMIL_DATA*, real_prec*, real_prec*);
  real_prec (*log_like)(struct HAMIL_DATA*, real_prec*);
  void (*grad_f_delta_x_comp)(struct HAMIL_DATA *, real_prec *, real_prec *, unsigned int);

  real_prec (*log_prior)(struct HAMIL_DATA*, real_prec*);
  void (*grad_log_prior)(struct HAMIL_DATA*, real_prec*, real_prec*);
  // EGP end function pointers
};

struct NUMERICAL
{	
  // EGP added:
  // some flags
  bool random_test;
  bool random_test_rsd;

  int window_type;

  int data_model;
  bool negative_obs;
  int likelihood;
  int prior;
  int sfmodel;
  bool rsd_model;
  // log file for performance statistics from HamiltonianMC
  ofstream performance_log;
  // pseudo-timestep parameters
  real_prec N_eps_fac;
  real_prec eps_fac;
  real_prec eps_fac_target;
  real_prec eps_fac_initial;
  real_prec eps_fac_power;
  real_prec s_eps_total_fac;
  real_prec s_eps_total_scaling;
  int s_eps_total_Nx_norm;
  ULONG s_eps_total;
  // initial guess:
  int initial_guess;
  string initial_guess_file;
  int initial_guess_smoothing_type;
  real_prec initial_guess_smoothing_scale;
  // mass stuff
  // bool mass_fs;
  int mass_type;
  ULONG massnum_init;
  ULONG massnum_burn;
  // output
  unsigned outnum;
  unsigned outnum_ps;
  // restarting
  unsigned start_at;
  // testing:
  real_prec mass_factor;
  int calc_h;
  // EGP: end added

  ULONG seed;

  ULONG N_bin;

  int inputmode;
  
  int INV_SUCCESS;

  int filter;

  bool PRECON;

  string codename;
  string fnamePS;
  string dataFileName;
  string fname_ps_sample;
  string fname_PROT_SPEC;
  string fname_PROT_SIGNAL_BOX;
  string fname_PROT_SIGNAL_SLICE;
  
  string dir;

  int mk;

  unsigned N1;			/* Number of cells per axis */
  unsigned N2;			/* Number of cells per axis */
  unsigned N3;			/* Number of cells per axis */
  ULONG N; // total cell number
  ULONG Nhalf;

  bool readPS;
  
  real_prec code_norm;
  
  real_prec slength;			

  ULONG iGibbs;
  ULONG N_Gibbs;
  ULONG rejections;
  ULONG total_steps_lim;
  ULONG count_attempts;  // counter for total accepted + rejected steps
                         // used for epsilon/acceptance tables

  // stuff for self-adjusting epsilon
  int eps_fac_update_type;
  unsigned N_a_eps_update;
  real_prec acc_min;
  real_prec acc_max;
  int eps_down_smooth;
  real_prec eps_up_fac;
  vector<bool> acc_flag_N_a;
  vector<real_prec> epsilon_N_a;
  // related: circular buffer for keeping track of "recent" acceptance rate
  vector<short> acc_recent;

  real_prec xllc;
  real_prec yllc;
  real_prec zllc;
  
  // Redshift space distortions
  real_prec xobs;
  real_prec yobs;
  real_prec zobs;
  bool planepar;
  bool periodic;
  
  real_prec L1;			/* Length of axis */
  real_prec L2;			/* Length of axis */
  real_prec L3;			/* Length of axis */

  real_prec vol;
  
  real_prec d1;			/* grid spacing x-direction */
  real_prec d2;			/* grid spacing y-direction */
  real_prec d3;			/* grid spacing z-direction */


  /// TESTING ///
  real_prec grad_psi_prior_factor;
  real_prec grad_psi_likeli_factor;
  bool grad_psi_prior_conjugate;
  bool grad_psi_likeli_conjugate;
  bool grad_psi_prior_times_i;
  bool grad_psi_likeli_times_i;
  bool div_dH_by_N;
  real_prec deltaQ_factor;

  int particle_kernel;
  real_prec particle_kernel_h_rel;
  real_prec particle_kernel_h;

  int N_cells;
  std::vector<int> kernel_cells_i, kernel_cells_j, kernel_cells_k;

  // FFTW
  real_prec *in_r2c;
  real_prec *out_c2r;
  complex_prec *in_c2r;
  complex_prec *out_r2c;
  struct plan_pkg *R2Cplan;
  struct plan_pkg *C2Rplan;
  // these two could be used with a common complex array to reduce copying
  // struct plan_pkg *R2Bplan; 
  // struct plan_pkg *B2Rplan;
  // these two could be used with both common complex and real arrays
  // struct plan_pkg *A2Bplan; 
  // struct plan_pkg *B2Aplan;

  bool correct_delta;

  // ~NUMERICAL() {
  // }
};



struct DATA
{
  struct COSMOLOGY *cosmology;		/* cosmological model */
  struct OBSERVATIONAL *observational;	/* numerical parameters */
  struct NUMERICAL *numerical;		/* numerical parameters */
  struct CURSES_STRUCT *curses;
  
  DATA()   // Example of a constructor used in a structure.
  {
    cosmology=new COSMOLOGY;
    observational=new OBSERVATIONAL;
    numerical=new NUMERICAL;
    // curses = new CURSES_STRUCT;
  }

  ~DATA() {
    delete cosmology;
    delete observational;
    delete numerical;
  }
};
