/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <fftw3.h>
#include <sstream>
#include <climits>

#include <gsl/gsl_integration.h>

#include "define_opt.h"

#include "fftwrapper.h"
#include "ini_reader.hpp"

#include "struct_main.h"
#include "cosmo.h"

#include "HMC_models.h"
#include "hmc/prior/gaussian.hpp"
#include "hmc/likelihood/gaussian_random_field.hpp"
#include "hmc/likelihood/poissonian.hpp"
#include "hmc/likelihood/gaussian_independent.hpp"
#include "hmc/likelihood/lognormal_independent.hpp"

using namespace std;

int cmbcosm=3;// 1: WMAP3, 2: WMAP7 mean, 3: WMAP7 maximum. take care!
//power spectrum is calculated only approximately with fitting formulae!

void INIT_PARAMS(struct DATA *data) {
  data->numerical->codename=string("BARCODE");

  string codename=data->numerical->codename;
  cout << " >>> "<<codename<<": initialising settings...."<<endl<<endl;

  data->numerical->rejections = 0; // EGP: start counting

  char fname_PAR[100]="input.par";
  parameter_inifile params (fname_PAR);

  auto inputmode=params.find<int>("inputmode");

  string dataFileName=params.find<string>("file");
  cout << "file: " << dataFileName << endl;

  // EGP: this was lost, causing fuckups:
  data->numerical->INV_SUCCESS = 0;

  // EGP: added following switches
  auto random_test = params.find<bool>("random_test");
  auto random_test_rsd = params.find<bool>("random_test_rsd");

  auto window_type = params.find<int>("window_type");

  auto data_model = params.find<int>("data_model");
  int negative_obs = params.find<bool>("negative_obs");

  auto likelihood = params.find<int>("likelihood");
  auto prior = params.find<int>("prior");

  auto sfmodel = params.find<int>("sfmodel");
  auto rsd_model = params.find<bool>("rsd_model");

  if (((data_model == 1) && (likelihood != 2)) ||
      ((likelihood == 2) && (data_model != 1)) ) {
    throw runtime_error("Error: incompatible data_model and likelihood in input.par! Logarithmic and log-normal must go together, or you must choose other models.");
  }
  // EGP: end added switches

  // EGP: pseudo-timestep parameters
  auto N_eps_fac = params.find<real_prec>("N_eps_fac");
  auto eps_fac = params.find<real_prec>("eps_fac");
  auto eps_fac_initial = params.find<real_prec>("eps_fac_initial");
  auto eps_fac_power = params.find<real_prec>("eps_fac_power");
  auto s_eps_total_fac = params.find<real_prec>("s_eps_total_fac");
  auto s_eps_total_scaling = params.find<real_prec>("s_eps_total_scaling");
  auto s_eps_total_Nx_norm = params.find<int>("s_eps_total_Nx_norm");

  // EGP: adjusting epsilon parameters
  auto eps_fac_update_type = params.find<int>("eps_fac_update_type");
  auto N_a_eps_update = params.find<unsigned>("N_a_eps_update");
  auto acc_min = params.find<real_prec>("acc_min");
  auto acc_max = params.find<real_prec>("acc_max");
  auto eps_down_smooth = params.find<int>("eps_down_smooth");
  auto eps_up_fac = params.find<real_prec>("eps_up_fac");

  // EGP: initial guess
  auto initial_guess = params.find<int>("initial_guess");
  string initial_guess_file = params.find<string>("initial_guess_file");
  auto initial_guess_smoothing_type = params.find<int>("initial_guess_smoothing_type");
  auto initial_guess_smoothing_scale = params.find<real_prec>("initial_guess_smoothing_scale");

  // mass type
  // bool mass_fs = params.find<bool>("mass_fs");
  auto mass_type = params.find<int>("mass_type");
  auto massnum_init = params.find<ULONG>("massnum_burn");
  auto massnum_burn = params.find<ULONG>("massnum_post");

  // output
  auto outnum = params.find<unsigned>("outnum");
  auto outnum_ps = params.find<unsigned>("outnum_ps");

  string dir=params.find<string>("dir");

  auto N1=params.find<unsigned>("Nx");
  unsigned N2=N1;
  unsigned N3=N1;

  auto L1=params.find<real_prec>("Lx");
  real_prec L2=L1;
  real_prec L3=L1;

  auto xllc=params.find<real_prec>("xllc");
  auto yllc=params.find<real_prec>("yllc");
  auto zllc=params.find<real_prec>("zllc");

  auto xobs=params.find<real_prec>("xobs");
  auto yobs=params.find<real_prec>("yobs");
  auto zobs=params.find<real_prec>("zobs");
  auto planepar = params.find<bool>("planepar");
  auto periodic = params.find<bool>("periodic");

  auto seed=params.find<ULONG>("seed");

  auto N_bin=params.find<ULONG>("N_bin");
  auto N_Gibbs=params.find<ULONG>("N_Gibbs"); // EGP: added
  auto total_steps_lim = params.find<ULONG>("total_steps_lim");

  auto masskernel=params.find<int>("masskernel");

  auto z=params.find<real_prec>("z");
  auto ascale = static_cast<real_prec>(1./(1.+z));

  auto slength=params.find<real_prec>("slength");

  auto readPS=params.find<bool>("readPS");

  auto sigma_fac = params.find<real_prec>("sigma_fac");
  auto sigma_min = params.find<real_prec>("sigma_min");
  auto delta_min = params.find<real_prec>("delta_min");

  auto correct_delta = params.find<bool>("correct_delta");

  string fnamePS=params.find<string>("fnamePS");
  cout << "fnamePS: " << fnamePS << endl;

  data->cosmology->z = z;
  cout<<"---> attention: redshift= "<<data->cosmology->z<<endl;
  data->cosmology->ascale = ascale;

  data->numerical->seed = seed;
  cout<<"---> attention: seed= "<<data->numerical->seed<<endl;

  data->numerical->slength = slength;
  cout<<"---> attention: slength= "<<data->numerical->slength<<endl;

  data->numerical->inputmode = inputmode;

  data->numerical->random_test = random_test;
  data->numerical->random_test_rsd = random_test_rsd;

  data->numerical->window_type = window_type;

  // EGP: absolute! Not w.r.t. dir.
  data->numerical->dataFileName = dataFileName;

  data->numerical->dir = dir;

  data->numerical->N_bin = N_bin;
  data->numerical->N_Gibbs = N_Gibbs; // EGP: added
  if (total_steps_lim > 0) {
    data->numerical->total_steps_lim = total_steps_lim;
  } else {
    data->numerical->total_steps_lim = ULONG_MAX;
  }
  data->numerical->count_attempts = 0;  // initialize at zero
  // N.B.: when restarting, this will also stay at zero! In barcoderunner.cc,
  // the performance log is read and this could also be used to reset this
  // counter to its pre-restart value. However, after a second restart, this
  // will no longer be accurate, because the part of the log that was written
  // after the attempt in the first run from which the restart will continue
  // is not properly taken into account.

  data->numerical->readPS = readPS;

  // EGP: also absolute!
  data->numerical->fnamePS = fnamePS;

  data->numerical->mk = masskernel;	

  data->numerical->data_model = data_model;
  data->numerical->negative_obs = negative_obs;
  data->numerical->sfmodel = sfmodel;
  data->numerical->rsd_model = rsd_model;
  data->numerical->likelihood = likelihood;
  data->numerical->prior = prior;

  data->numerical->initial_guess = initial_guess;
  data->numerical->initial_guess_file = initial_guess_file;
  data->numerical->initial_guess_smoothing_type = initial_guess_smoothing_type;
  data->numerical->initial_guess_smoothing_scale = initial_guess_smoothing_scale;

  data->numerical->N1 = N1;		/* Number of cells per axis */
  cout<<"---> attention: N1= "<<data->numerical->N1<<endl;
  data->numerical->N2 = N2;		/* Number of cells per axis */
  data->numerical->N3 = N3;		/* Number of cells per axis */

  ULONG N = N1 * N2 * N3;
  data->numerical->N = N;
  data->numerical->Nhalf = N1 * N2 * (N3/2 + 1);

  data->numerical->L1 = L1;		/* Length of axis */
  data->numerical->L2 = L2;		/* Length of axis */
  data->numerical->L3 = L3;		/* Length of axis */

  data->numerical->vol = L1 * L2 * L3;

  data->numerical->xllc = xllc;
  data->numerical->yllc = yllc;
  data->numerical->zllc = zllc;

  data->numerical->xobs = xobs;
  data->numerical->yobs = yobs;
  data->numerical->zobs = zobs;
  data->numerical->planepar = planepar;
  data->numerical->periodic = periodic;

  data->numerical->d1=L1/real_prec(N1);		/* grid spacing x-direction */
  data->numerical->d2=L2/real_prec(N2);		/* grid spacing y-direction */
  data->numerical->d3=L3/real_prec(N3);		/* grid spacing z-direction */

  // "Observational" parameters
  data->observational->sigma_fac = sigma_fac;
  data->observational->sigma_min = sigma_min;
  data->observational->delta_min = delta_min;

  // EGP: pseudo-timestep parameters
  // step number
  data->numerical->N_eps_fac = N_eps_fac;
  // step size
  data->numerical->eps_fac_update_type = eps_fac_update_type;
  if (eps_fac > 0)
  {
    data->numerical->eps_fac_target = eps_fac;
  }
  else
  {
    // power law fit (heuristic (calibration/epsilon/fit_powerlaw.py))
    data->numerical->eps_fac_target = 2.38902581 * pow(N, -0.57495347);
    cout << "  eps_fac_target determined to be " << data->numerical->eps_fac_target << endl;
  }
  if (eps_fac_initial > 0)
  {
    data->numerical->eps_fac_initial = eps_fac_initial;
  }
  else
  {
    data->numerical->eps_fac_initial = data->numerical->eps_fac_target;
  }
  // set actual initial eps_fac for running code:
  switch (eps_fac_update_type) {
    case 0:
      data->numerical->eps_fac = data->numerical->eps_fac_target;
      break;
    case 1:
      data->numerical->eps_fac = data->numerical->eps_fac_initial;
      break;
    case 2:
    case 3:
      if (eps_fac > 0) {
        data->numerical->eps_fac = data->numerical->eps_fac_target;
      } else {
        data->numerical->eps_fac = 2.;
      }
      break;
  }
  
  // other eps parameters (updating scheme):
  data->numerical->eps_fac_power = eps_fac_power;
  data->numerical->s_eps_total_fac = s_eps_total_fac;
  data->numerical->s_eps_total_scaling = s_eps_total_scaling;
  data->numerical->s_eps_total_Nx_norm = s_eps_total_Nx_norm;
  real_prec grondtal = static_cast<real_prec>(N)/static_cast<real_prec>(pow(s_eps_total_Nx_norm,3));
  data->numerical->s_eps_total = static_cast<ULONG>(ceil( s_eps_total_fac * pow(grondtal, s_eps_total_scaling) ));
  if (data->numerical->s_eps_total <= 0)
    data->numerical->s_eps_total = 1;
  cout << "  s_eps_total determined to be " << data->numerical->s_eps_total << endl;

  // eps_fac_update_type 2
  data->numerical->N_a_eps_update = N_a_eps_update;
  data->numerical->acc_min = acc_min;
  data->numerical->acc_max = acc_max;
  data->numerical->eps_down_smooth = eps_down_smooth;
  data->numerical->eps_up_fac = eps_up_fac;
  data->numerical->acc_flag_N_a = vector<bool>(N_a_eps_update);
  data->numerical->epsilon_N_a = vector<real_prec>(N_a_eps_update);
  // fill initially with initial eps_fac, so that there are no zeroes on restart
  // (which would break the update function, which searches a minimum value,
  // which should never be zero!)
  fill(data->numerical->epsilon_N_a.begin(), data->numerical->epsilon_N_a.end(),
       data->numerical->eps_fac);

  data->numerical->acc_recent = vector<short>(N_a_eps_update);


  // Mass type parameters
  // data->numerical->mass_fs = mass_fs;
  data->numerical->mass_type = mass_type;

  // How often recompute the mass (before and after "burn-in")
  if (massnum_init > 0) {
    data->numerical->massnum_init = massnum_init;
  } else {
    data->numerical->massnum_init = N_Gibbs;
  }
  if (massnum_burn > 0) {
    data->numerical->massnum_burn = massnum_burn;
  } else {
    data->numerical->massnum_burn = N_Gibbs;
  }

  // output
  data->numerical->outnum = outnum;
  data->numerical->outnum_ps = outnum_ps;

  // Testing:
  auto mass_factor = params.find<real_prec>("mass_factor");
  data->numerical->mass_factor = mass_factor;
  data->numerical->grad_psi_prior_factor = params.find<real_prec>("grad_psi_prior_factor");
  data->numerical->grad_psi_likeli_factor = params.find<real_prec>("grad_psi_likeli_factor");
  data->numerical->grad_psi_prior_conjugate = params.find<bool>("grad_psi_prior_conjugate");
  data->numerical->grad_psi_likeli_conjugate = params.find<bool>("grad_psi_likeli_conjugate");
  data->numerical->grad_psi_prior_times_i = params.find<bool>("grad_psi_prior_times_i");
  data->numerical->grad_psi_likeli_times_i = params.find<bool>("grad_psi_likeli_times_i");
  data->numerical->div_dH_by_N = params.find<bool>("div_dH_by_N");
  data->numerical->calc_h = params.find<int>("calc_h");

  data->numerical->deltaQ_factor = params.find<real_prec>("deltaQ_factor");


  data->numerical->particle_kernel = params.find<int>("particle_kernel");

  data->numerical->particle_kernel_h_rel = params.find<real_prec>("particle_kernel_h_rel");
  if (data->numerical->particle_kernel_h_rel != 1.)
  {
    cout << "!!! Warning: particle_kernel_h_rel other than 1 for the SPH cubic spline kernel will cause problems! If it is smaller, not enough contributions from surrounding cells will be included. If it is larger, it will smooth the field too much, wiping out small scale structure." << endl;

#ifdef DEBUG
    char type;
    do
    {
      cout << "Do you want to continue with particle_kernel_h_rel = " << data->numerical->particle_kernel_h_rel << "? [y/n]" << endl;
      cin >> type;
    }
    while( !cin.fail() && type!='y' && type!='n' );

    if (type == 'n') {
      throw runtime_error("entered n, aborting.");
    }
#endif  // DEBUG
  }
  if (data->numerical->particle_kernel_h_rel > static_cast<real_prec>(data->numerical->N1)/4) {
    throw runtime_error("A particle_kernel_h_rel of more than Nx/4 cells does not make sense for SPH kernels, since the diameter of the kernel will be larger than the width of the box!\nAborting.");
  }

  real_prec d1 = data->numerical->d1, d2 = data->numerical->d2, d3 = data->numerical->d3;
  real_prec avg_cell_size = (d1 + d2 + d3)/3.;
  data->numerical->particle_kernel_h = data->numerical->particle_kernel_h_rel * avg_cell_size;

  // prepare SPH kernel cells vectors
  int particle_kernel_type = data->numerical->particle_kernel;
  real_prec particle_kernel_h = data->numerical->particle_kernel_h;
  int N_cells = SPH_kernel_3D_cells_count(particle_kernel_type,
                                          particle_kernel_h, d1, d2, d3);
  std::vector<int> cells_i, cells_j, cells_k;
  SPH_kernel_3D_cells(particle_kernel_type, particle_kernel_h, d1, d2, d3,
                      cells_i, cells_j, cells_k);

  // std::cerr << "cells_ijk:" << std::endl;
  // for (int ix = 0; ix < N_cells; ++ix)
  //   std::cerr << cells_i[ix] << " " << cells_j[ix] << " " << cells_k[ix] << std::endl;
  // std::cerr << endl;

  // std::vector< std::pair<int, int> > ij_out;
  // std::vector<int> k_begin, k_last;
  // SPH_kernel_3D_cells_hull_1(cells_i, cells_j, cells_k, ij_out, k_begin, k_last);

  // std::cerr << "cells_ij{k_b}{k_l}:" << std::endl;
  // int new_N = ij_out.size();
  // std::cout << new_N << " " << k_begin.size() << " " << k_last.size() << std::endl;
  // for (int ix = 0; ix < new_N; ++ix)
  //   std::cerr << ij_out[ix].first << " " << ij_out[ix].second << " "
  //            << k_begin[ix] << " " << k_last[ix] << std::endl;

  // exit(1);

  data->numerical->N_cells = N_cells;
  data->numerical->kernel_cells_i = cells_i;
  data->numerical->kernel_cells_j = cells_j;
  data->numerical->kernel_cells_k = cells_k;

  data->numerical->correct_delta = correct_delta;

  cout<<endl;
}

void INIT_FFTW(struct DATA *data, real_prec *in_r2c, real_prec *out_c2r, complex_prec *in_c2r, complex_prec *out_r2c) {
  struct NUMERICAL *n = data->numerical;
  n->in_r2c = in_r2c;
  n->out_c2r = out_c2r;
  n->in_c2r = in_c2r;
  n->out_r2c = out_r2c;
  n->R2Cplan = new plan_pkg(n->N1, n->N2, n->N3, in_r2c, out_r2c);
  n->C2Rplan = new plan_pkg(n->N1, n->N2, n->N3, in_c2r, out_c2r);
}


/* --- function init_cosmology --- */
int INIT_COSMOLOGY(struct COSMOLOGY *c, const std::string& codename)
{
  //struct COSMOLOGY *c = data->cosmology;

  //string codename=data->numerical->codename;
  cout<<"\n >>> starting "<<codename<<"  ...\n"; 
  cout<<endl;
  cout << " >>> "<<codename<<" initialising cosmological model...."<<endl<<endl;

  //real_prec Lx=data->numerical->L1;
  //real_prec Ly=data->numerical->L2;
  //real_prec Lz=data->numerical->L3;

  //int Nx=data->numerical->N1;
  //int Ny=data->numerical->N2;
  //int Nz=data->numerical->N3;

  //EGP  ULONG  N=Nx*Ny*Nz;

  c->omega_r = 0.0;		/* negligible radiation density */
  c->omega_k = 0.0;		/* curvature - flat prior for everything! */

  switch (cmbcosm)
  {	
    //WMAP3      
    case 1:
      {
        c->omega_m = static_cast<real_prec>(0.25);//WMPA3
        c->omega_b = static_cast<real_prec>(0.0456);
        c->omega_q = static_cast<real_prec>(1.)-c->omega_m;
        c->w = -static_cast<real_prec>(1.0);
        c->n_s = static_cast<real_prec>(1.0);//WMPA3
        c->wprime = static_cast<real_prec>(0.0);
        c->sigma8 = static_cast<real_prec>(0.9);//WMPA3
        c->rsmooth = static_cast<real_prec>(8.0);
        c->h = static_cast<real_prec>(0.73); //WMPA3
        c->beta = static_cast<real_prec>(3.0 / 2.0);
      }	
      break;
      //WMAP7  mean 
    case 2:
      {
        c->omega_m = static_cast<real_prec>(0.272);//WMPA7
        c->omega_b = static_cast<real_prec>(0.0456);
        c->omega_q = static_cast<real_prec>(1.-c->omega_m);
        c->w = static_cast<real_prec>(-1.0);
        c->n_s = static_cast<real_prec>(0.963);//WMPA7
        c->wprime = static_cast<real_prec>(0.0);
        c->sigma8 = static_cast<real_prec>(0.809);//WMPA7
        c->rsmooth = static_cast<real_prec>(8.0);
        c->h = static_cast<real_prec>(0.704);//WMPA7
        c->beta = static_cast<real_prec>(3.0 / 2.0);
      }
      break;
      //WMAP7  max
    case 3:
      {
        c->omega_m = static_cast<real_prec>(0.272);//WMPA7
        c->omega_b = static_cast<real_prec>(0.046);
        c->omega_q = static_cast<real_prec>(1.-c->omega_m);
        c->w = static_cast<real_prec>(-1.0);
        c->n_s = static_cast<real_prec>(0.961);//WMPA7
        c->wprime = static_cast<real_prec>(0.0);
        c->sigma8 = static_cast<real_prec>(0.807);//WMPA7
        c->rsmooth = static_cast<real_prec>(8.0);
        c->h = static_cast<real_prec>(0.702);//WMPA7
        c->beta = static_cast<real_prec>(3.0 / 2.0);
      }
      break;
    case 4:
      {
        //WMAP9 +eCMB+BAO+H0 best fit
        c->omega_m = static_cast<real_prec>(0.28645);//WMPA9
        c->omega_b = static_cast<real_prec>(0.04628);//WMPA9
        c->omega_q = static_cast<real_prec>(1.-c->omega_m);
        c->w = static_cast<real_prec>(-1.0);
        c->n_s = static_cast<real_prec>(0.972);//WMPA9
        c->wprime = static_cast<real_prec>(0.0);
        c->sigma8 = static_cast<real_prec>(0.82);//WMPA9
        c->rsmooth = static_cast<real_prec>(8.0);
        c->h = static_cast<real_prec>(0.6932);//WMPA9
        c->beta = static_cast<real_prec>(3.0 / 2.0);
      }
      break;
  }

  // growth factors
  real_prec Omega_C = num_1 - c->omega_m - c->omega_q;

  auto H0 = static_cast<real_prec>(100. * c->h *cgs_km/cgs_Mpc/cgs_sec);
  real_prec H = H0 * sqrt(c->omega_m/c->ascale/c->ascale/c->ascale + c->omega_q + Omega_C/c->ascale/c->ascale);
  real_prec Omega=c->omega_m / (((H/H0)*(H/H0))*(c->ascale*c->ascale*c->ascale));

  real_prec D1 = D_growth(c->ascale, c->omega_m, c->omega_q, c->h);
  c->D1 = D1;
  cout<<"D1 = "<<D1<<endl;
  auto fD2 = static_cast<real_prec>(pow(Omega,-1./143.));
  auto D2 = static_cast<real_prec>(-3./7.*D1*D1*fD2);
  c->D2=D2;
  cout<<"D2 = "<<D2<<endl;

  return(0);
}

void set_likelihood_functions(struct DATA *data)
{
  switch (data->numerical->likelihood)
  {
    case 0: // Poissonian 
      data->observational->log_like = poissonian_likelihood_log_like;
      data->observational->partial_f_delta_x_log_like = poissonian_likelihood_partial_f_delta_x_log_like;
      data->observational->grad_f_delta_x_comp = poissonian_likelihood_grad_f_delta_x_comp;
      break;
    case 1: // Gaussian 
      data->observational->log_like = gaussian_likelihood_log_like;
      data->observational->partial_f_delta_x_log_like = gaussian_likelihood_partial_f_delta_x_log_like;
      data->observational->grad_f_delta_x_comp = gaussian_likelihood_grad_f_delta_x_comp;
      break;
    case 2: // log-normal 
      data->observational->log_like = lognormal_likelihood_log_like;
      data->observational->partial_f_delta_x_log_like = lognormal_likelihood_partial_f_delta_x_log_like;
      data->observational->grad_f_delta_x_comp = lognormal_likelihood_grad_f_delta_x_comp;
      break;
    case 3: // Gaussian Random Field (only prior, effectively)
      data->observational->log_like = grf_likelihood_log_like;
      data->observational->partial_f_delta_x_log_like = grf_likelihood_partial_f_delta_x_log_like;
      data->observational->grad_f_delta_x_comp = grf_likelihood_grad_f_delta_x_comp;
      break;
  }
}

int INIT_OBSERVATIONAL(struct DATA *data,real_prec *POWER,real_prec *SIGNAL,real_prec *SIGNALX,real_prec *window,real_prec *NOISE_SF,real_prec *NOBS,real_prec *CORRF)
{
  string codename=data->numerical->codename;
  cout << " >>> "<<codename<<": initialising observational settings...."<<endl<<endl;

  data->observational->Power=POWER;		
  data->observational->signal=SIGNAL;	
  data->observational->signalX=SIGNALX;	
  data->observational->nobs=NOBS;	
  data->observational->window=window;	
  data->observational->noise_sf=NOISE_SF;	
  data->observational->corrf=CORRF;	

  data->observational->biasP=static_cast<real_prec>(1.);
  data->observational->biasE=static_cast<real_prec>(1.);
  data->observational->mu=static_cast<real_prec>(1.);
  data->observational->Nmean_Gal=static_cast<real_prec>(1.);
  data->observational->rho_c = static_cast<real_prec>(1.);

  set_likelihood_functions(data);

  switch (data->numerical->prior)
  {
    case 0: // Gaussian
      data->observational->log_prior = prior_gaussian_log_prior;
      data->observational->grad_log_prior = prior_gaussian_grad_log_prior;
      break;
  }

  return (0);
}
