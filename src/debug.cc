/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include "define_opt.h"
#include "struct_main.h"
#include "struct_hamil.h"
#include "math_funcs.h" // mean_arr, etc.
#include "IOfunctionsGen.h" // dumps
#include "IOfunctions.h" // dump_measured_spec
#include "field_statistics.h" // measure_spectrum

void debug_array_statistics(real_prec *array, ULONG size, const std::string &name) {
  cout << endl << "DEBUG -- mean " << name << ": " << mean_arr(size, array) << endl;
  cout <<         "DEBUG -- min "  << name << ": " << min_arr(size, array)  << " max " << name << ": " << max_arr(size, array) << endl;
}

void debug_scalar_dump(real_prec *array, unsigned int N1, unsigned int N2, unsigned int N3, real_prec L1, real_prec L2,
                       real_prec L3, ULONG N_bin, const std::string &fname) {
  dump_scalar(array, N1, N2, N3, fname);

  real_prec power[N_bin];
  real_prec kmode[N_bin];

  measure_spectrum(N1, N2, N3, L1, L2, L3, array, kmode, power, N_bin);
  dump_measured_spec(kmode, power, fname + std::string("_spec.dat"), N_bin);
}

void print_data(struct DATA *d) {
  using namespace std;

  struct NUMERICAL *n = d->numerical;

  cout << "Numerical\n=========" << endl;

  cout << "random_test:                       " << n->random_test << endl;
  cout << "random_test_rsd:                   " << n->random_test_rsd << endl;

  cout << "window_type:                       " << n->window_type << endl;

  cout << "data_model:                        " << n->data_model << endl;
  cout << "likelihood:                        " << n->likelihood << endl;
  cout << "prior:                             " << n->prior << endl;
  cout << "sfmodel:                           " << n->sfmodel << endl;
  cout << "rsd_model:                         " << n->rsd_model << endl;
  
  cout << "N_eps_fac:                         " << n->N_eps_fac << endl;
  cout << "eps_fac:                           " << n->eps_fac << endl;
  cout << "eps_fac_target:                    " << n->eps_fac_target << endl;
  cout << "eps_fac_initial:                   " << n->eps_fac_initial << endl;
  cout << "eps_fac_power:                     " << n->eps_fac_power << endl;
  cout << "s_eps_total_fac:                   " << n->s_eps_total_fac << endl;
  cout << "s_eps_total_scaling:               " << n->s_eps_total_scaling << endl;
  cout << "s_eps_total_Nx_norm:               " << n->s_eps_total_Nx_norm << endl;
  cout << "s_eps_total:                       " << n->s_eps_total << endl;

  // initial guess:
  cout << "initial_guess:                     " << n->initial_guess << endl;
  cout << "initial_guess_file:                " << n->initial_guess_file << endl;
  cout << "initial_guess_smoothing_type:      " << n->initial_guess_smoothing_type << endl;
  cout << "initial_guess_smoothing_scale:     " << n->initial_guess_smoothing_scale << endl;
  // mass stuff
  cout << "mass_type:                         " << n->mass_type << endl;
  cout << "massnum_init:                      " << n->massnum_init << endl;
  cout << "massnum_burn:                      " << n->massnum_burn << endl;
  // output
  cout << "outnum:                            " << n->outnum << endl;
  cout << "outnum_ps:                         " << n->outnum_ps << endl;
  // restarting
  cout << "start_at:                          " << n->start_at << endl;
  // testing:
  cout << "mass_factor:                       " << n->mass_factor << endl;
  cout << "calc_h:                            " << n->calc_h << endl;
  // EGP: end added

  cout << "seed:                              " << n->seed << endl;

  cout << "N_bin:                             " << n->N_bin << endl;

  cout << "inputmode:                         " << n->inputmode << endl;
  
  cout << "INV_SUCCESS:                       " << n->INV_SUCCESS << endl;

  cout << "filter:                            " << n->filter << endl;

  cout << "PRECON:                            " << n->PRECON << endl;

  cout << "codename:                          " << n->codename << endl;
  cout << "fnamePS:                           " << n->fnamePS << endl;
  cout << "dataFileName:                      " << n->dataFileName << endl;
  cout << "fname_ps_sample:                   " << n->fname_ps_sample << endl;
  cout << "fname_PROT_SPEC:                   " << n->fname_PROT_SPEC << endl;
  cout << "fname_PROT_SIGNAL_BOX:             " << n->fname_PROT_SIGNAL_BOX << endl;
  cout << "fname_PROT_SIGNAL_SLICE:           " << n->fname_PROT_SIGNAL_SLICE << endl;
  
  cout << "dir:                               " << n->dir << endl;

  cout << "mk:                                " << n->mk << endl;

  cout << "N1:                                " << n->N1 << endl;			/* Number of cells per axis */
  cout << "N2:                                " << n->N2 << endl;			/* Number of cells per axis */
  cout << "N3:                                " << n->N3 << endl;			/* Number of cells per axis */
  cout << "N:                                 " << n->N << endl; // total cell number
  cout << "Nhalf:                             " << n->Nhalf << endl;

  cout << "readPS:                            " << n->readPS << endl;
  
  cout << "code_norm:                         " << n->code_norm << endl;
  
  cout << "slength:                           " << n->slength << endl;			

  cout << "iGibbs:                            " << n->iGibbs << endl;
  cout << "N_Gibbs:                           " << n->N_Gibbs << endl;
  cout << "rejections:                        " << n->rejections << endl;
  cout << "total_steps_lim:                   " << n->total_steps_lim << endl;

  cout << "xllc:                              " << n->xllc << endl;			
  cout << "yllc:                              " << n->yllc << endl;			
  cout << "zllc:                              " << n->zllc << endl;			
  
  // Redshift space distortions
  cout << "xobs:                              " << n->xobs << endl;			
  cout << "yobs:                              " << n->yobs << endl;			
  cout << "zobs:                              " << n->zobs << endl;			
  cout << "planepar:                          " << n->planepar << endl;
  cout << "periodic:                          " << n->periodic << endl;
  
  cout << "L1:                                " << n->L1 << endl;			/* Length of axis */
  cout << "L2:                                " << n->L2 << endl;			/* Length of axis */
  cout << "L3:                                " << n->L3 << endl;			/* Length of axis */

  cout << "vol:                               " << n->vol << endl;
  
  cout << "d1:                                " << n->d1 << endl;			/* grid spacing x-direction */
  cout << "d2:                                " << n->d2 << endl;			/* grid spacing y-direction */
  cout << "d3:                                " << n->d3 << endl;			/* grid spacing z-direction */

  cout << "grad_psi_prior_factor:             " << n->grad_psi_prior_factor << endl;
  cout << "grad_psi_likeli_factor:            " << n->grad_psi_likeli_factor << endl;
  cout << "grad_psi_prior_conjugate:          " << n->grad_psi_prior_conjugate << endl;
  cout << "grad_psi_likeli_conjugate:         " << n->grad_psi_likeli_conjugate << endl;
  cout << "grad_psi_prior_times_i:            " << n->grad_psi_prior_times_i << endl;
  cout << "grad_psi_likeli_times_i:           " << n->grad_psi_likeli_times_i << endl;
  cout << "div_dH_by_N:                       " << n->div_dH_by_N << endl;
  cout << "deltaQ_factor:                     " << n->deltaQ_factor << endl;

  cout << "particle_kernel:                   " << n->particle_kernel << endl;
  cout << "particle_kernel_h_rel:             " << n->particle_kernel_h_rel << endl;
  cout << "particle_kernel_h:                 " << n->particle_kernel_h << endl;

  struct COSMOLOGY *c = d->cosmology;

  cout << endl << "Cosmology\n=========" << endl;

  cout << "omega_r:                           " << c->omega_r << endl;			/* radiation density */
  cout << "omega_k:                           " << c->omega_k << endl;			/* curvature density */
  cout << "omega_m:                           " << c->omega_m << endl;			/* matter density */
  cout << "omega_b:                           " << c->omega_b << endl;			/* baryon density */
  cout << "omega_q:                           " << c->omega_q << endl;			/* dark energy density */
  cout << "w:                                 " << c->w << endl;
  cout << "wprime:                            " << c->wprime << endl;		        /* de eos parameter w and first time derivative wprime=dw/da */
  cout << "sigma8:                            " << c->sigma8 << endl;			/* mass fluctuation on 8Mpc/h scales */
  cout << "rsmooth:                           " << c->rsmooth << endl;			/* smoothing radius: rsmooth = 8 yields sigma8 */
  cout << "n_s:                               " << c->n_s << endl;			/* index of primordial dm power spectrum */
  cout << "h:                                 " << c->h << endl;			        /* Hubble parameter */
  cout << "beta:                              " << c->beta << endl;			/* lens galaxy distribution parameter 1 */
  cout << "A_spec:                            " << c->A_spec << endl;			/* Normalization of the power spectrum */  
  cout << "z:                                 " << c->z << endl;  
  cout << "ascale:                            " << c->ascale << endl;
  cout << "mu:                                " << c->mu << endl;

  cout << "D1:                                " << c->D1 << endl;
  cout << "D2:                                " << c->D2 << endl;

  struct OBSERVATIONAL *o = d->observational;

  cout << endl << "Observational\n=============" << endl;

  cout << "biasP:                             " << o->biasP << endl;
  cout << "biasE:                             " << o->biasE << endl;
  cout << "Nmean_Gal:                         " << o->Nmean_Gal << endl;
  cout << "rho_c:                             " << o->rho_c << endl;
  cout << "mu:                                " << o->mu << endl;

  cout << "sigma_fac:                         " << o->sigma_fac << endl;
  cout << "sigma_min:                         " << o->sigma_min << endl;
  cout << "delta_min:                         " << o->delta_min << endl;
}

void print_hamil_data(struct HAMIL_DATA *hd) {
  using namespace std;

  cout << "Basic\n=====" << endl;

  cout << "sfmodel:                           " << hd->sfmodel << endl;
  cout << "rsd_model:                         " << hd->rsd_model << endl;
  
  // set scalars
  cout << "Nmean_Gal:                         " << hd->Nmean_Gal << endl;
  cout << "rho_c:                             " << hd->rho_c << endl;
  cout << "code_norm:                         " << hd->code_norm << endl;
  cout << "mu:                                " << hd->mu << endl;
  
  cout << "sigma_fac:                         " << hd->sigma_fac << endl;
  cout << "sigma_min:                         " << hd->sigma_min << endl;
  cout << "delta_min:                         " << hd->delta_min << endl;

  cout << "biasP:                             " << hd->biasP << endl;
  cout << "biasE:                             " << hd->biasE << endl;
  
  cout << "z:                                 " << hd->z << endl;
  cout << "ascale:                            " << hd->ascale << endl;

  cout << "D1:                                " << hd->D1 << endl;
  cout << "D2:                                " << hd->D2 << endl;
 
  cout << "OM:                                " << hd->OM << endl;
  cout << "OL:                                " << hd->OL << endl;
  cout << "h:                                 " << hd->h << endl;

  struct HAMIL_NUMERICAL *n = hd->numerical;

  cout << endl << "Numerical\n=========" << endl;

  cout << "N1:                                " << n->N1 << endl;				/* Number of cells per axis */
  cout << "N2:                                " << n->N2 << endl;				/* Number of cells per axis */
  cout << "N3:                                " << n->N3 << endl;				/* Number of cells per axis */
  cout << "N:                                 " << n->N << endl;
  cout << "Nhalf:                             " << n->Nhalf << endl;
  
  cout << "mk:                                " << n->mk << endl;				

  cout << "L1:                                " << n->L1 << endl;				
  cout << "L2:                                " << n->L2 << endl;				
  cout << "L3:                                " << n->L3 << endl;				

  cout << "vol:                               " << n->vol << endl;

  cout << "d1:                                " << n->d1 << endl;				
  cout << "d2:                                " << n->d2 << endl;				
  cout << "d3:                                " << n->d3 << endl;				

  cout << "min1:                              " << n->min1 << endl;				
  cout << "min2:                              " << n->min2 << endl;				
  cout << "min3:                              " << n->min3 << endl;				

  // Redshift space distortions
  cout << "xobs:                              " << n->xobs << endl;			
  cout << "yobs:                              " << n->yobs << endl;			
  cout << "zobs:                              " << n->zobs << endl;			
  cout << "planepar:                          " << n->planepar << endl;
  cout << "periodic:                          " << n->periodic << endl;
  
  cout << "N_bin:                             " << n->N_bin << endl;				

  cout << "iGibbs:                            " << n->iGibbs << endl;
  cout << "rejections:                        " << n->rejections << endl;
  cout << "total_steps_lim:                   " << n->total_steps_lim << endl;

  // EGP: some extra parameters for statistics on performance:
  cout << "accepted:                          " << n->accepted << endl;
  cout << "epsilon:                           " << n->epsilon << endl;
  cout << "Neps:                              " << n->Neps << endl;

  cout << "method:                            " << n->method << endl;			/* which method */
  cout << "filter:                            " << n->filter << endl;			/* which filter method */
  cout << "itmax:                             " << n->itmax << endl;			/* maximal number of iterations */
  cout << "i_reset:                           " << n->i_reset << endl;			/* number to reset */
  cout << "INV_SUCCESS:                       " << n->INV_SUCCESS << endl;		/* indicate inversion success */

  cout << "epsi:                              " << n->epsi << endl;			/* precision */
 
  cout << "kth:                               " << n->kth << endl;

  // EGP: pseudo-timestep parameters
  cout << "N_eps_fac:                         " << n->N_eps_fac << endl;
  cout << "eps_fac:                           " << n->eps_fac << endl;
  cout << "eps_fac_target:                    " << n->eps_fac_target << endl;
  cout << "eps_fac_power:                     " << n->eps_fac_power << endl;
  cout << "s_eps_total:                       " << n->s_eps_total << endl;
  
  // mass type
  cout << "mass_fs:                           " << n->mass_fs << endl;
  cout << "mass_rs:                           " << n->mass_rs << endl;
  cout << "mass_type:                         " << n->mass_type << endl;
  cout << "massnum_init:                      " << n->massnum_init << endl;
  cout << "massnum_burn:                      " << n->massnum_burn << endl;

  // EGP: testing:
  cout << "mass_factor:                       " << n->mass_factor << endl;

  /// TESTING ///
  cout << "grad_psi_prior_factor:             " << n->grad_psi_prior_factor << endl;
  cout << "grad_psi_likeli_factor:            " << n->grad_psi_likeli_factor << endl;
  cout << "grad_psi_prior_conjugate:          " << n->grad_psi_prior_conjugate << endl;
  cout << "grad_psi_likeli_conjugate:         " << n->grad_psi_likeli_conjugate << endl;
  cout << "grad_psi_prior_times_i:            " << n->grad_psi_prior_times_i << endl;
  cout << "grad_psi_likeli_times_i:           " << n->grad_psi_likeli_times_i << endl;
  cout << "div_dH_by_N:                       " << n->div_dH_by_N << endl;
  cout << "calc_h:                            " << n->calc_h << endl;
  cout << "deltaQ_factor:                     " << n->deltaQ_factor << endl;

  cout << "particle_kernel:                   " << n->particle_kernel << endl;
  cout << "particle_kernel_h:                 " << n->particle_kernel_h << endl;
}
