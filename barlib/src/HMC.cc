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
#include <cassert>
#include <iostream>

#include "fftw_array.h"

#include "convolution.hpp" // convolve
#include "IOfunctionsGen.h"

#include "HMC_help.h"
#include "HMC_momenta.h"
#include "HMC_mass.h"
#include "hmc/likelihood/gaussian_random_field.hpp"
#include "hmc/leapfrog/time_step.hpp" // update_eps_fac, update_epsilon_acc_rate_tables

#include "convenience.h"

#include "debug.h"

using namespace std;

#ifdef MASKING
bool masking = true;
#else
bool masking = false;
#endif  // ifdef MASKING


void write_to_performance_log(struct DATA *data, struct HAMIL_DATA *hd) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  ofstream &plog = data->numerical->performance_log;

  assert(plog.is_open());

  plog << n->accepted << "\t";
  plog << n->epsilon << "\t";
  plog << n->Neps << "\t";
  plog << n->dH << "\t";
  plog << n->dK << "\t";
  plog << n->dE << "\t";
  plog << n->dprior << "\t";
  plog << n->dlikeli << "\t";
  plog << n->psi_prior_i << "\t";
  plog << n->psi_prior_f << "\t";
  plog << n->psi_likeli_i << "\t";
  plog << n->psi_likeli_f << "\t";
  plog << n->H_kin_i << "\t";
  plog << n->H_kin_f << endl;
}


// for Gaussian momentum distribution //
real_prec kinetic_term(struct HAMIL_DATA *hd, real_prec *momenta, struct DATA *data) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  // H_kin = 1/2 * p * M^-1 * p
  fftw_array<real_prec> dummy(n->N);

  if (n->mass_fs) {
    if (masking == true) {
      // X^{-1}S^(1/2)p
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
      for (ULONG i = 0; i < n->N; i++) {
        real_prec corr = hd->corrf[i];
        real_prec winprime = num_1 + corr;
        dummy[i] = momenta[i] / winprime;
      }
      convolveInvCorrFuncWithSignal(hd, dummy, dummy, hd->mass_f);
    } else {
      convolveInvCorrFuncWithSignal(hd, momenta, dummy, hd->mass_f);
    }
  } else {
    fillZero(dummy, n->N);  // prepare for addition to real part
  }

  if (n->mass_rs) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
    for (ULONG i = 0; i < n->N; i++) {
      real_prec invM = 0.;
      if (hd->mass_r[i] > 0.0)
        invM = static_cast<real_prec>(1. / hd->mass_r[i]);

      dummy[i] += invM * momenta[i];
    }
  }

  // 5) 1/2*p*IFT[ 1/M*FT[p] ]

  bool flag_write = false;
  real_prec value = 0.;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:value)
#endif  // MULTITHREAD
  for (ULONG i = 0; i < n->N; i++) {
    value += num_0_5*momenta[i]*dummy[i];

    if (momenta[i] > 1e100 && flag_write == false) {
      cout << "aaargggh! too much kinetic energy!" << endl;
      flag_write = true;
    }
  }
  // cout << endl;
  // cout << "... kinetic energy: " << value << endl;
  wprintw(data->curses->table, "%5.0e ", value);

  return(value);
}


real_prec psi(struct HAMIL_DATA *hd, real_prec *signal, struct DATA *data) {
  // Psi = 1/2* s_k * S^-1 * s_k -log L
  // The "potential energy" analogue of the Hamiltonian MC algorithm.
  struct HAMIL_NUMERICAL *n = hd->numerical;

  // real_prec psi_prior = prior_gaussian_log_prior(hd, signal);
  real_prec psi_prior = hd->log_prior(hd, signal);
  real_prec psi_likelihood = hd->log_like(hd, signal);  // care
  real_prec psi_total = psi_prior + psi_likelihood;

  // cout << endl << "... energy of the prior: " << psi_prior << " energy of the likelihood: " << psi_likelihood << endl;
  wprintw(data->curses->table, "%5.0e ", psi_prior);
  wprintw(data->curses->table, "%5.0e ", psi_likelihood);

  // for performance log:
  n->psi_prior = psi_prior;
  n->psi_likeli = psi_likelihood;

  return(psi_total);
}


void gradient_psi(struct HAMIL_DATA *hd, real_prec *signal, struct DATA *data) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  fftw_array<real_prec> grad_psi_prior(n->N), grad_psi_likeli(n->N);

  // prior_gaussian_grad_log_prior(hd, signal, grad_psi_prior);
  hd->grad_log_prior(hd, signal, grad_psi_prior);

#ifdef DEBUG
  debug_array_statistics(grad_psi_prior, n->N, "gradPsi prior");
#endif // DEBUG

  // IFT[ 1/P*FT[s] ] + dPsiL/ds //

  if (data->numerical->likelihood == 3)  // Gaussian random field
    grf_likelihood_grad_log_like(hd, signal, grad_psi_likeli);
  else
    likelihood_grad_log_like(hd, signal, grad_psi_likeli);

#ifdef DEBUG
  debug_array_statistics(grad_psi_likeli, n->N, "gradPsi likelihood");
#endif // DEBUG


  //// TESTING ////
  multiply_factor_array(n->grad_psi_prior_factor, grad_psi_prior,
                        grad_psi_prior, n->N);
  multiply_factor_array(n->grad_psi_likeli_factor, grad_psi_likeli,
                        grad_psi_likeli, n->N);
  if (n->grad_psi_prior_conjugate || n->grad_psi_prior_times_i) {
    // fftw_array<complex_prec> AUX(n->N), dummyC(n->N);
    // complexify_array(grad_psi_prior, dummyC, n->N);
    // FFT3d(n->N1, n->N2, n->N3, true, dummyC, AUX);
    fftw_array<complex_prec> grad_psi_prior_C(n->Nhalf);
    fftR2C(n->N1, n->N2, n->N3, grad_psi_prior, grad_psi_prior_C);
    if (n->grad_psi_prior_conjugate)
      conjugate_array(grad_psi_prior_C, grad_psi_prior_C, n->Nhalf);
    if (n->grad_psi_prior_times_i)
      times_i_array(grad_psi_prior_C, grad_psi_prior_C, n->Nhalf);
    // FFT3d(n->N1, n->N2, n->N3, false, AUX, dummyC);
    // real_part_array(dummyC, grad_psi_prior, n->N);
    fftC2R(n->N1, n->N2, n->N3, grad_psi_prior_C, grad_psi_prior);
  }
  if (n->grad_psi_likeli_conjugate || n->grad_psi_likeli_times_i) {
    // fftw_array<complex_prec> AUX(n->N), dummyC(n->N);
    // complexify_array(grad_psi_likeli, dummyC, n->N);
    // FFT3d(n->N1, n->N2, n->N3, true, dummyC, AUX);
    fftw_array<complex_prec> grad_psi_likeli_C(n->Nhalf);
    fftR2C(n->N1, n->N2, n->N3, grad_psi_likeli, grad_psi_likeli_C);
    if (n->grad_psi_likeli_conjugate)
      conjugate_array(grad_psi_likeli_C, grad_psi_likeli_C, n->N);
    if (n->grad_psi_likeli_times_i)
      times_i_array(grad_psi_likeli_C, grad_psi_likeli_C, n->N);
    // FFT3d(n->N1, n->N2, n->N3, false, AUX, dummyC);
    // real_part_array(dummyC, grad_psi_likeli, n->N);
    fftC2R(n->N1, n->N2, n->N3, grad_psi_likeli_C, grad_psi_likeli);
  }
  //// END TESTING ////


  sum_arrays(grad_psi_prior, grad_psi_likeli, hd->gradpsi, n->N);
}


real_prec delta_Hamiltonian(struct HAMIL_DATA *hd, real_prec *signali, real_prec *momentai, real_prec *signalf,
                            real_prec *momentaf, struct DATA *data) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  // dH = H(sf,pf) - H(si,pi)

  real_prec Hkini = kinetic_term(hd, momentai, data);
  real_prec Hpsii = psi(hd, signali, data);

  // save initial psi's for performance log
  n->psi_prior_i = n->psi_prior;
  n->psi_likeli_i = n->psi_likeli;
  n->H_kin_i = Hkini;

  real_prec Hami = Hkini + Hpsii;

  real_prec Hkinf = kinetic_term(hd, momentaf, data);
  real_prec Hpsif = psi(hd, signalf, data);

  // more performance log stuff:
  n->dprior = n->psi_prior - n->psi_prior_i;
  n->dlikeli = n->psi_likeli - n->psi_likeli_i;

  real_prec Hamf = Hkinf + Hpsif;

  real_prec dHam = Hamf - Hami;
  if (n->div_dH_by_N) {
    // test: divide dH by N to scale for cell number & make the code accept more
    dHam /= static_cast<real_prec>(n->N);
  }

  // save final stuff in hd for output to performance log:
  n->dH = dHam;
  n->dK = Hkinf - Hkini;
  n->dE = Hpsif - Hpsii;
  n->psi_prior_f = n->psi_prior;
  n->psi_likeli_f = n->psi_likeli;
  n->H_kin_f = Hkinf;

  return(dHam);
}


void Hamiltonian_EoM(struct HAMIL_DATA *hd, real_prec *signali,
                     real_prec *momentai, real_prec *signalf,
                     real_prec *momentaf, gsl_rng * seed, struct DATA *data) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  // tau=n*e
  // 1) pi(t+e/2)=pi(t)-e/2*gradPsi(si(t))
  // 2) si(t+e)=si(t)+e/Mi*pi(t+e/2)
  // 3) pi(t+e)=pi(t+e/2)-e/2*gradPsi(si(t+e))

  n->Neps = static_cast<ULONG>(n->N_eps_fac*(gsl_rng_uniform(seed))) + 1;
  n->epsilon = static_cast<real_prec>(n->eps_fac*gsl_rng_uniform(seed));

  if (n->epsilon > 2.)
    n->epsilon = 2.;

  // cout << endl;
  // cout << "----> Number of leapfrog iterations: " << n->Neps <<endl;
  // cout << "----> Step size of the leapfrog scheme: " << n->epsilon <<endl;
  wprintw(data->curses->table, "%5.0e ", n->epsilon);
  wprintw(data->curses->table, "%4lu ", n->Neps);
  wrefresh(data->curses->table);

  fftw_array<real_prec> dummy(n->N);

  copyArray(signali, signalf, n->N);
  copyArray(momentai, momentaf, n->N);

  // 0) determine gradient_psi at t = 0
  fillZero(hd->gradpsi, n->N);
  gradient_psi(hd, signalf, data);

  // cout << "... Leaping " << flush;

  for (ULONG jj = 0; jj < n->Neps; jj++) {
    wprintw(data->curses->status, "\nLeap-frogging ... %lu/%lu", jj+1, n->Neps);
    wrefresh(data->curses->status);
    // cout << "." << flush;

    // 1) pi(t+e/2)=pi(t)-e/2*gradPsi(si(t))
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
    for (ULONG i = 0; i < n->N; i++)
      momentaf[i] -= num_0_5 * n->epsilon * hd->gradpsi[i];


    // 2a) M^-1 * p(t+e/2)
    if (n->mass_fs) {
#ifdef MASKING
      {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
        for (ULONG i = 0; i < n->N; i++)
          dummy[i]=momentaf[i]/(num_1 + hd->corrf[i]);
        convolveInvCorrFuncWithSignal(hd, dummy, dummy, hd->mass_f);
      }
#else  // MASKING
      {
        convolveInvCorrFuncWithSignal(hd, momentaf, dummy, hd->mass_f);
      }
#endif  // MASKING
    } else {  // mass_fs == false
      fillZero(dummy, n->N);
    }

    if (n->mass_rs) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
      for (ULONG i = 0; i < n->N; i++) {
        real_prec invM = 0.;
        if (hd->mass_r[i] > 0.0)
          invM = num_1/hd->mass_r[i];
        dummy[i] += momentaf[i]*invM;
      }
    }

#ifdef DEBUG
    debug_array_statistics(dummy, n->N, "shift / epsilon");
#endif // DEBUG


    // 2b) si(t+e)=si(t)+e/Mi*pi(t+e/2)
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
    for (ULONG i = 0; i < n->N; i++)
      signalf[i] += n->epsilon*dummy[i];


    // 3a) determine gradient_psi at t+e
    fillZero(hd->gradpsi, n->N);
    gradient_psi(hd, signalf, data);


    // 3b) pi(t+e)=pi(t+e/2)-e/2*gradPsi(si(t+e))
#ifdef MULTITHREAD
#pragma omp parallel for
#endif  // MULTITHREAD
    for (ULONG i = 0; i < n->N; i++)
      momentaf[i] -= num_0_5 * n->epsilon * hd->gradpsi[i];

#ifdef DEBUG
    debug_array_statistics(momentaf, n->N, "momf");
#endif // DEBUG

    // check for too high values of momentum; it sometimes goes crazy out of
    // control and at ~1e308 a SIGFPE stops the program
    if (abs(momentaf[0]) > 1e50) {
      wprintw(data->curses->message, "\nLeap-frogging ... stopped at %lu/%lu, momentum too high (momenta[0] = %e)", jj+1, n->Neps, momentaf[0]);
      wrefresh(data->curses->message);
      jj = n->Neps;  // stop the loop, this EoM will lead to nothing
    }
  }
  // cout<<endl;
  // cout <<"----> maximum value in the reconstruction: "<<max_arr(n->N,signalf)<<endl;
  data->numerical->count_attempts++;
}


void HamiltonianMC(struct HAMIL_DATA *hd, gsl_rng * seed, struct DATA *data) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  struct NUMERICAL *dn = data->numerical;
  fftw_array<real_prec> momentai(n->N), momentaf(n->N), signali(n->N),
                        signalf(n->N);

  bool reach_acceptance = false;

  // define the starting signal estimate //
  ULONG iter = 0;

  // real_prec maxs = max_arr(n->N, hd->x);
  // cout << "maximum of the guess signal: "<<maxs<<endl;

  // calculate Hamiltonian masses //
  ULONG massnum;
  if (n->iGibbs > n->massnum_burn)
    massnum = n->massnum_burn;
  else
    massnum = n->massnum_init;

  // if (massnum == 1)
  // {
  //   cout<<"---->compute mass for the momenta ... "<<endl;
  //   Hamiltonian_mass(hd,hd->x,seed, data);
  // }
  // else
  // {
  if (0 == n->iGibbs % massnum || n->iGibbs == 1) {
    // cout << "---->compute mass for the momenta ... " << endl;
    Hamiltonian_mass(hd, hd->x, data);
    string fname = dn->dir + string("auxmass_r");
    if (n->mass_rs) {
      if (contains_nan(hd->mass_r, n->N)) {
        throw runtime_error("auxmass_r contains a NaN! aborting.");
      }
      write_array(fname, hd->mass_r, n->N1, n->N2, n->N3);
    }
    fname = dn->dir + string("auxmass_f");
    if (n->mass_fs) {
      write_array(fname, hd->mass_f, n->N1, n->N2, n->N3);
    }
  } else {
    string fname = dn->dir + string("auxmass_r");
    if (n->mass_rs) {
      read_array(fname, hd->mass_r, n->N1, n->N2, n->N3);
    }
    fname = dn->dir + string("auxmass_f");
    if (n->mass_fs) {
      read_array(fname, hd->mass_f, n->N1, n->N2, n->N3);
    }
  }
  // }

  // cout << "---->start Hamiltonian sampling ... " << endl;
  wprintw(data->curses->status, "starting Hamiltonian sampling");
  wrefresh(data->curses->status);

  int naccepted = 0;
  while (iter < n->itmax && reach_acceptance == false) {
    // print sample number
    wprintw(data->curses->table, "%6lu ", n->iGibbs);
    wrefresh(data->curses->table);

    iter++;
    // cout<<endl
    // cout << "----> Number of Metropolis-Hastings sample candidate: " << iter <<endl;
    // print sample-candidate number
    wprintw(data->curses->table, "%4lu ", iter);
    wrefresh(data->curses->table);

    // 0) initiliaze
    // define the starting position
    copyArray(hd->x, signali, n->N);

    // 1) step 1: new momenta
    // draw momenta
    draw_momenta(hd, seed, momentai, data);
    // quick_dump_scalar(momentai, n->N1, "initial_momentum", iter, true);

    // 2) step 2: Eq of motion
    update_eps_fac(hd, data);
    // calculate next positions with the Equation of motions
    Hamiltonian_EoM(hd, signali, momentai, signalf, momentaf, seed, data);

    // 3) step 3: acceptance of the new position
    // calculate the Hamiltonian difference
    real_prec dH = delta_Hamiltonian(hd, signali, momentai, signalf, momentaf, data);

    // calculate the acceptance probability
    real_prec p_acceptance = 1.;

    if (dH < 0.0)
      p_acceptance = 1.0;
    else if (exp(-dH) < 1.0)
      p_acceptance = exp(-dH);

    // cout<<endl;
    // cout<<"----> ACHTUNG dH: "<< dH <<endl;
    // cout<<"----> ACHTUNG acceptance probability: "<< p_acceptance <<endl;
    wprintw(data->curses->table, "%6.0e ", dH);
    wprintw(data->curses->table, "%5.0e ", p_acceptance);
    wrefresh(data->curses->table);

    // draw acceptance
    if (p_acceptance >= 1.0) {
      reach_acceptance = true;
    } else {
      auto u = static_cast<real_prec>(gsl_rng_uniform(seed));

      if (u < p_acceptance)
        reach_acceptance = true;
      else
        reach_acceptance = false;
    }

    if (reach_acceptance == true) {
      naccepted++;
      // cout<<"----> accepted! ;-) # "<<naccepted<<endl;
      wprintw(data->curses->table, "y ");
    } else {
      wprintw(data->curses->table, "n ");
    }
    wrefresh(data->curses->table);

    if (reach_acceptance)
      copyArray(signalf, hd->x, n->N);

    if (!reach_acceptance)
      n->rejections++;

    // EGP: save in hd for performance log:
    n->accepted = reach_acceptance;

    write_to_performance_log(data, hd);
    update_epsilon_acc_rate_tables(data, hd);

    // ULONG total_steps = n->iGibbs + n->rejections + dn->rejections;
    ULONG total_steps;
    if (reach_acceptance) {
      total_steps = n->iGibbs + n->rejections + dn->rejections + 1;
    } else {
      total_steps = n->iGibbs + n->rejections + dn->rejections;
    }

    // compute recent acceptance rate
    // -1: index, not number
    ULONG ix_acc = (dn->count_attempts-1) % dn->N_a_eps_update;
    if (reach_acceptance) {
      dn->acc_recent[ix_acc] = 1;
    } else {
      dn->acc_recent[ix_acc] = 0;
    }
    real_prec acc_N_a = 0;
    for (unsigned i = 0; i < dn->N_a_eps_update; ++i) {
      acc_N_a += static_cast<real_prec>(dn->acc_recent[i]);
    }
    acc_N_a /= static_cast<real_prec>(dn->N_a_eps_update);
    wprintw(data->curses->table, "%4.2f", acc_N_a);
    wrefresh(data->curses->table);

    if (total_steps >= n->total_steps_lim) {
      throw runtime_error("ABORTING: total steps exceeds total_steps_lim");
    }

    wprintw(data->curses->table, "\n");
    wrefresh(data->curses->table);
  }

  // set inversion success
  if (reach_acceptance == true)
    n->INV_SUCCESS = 1;
  else
    n->INV_SUCCESS = 0;

  // cout<<endl;
}
