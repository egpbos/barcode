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
#include <cassert>
#include <numeric>  // accumulate, partial_sum
#include <algorithm>  // find, find_if, copy, count_if
#include <iostream>

#include "fftw_array.h"

#include "math_funcs.h" // max_arr, power_mean
#include "convolution.hpp" // convolve
#include "IOfunctionsGen.h"

#include "HMC_help.h"
#include "HMC_momenta.h"
#include "HMC_mass.h"
#include "hmc/likelihood/gaussian_random_field.hpp"

#include "convenience.h"

#include "debug.h"

using namespace std;

#ifdef MASKING
bool masking = true;
#else
bool masking = false;
#endif  // ifdef MASKING


void PROTOCOL_HMC(real_prec dH, real_prec Pa, struct DATA *data) {
  string outputFileName= data->numerical->dir + string("HMC.prt");

  ofstream outStream(outputFileName.data());
  assert(outStream.is_open());

  outStream << dH  <<endl;
  outStream << Pa  <<endl;
  outStream.close();
}


void PROTOCOL_NREJ(ULONG nrej, struct DATA *data) {
  string outputFileName= data->numerical->dir + string("NREJ.prt");

  ofstream outStream(outputFileName.data());
  assert(outStream.is_open());

  outStream << nrej-1  <<endl;
  outStream.close();
}



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

// TODO: MOVE THIS TO TEMPLATE/HEADER FILE
// From stackoverflow.com/a/12399290
template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {
  // initialize original index locations
  vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

template <typename T, typename T_other>
vector<T> sort_vector_by_other(const vector<T> &sortee,
                               const vector<T_other> &other) {
  if (sortee.size() != other.size()) {
    throw runtime_error("sort_vector_by_other: vector sizes must be equal!");
  }
  auto ix_sort = sort_indexes(other);
  vector<T> sorted;
  sorted.reserve(ix_sort.size());
  for (auto i: ix_sort) {
    sorted.push_back(sortee[i]);
  }
  return (sorted);
}

template <typename InputIterator, typename OutputIterator>
OutputIterator cumulative_moving_average(InputIterator begin, InputIterator end,
                                         OutputIterator result) {
  partial_sum(begin, end, result);
  int i = 1;  // start at one, not zero (number, not index)
  while (++begin != end) {
    ++result;
    *result = *result / static_cast<real_prec>(++i);  // average up to that ix
  }
  return ++result;
}

template <typename InputIterator, typename OutputIterator>
OutputIterator stl_smooth(InputIterator begin, InputIterator end,
                          OutputIterator result, int smooth_size) {
  for (auto ix = begin; ix != end; ++ix) {
    auto ix_min = ix - smooth_size;
    auto ix_max = ix + smooth_size + 1;  // +1: accumulate is not inclusive last
    if (ix_min < begin) ix_min = begin;
    if (ix_max > end)   ix_max = end;
    *result++ = accumulate(ix_min, ix_max, 0.0)
                / static_cast<real_prec>(ix_max - ix_min);
  }
  return result;
}

real_prec bool_mean(const vector<bool>& input) {
  auto result = static_cast<real_prec>(std::count_if(input.begin(), input.end(), [](bool i){return i;}));
  result /= static_cast<real_prec>(input.size());
  return result;
}

vector<real_prec> real_from_bool(const vector<bool>& input) {
  vector<real_prec> output;
  output.reserve(input.size());
  for (auto i : input) {
    output.push_back(static_cast<real_prec>(i));
  }
  return output;
}

string update_eps_fac_acceptance_rate_downwards(struct HAMIL_DATA *hd,
                                                struct DATA *data) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  struct NUMERICAL *dn = data->numerical;

  real_prec alpha_N_a = bool_mean(dn->acc_flag_N_a);
  real_prec acc_target = (dn->acc_max + dn->acc_min)/2.;

  string message = "\nadjusted eps_fac downwards to %f";
  // make vector of acceptance flags sorted by epsilon:
  auto a_sort_bool = sort_vector_by_other(dn->acc_flag_N_a, dn->epsilon_N_a);
  auto a_sort = real_from_bool(a_sort_bool);

  // turn it into cumulative moving average of acceptance rate (sort by eps)
  cumulative_moving_average(a_sort.begin(), a_sort.end(), a_sort.begin());

  // smooth it a bit, due to the discrete acc_flags it will be too bumpy
  vector<real_prec> a_sm(dn->N_a_eps_update);
  stl_smooth(a_sort.begin(), a_sort.end(), a_sm.begin(), dn->eps_down_smooth);

  // loop until an element larger than acc_target is found, just to be sure we
  // start at a peak; in some cases there might be zeroes at the beginning of
  // the acceptance_rate array, which will cause an initial dip, which must not
  // be mistaken for the true location where the line goes below the target
  // acceptance rate
  auto ix_max_a_sm = max_element(a_sm.begin(), a_sm.end());
  // check if indeed is larger than acc_target (if not, we can't continue):
  if (*ix_max_a_sm > acc_target) {
    // then find where the line goes below the target acceptance rate.
    auto ix_target = find_if(ix_max_a_sm, a_sm.end(), [&] (real_prec acc_sm) {
                              return (acc_sm < acc_target);
                            });
    // if acc_target not crossed, we have some strange special case:
    if (ix_target == a_sm.end()) {
      stringstream ss, debugss;
      ss << "\neps_fac stays at %f (special: alpha_N_a=" << alpha_N_a <<
            ", a_t=" << acc_target << ", a_sm_max=" << *ix_max_a_sm << ")";
      message = ss.str();
    } else {
      // adjust downward the way we wanted
      vector<real_prec> eps_sort(dn->N_a_eps_update);
      copy(dn->epsilon_N_a.begin(), dn->epsilon_N_a.end(), eps_sort.begin());
      sort(eps_sort.begin(), eps_sort.end());
      unsigned ix_eps = static_cast<unsigned>(ix_target - a_sm.begin());
      n->eps_fac = eps_sort[ix_eps];
    }
  } else {
    // Can't continue, so we just take a guess for the next epsilon.
    if (alpha_N_a == 0.) {
      // If there was no accepted step, set to the lowest epsilon we tried:
      auto ix_min_eps = min_element(dn->epsilon_N_a.begin(),
                                    dn->epsilon_N_a.end());
      n->eps_fac = *ix_min_eps;
    } else {
      // Otherwise, just divide by three (not two, that's not enough most
      // of the time that you can continue...)
      n->eps_fac /= 3.;
    }
  }
  if (n->eps_fac == 0.) {
    throw runtime_error("In update_eps_fac_acceptance_rate_downwards: "
                           "epsilon became zero, shouldn't happen!");
  }
  return (message);
}

void update_eps_fac_acceptance_rate(struct HAMIL_DATA *hd, struct DATA *data) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  struct NUMERICAL *dn = data->numerical;
  // ULONG total_steps = n->iGibbs + n->rejections + data->numerical->rejections;
  // modulo (total - 1), because this is the step we haven't finished yet
  // if (((total_steps - 1)% dn->N_a_eps_update == 0) && (total_steps > 1)) {

  // N.B.: the above (with -1) was only for total_steps, not for count_attempts!
  if ((dn->count_attempts % dn->N_a_eps_update == 0)
      && (dn->count_attempts > 0)) {
    string message;  // to send to message window at the end
    // recent (N_a tries) acceptance rate
    real_prec alpha_N_a = bool_mean(dn->acc_flag_N_a);

    if (alpha_N_a < dn->acc_min) {
      // adjust downwards
      message = update_eps_fac_acceptance_rate_downwards(hd, data);
    } else if (alpha_N_a > dn->acc_max) {
      // adjust upwards
      real_prec acc_target = (dn->acc_max + dn->acc_min)/2.;
      n->eps_fac *= dn->eps_up_fac * (alpha_N_a / acc_target);
      message = "\nadjusted eps_fac upwards to %f";
    } else {
      // don't adjust
      message = "\nnot adjusting eps_fac, stays at %f";
    }
    wprintw(data->curses->message, message.c_str(), n->eps_fac);
    wrefresh(data->curses->message);
  }
}

void update_eps_fac_acceptance_rate_fast_initial(struct HAMIL_DATA *hd,
                                                 struct DATA *data) {
  struct HAMIL_NUMERICAL *n = hd->numerical;

  if ((n->iGibbs == 1) && (n->rejections > 0)) {
    n->eps_fac /= 2.;
    string message = "\nadjusted eps_fac downwards to %f";
    wprintw(data->curses->message, message.c_str(), n->eps_fac);
    wrefresh(data->curses->message);
  } else {
    update_eps_fac_acceptance_rate(hd, data);
  }
}

void update_eps_fac(struct HAMIL_DATA *hd, struct DATA *data) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  struct NUMERICAL *dn = data->numerical;
  switch (data->numerical->eps_fac_update_type) {
    case 0:
      // don't update
      break;
    case 1:
    {
      // update eps_fac every s_eps_total steps
      // ULONG total_steps = n->iGibbs + n->rejections
                          // + data->numerical->rejections;
      // modulo (total - 1), because this is the step we haven't finished yet
      // if (((total_steps-1) % n->s_eps_total == 0) && (total_steps > 1)) {

      // N.B.: the above (-1) was only for total_steps, not for count_attempts!
      if ((dn->count_attempts % n->s_eps_total == 0)
          && (dn->count_attempts > 0)) {
        n->eps_fac = power_mean(n->eps_fac, n->eps_fac_target,
                                n->eps_fac_power);
        cout << "  updating eps_fac to " << n->eps_fac << endl;
      }
      break;
    }
    case 2:
      // adjust eps_fac based on recent acceptance rate
      update_eps_fac_acceptance_rate(hd, data);
      break;
    case 3:
      // adjust eps_fac based on recent acceptance rate
      // version 2: with extra fast initial phase until first accepted step
      update_eps_fac_acceptance_rate_fast_initial(hd, data);
      break;
  }
}

void update_epsilon_acc_rate_tables(struct DATA *data, struct HAMIL_DATA *hd) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  struct NUMERICAL *dn = data->numerical;
  // ULONG total_steps = n->iGibbs + n->rejections + data->numerical->rejections;
  // ULONG ix_table = (total_steps-1)  // -1: index, not number
  ULONG ix_table = (dn->count_attempts-1)  // -1: index, not number
                   % data->numerical->N_a_eps_update;
  data->numerical->acc_flag_N_a[ix_table] = n->accepted;
  data->numerical->epsilon_N_a[ix_table] = n->epsilon;
  stringstream s;
  s << setw(9) << "\n" << data->numerical->epsilon_N_a[0];
  for (unsigned i = 1; (i < 8) && (i < data->numerical->N_a_eps_update); ++i) {
    s << " " << data->numerical->epsilon_N_a[i];
  }
  wprintw(data->curses->debug, s.str().c_str());
  wrefresh(data->curses->debug);
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

#ifdef DEBUG
    PROTOCOL_HMC(dH, p_acceptance, data);
#endif  // DEBUG

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
#ifdef DEBUG
  PROTOCOL_NREJ(iter, data);
#endif  // DEBUG

  // set inversion success
  if (reach_acceptance == true)
    n->INV_SUCCESS = 1;
  else
    n->INV_SUCCESS = 0;

  // cout<<endl;
}
