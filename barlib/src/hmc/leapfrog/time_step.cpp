/*
 * Barcode
 * Copyright E.G.P. Bos
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <ncurses.h>

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>  // find_if, sort, count_if, max_element, min_element, copy
#include <iomanip> // setw

#include "struct_main.h"
#include "struct_hamil.h"
#include "math_funcs.h" // power_mean

#include "hmc/leapfrog/time_step.hpp"


real_prec bool_mean(const std::vector<bool>& input) {
  auto result = static_cast<real_prec>(std::count_if(input.begin(), input.end(), [](bool i){return i;}));
  result /= static_cast<real_prec>(input.size());
  return result;
}

std::vector<real_prec> real_from_bool(const std::vector<bool>& input) {
  std::vector<real_prec> output;
  output.reserve(input.size());
  for (auto i : input) {
    output.push_back(static_cast<real_prec>(i));
  }
  return output;
}


std::string update_eps_fac_acceptance_rate_downwards(struct HAMIL_DATA *hd,
                                                     struct DATA *data) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  struct NUMERICAL *dn = data->numerical;

  real_prec alpha_N_a = bool_mean(dn->acc_flag_N_a);
  real_prec acc_target = (dn->acc_max + dn->acc_min)/2.;

  std::string message = "\nadjusted eps_fac downwards to %f";
  // make vector of acceptance flags sorted by epsilon:
  auto a_sort_bool = sort_vector_by_other(dn->acc_flag_N_a, dn->epsilon_N_a);
  auto a_sort = real_from_bool(a_sort_bool);

  // turn it into cumulative moving average of acceptance rate (sort by eps)
  cumulative_moving_average(a_sort.begin(), a_sort.end(), a_sort.begin());

  // smooth it a bit, due to the discrete acc_flags it will be too bumpy
  std::vector<real_prec> a_sm(dn->N_a_eps_update);
  stl_smooth(a_sort.begin(), a_sort.end(), a_sm.begin(), dn->eps_down_smooth);

  // loop until an element larger than acc_target is found, just to be sure we
  // start at a peak; in some cases there might be zeroes at the beginning of
  // the acceptance_rate array, which will cause an initial dip, which must not
  // be mistaken for the true location where the line goes below the target
  // acceptance rate
  auto ix_max_a_sm = std::max_element(a_sm.begin(), a_sm.end());
  // check if indeed is larger than acc_target (if not, we can't continue):
  if (*ix_max_a_sm > acc_target) {
    // then find where the line goes below the target acceptance rate.
    auto ix_target = find_if(ix_max_a_sm, a_sm.end(), [&] (real_prec acc_sm) {
      return (acc_sm < acc_target);
    });
    // if acc_target not crossed, we have some strange special case:
    if (ix_target == a_sm.end()) {
      std::stringstream ss;
      ss << "\neps_fac stays at %f (special: alpha_N_a=" << alpha_N_a <<
         ", a_t=" << acc_target << ", a_sm_max=" << *ix_max_a_sm << ")";
      message = ss.str();
    } else {
      // adjust downward the way we wanted
      std::vector<real_prec> eps_sort(dn->N_a_eps_update);
      std::copy(dn->epsilon_N_a.begin(), dn->epsilon_N_a.end(), eps_sort.begin());
      std::sort(eps_sort.begin(), eps_sort.end());
      unsigned ix_eps = static_cast<unsigned>(ix_target - a_sm.begin());
      n->eps_fac = eps_sort[ix_eps];
    }
  } else {
    // Can't continue, so we just take a guess for the next epsilon.
    if (alpha_N_a == 0.) {
      // If there was no accepted step, set to the lowest epsilon we tried:
      auto ix_min_eps = std::min_element(dn->epsilon_N_a.begin(),
                                         dn->epsilon_N_a.end());
      n->eps_fac = *ix_min_eps;
    } else {
      // Otherwise, just divide by three (not two, that's not enough most
      // of the time that you can continue...)
      n->eps_fac /= 3.;
    }
  }
  if (n->eps_fac == 0.) {
    throw std::runtime_error("In update_eps_fac_acceptance_rate_downwards: "
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
    std::string message;  // to send to message window at the end
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
    std::string message = "\nadjusted eps_fac downwards to %f";
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
        std::cout << "  updating eps_fac to " << n->eps_fac << std::endl;
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
  std::stringstream s;
  s << std::setw(9) << "\n" << data->numerical->epsilon_N_a[0];
  for (unsigned i = 1; (i < 8) && (i < data->numerical->N_a_eps_update); ++i) {
    s << " " << data->numerical->epsilon_N_a[i];
  }
  wprintw(data->curses->debug, s.str().c_str());
  wrefresh(data->curses->debug);
}
