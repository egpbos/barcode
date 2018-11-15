/*
 * Barcode
 * Copyright E.G.P. Bos
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#ifndef BARCODE_TIME_STEP_HPP
#define BARCODE_TIME_STEP_HPP

#include <stdexcept> // runtime_error
#include <vector>
#include <algorithm> // sort
#include <numeric>  // accumulate, partial_sum

#include "define_opt.h"

void update_eps_fac(struct HAMIL_DATA *hd, struct DATA *data);
void update_epsilon_acc_rate_tables(struct DATA *data, struct HAMIL_DATA *hd);

// From stackoverflow.com/a/12399290
template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

template <typename T, typename T_other>
std::vector<T> sort_vector_by_other(const std::vector<T> &sortee,
                                    const std::vector<T_other> &other) {
  if (sortee.size() != other.size()) {
    throw std::runtime_error("sort_vector_by_other: vector sizes must be equal!");
  }
  auto ix_sort = sort_indexes(other);
  std::vector<T> sorted;
  sorted.reserve(ix_sort.size());
  for (auto i: ix_sort) {
    sorted.push_back(sortee[i]);
  }
  return (sorted);
}

template <typename InputIterator, typename OutputIterator>
OutputIterator cumulative_moving_average(InputIterator begin, InputIterator end,
                                         OutputIterator result) {
  std::partial_sum(begin, end, result);
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
    *result++ = std::accumulate(ix_min, ix_max, 0.0)
                / static_cast<real_prec>(ix_max - ix_min);
  }
  return result;
}

#endif //BARCODE_TIME_STEP_HPP
