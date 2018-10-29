/*
 * Barcode
 * Copyright E.G.P. Bos
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <cmath> // fmod

#include "pacman.hpp"

////////////////////////////////////////////////////////////////
// Pacman functions: periodic boundary conditions handling
////////////////////////////////////////////////////////////////

// EGP: convenience function for applying periodic boundary conditions on a
// single coordinate.
// Puts coordinate in [0,L) range; have to exclude L, because it is equal to 0!
void pacman_coordinate(real_prec *x, real_prec L) {
  if (*x < 0.) {
    *x = std::fmod(*x, L);
    *x += L;
  }
  if (*x >= L) {
    *x = std::fmod(*x, L);
  }
}

// Convenience function for applying periodic boundary conditions on
// coordinate *distances*. This differs from absolute coordinates, because
// distances can never be longer than half a boxsize; if it is longer, then
// the distance to reach it the other way around is shorter, and we want
// the shortest distance.
// N.B.: only for single coordinate components, not for full distances
// (so not for dr = sqrt(dx**2 + dy**2 + dz**2))!
// N.B. 2: assumes the differences are already in pacman coordinates, so
// not larger than the boxsize L. Use pacman_coordinate function first
// otherwise.
// So put d_x in [-L/2,L/2] range (note: unlike pacman_coordinate, this is a
// range with two inclusive boundaries).
void pacman_difference(real_prec *d_x, real_prec L) {
  if (*d_x > L/2)
    *d_x = L - *d_x;
  if (*d_x < -(L/2))
    *d_x = L + *d_x;
}

real_prec pacman_d_x_from_d_x(real_prec d_x, real_prec L) {
  if (d_x > L/2)
    d_x = L - d_x;
  if (d_x < -(L/2))
    d_x = L + d_x;
  return d_x;
}

// Put d_ix in [-N/2,N/2] range (again, two inclusive boundaries)
int pacman_d_ix_from_d_ix(int d_ix, int N) {
  if (d_ix > N/2)
    d_ix = N - d_ix;
  if (d_ix < -(N/2))
    d_ix = N + d_ix;
  return d_ix;
}
