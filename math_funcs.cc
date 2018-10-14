/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include "struct_main.h"

#include <math.h>
#include <stdlib.h>
#include <stdexcept>  // std::logic_error

#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>

#include <string>

#include "fftw_array.h"
#include "planck/paramfile.h"
#include "IOfunctionsGen.h"

#include "BarcodeException.h"

#include "convenience.h"

// EGP: for nth_element, max_element, min_element:
#include <algorithm>

// EGP: resolution independent random grid
#include "random.hpp"

/*
   F.S. Kitaura 2007-2013
   */

using namespace std;

////////////////////////////////////////////////////////////////
// General scalar math functions
////////////////////////////////////////////////////////////////

// EGP: added odd to check for odd-ness (instead of even-ness).
// Taken from http://forums.devshed.com/software-design-43/quick-algorithm-to-determine-odd-even-numbers-29843.html
int odd(int inputval)
{
  return inputval & 1;
}

real_prec factorial(int term)
{
  real_prec out=1.0;

  if (term>0)
    for (int i=1;i<=term;i++)
      out*=static_cast<real_prec>(i*1.0);
  else
    out=static_cast<real_prec>(1.);

  return out;
}

real_prec power_mean(real_prec x, real_prec y, real_prec p)
{
  real_prec mean;
  if (p == 0)
  {
    mean = sqrt(x*y);
  }
  else
  {
    mean = pow( ( pow(x, p) + pow(y, p) )/2, 1/p);
  }
  return mean;
}




////////////////////////////////////////////////////////////////
// Array functions / operators
// TODO: when the new array class (vector derived) is done, make
//       versions of these functions that take such a vector as
//       argument (by reference).
////////////////////////////////////////////////////////////////

real_prec min_arr ( ULONG factor, real_prec *in )
{
  real_prec *firstn=in;
  real_prec *lastn=in+factor;
  real_prec minn;
  //get min
  real_prec *min=min_element(firstn,lastn);
  minn = *min;
  return minn;
}

real_prec max_arr ( ULONG factor, real_prec *in )
{
  real_prec *firstn=in;
  real_prec *lastn=in+factor;
  real_prec max;
  //get max
  real_prec *maxn = max_element(firstn,lastn);
  max =*maxn;
  return max;
}

// EGP: added mean_arr
real_prec mean_arr ( ULONG size, real_prec *in )
{
  real_prec mean = 0.0;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:mean)
#endif // MULTITHREAD
  for(ULONG i=0; i<size; i++)
    mean+=in[i];
  return mean/static_cast<real_prec>(size);
}

// EGP: added median_arr
real_prec median_arr (ULONG size, real_prec *in)
{
  //real_prec * copy = (real_prec *) malloc(size * sizeof(real_prec));
  //real_prec * copy = reinterpret_cast<real_prec *>( malloc(size * sizeof(real_prec))) ); // EGP: fixes "old style cast" warning
  real_prec *copy = new real_prec[size]; // EGP: suggested by Johan; C++ way of doing it
  //fftw_array<real_prec> copy(size);
  copyArray(in, copy, size);

  real_prec median;

  if (odd(static_cast<int>(size)))
  {
    nth_element(copy, copy + size/2, copy + size);
    median = copy[size/2];
  }
  else
  {
    nth_element(copy, copy + size/2 - 1, copy + size);
    median = copy[size/2-1]/2.;
    median += *min_element(copy + size/2, copy + size)/2.;
  }

  //free( copy );
  delete [] copy; // EGP: this should be used when using new, not free!!!
  return median;
}

// EGP: added std_arr
real_prec std_arr (ULONG size, real_prec *in)
{
  real_prec s_sq = 0.;
  real_prec mean = mean_arr(size, in);
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:s_sq)
#endif // MULTITHREAD
  for (ULONG i = 0; i < size; i++)
    s_sq += pow(in[i] - mean, 2);
  s_sq /= static_cast<real_prec>(size-1); // "corrected sample standard deviation" (squared)

  return sqrt(s_sq);
}

void complexfactor_mult (ULONG factor, real_prec in_a, complex_prec *in_b, complex_prec *out )
{
  for (ULONG i = 0; i < factor; i++)
  {
    re(out[i]) = in_a*re(in_b[i]);
    im(out[i]) = in_a*im(in_b[i]);
  }
}




////////////////////////////////////////////////////////////////
// Pacman functions
////////////////////////////////////////////////////////////////

// EGP: convenience function for applying periodic boundary conditions on a
// single coordinate.
// Puts coordinate in [0,L) range; have to exclude L, because it is equal to 0!
void pacman_coordinate(real_prec *x, real_prec L) {
  if (*x < 0.) {
    *x = fmod(*x, L);
    *x += L;
  }
  if (*x >= L) {
    *x = fmod(*x, L);
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




////////////////////////////////////////////////////////////////
// CIC: helper stuff and interpolation
////////////////////////////////////////////////////////////////


// EGP: convenience function for getting the 8 nearest cells to a coordinate.
// Since these 8 form a cube, only two cells need to be returned, those that
// form a diagonal of the cube. cell_index1 is the one whose cell center
// coordinates xc_1 are lower than the coordinate x's, i.e. xc_1 < x, etc., and
// cell_index2 has the higher cell center coordinates xc_2, xc_2 > x.
// N.B.: cell_index arrays must have length 3!
void getCICcells(ULONG N1, ULONG N2, ULONG N3, real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec x, real_prec y, real_prec z, ULONG *cell_index1, ULONG *cell_index2)
{
  real_prec xpos, ypos, zpos;
  
  xpos = x - num_0_5*d1;
  ypos = y - num_0_5*d2;
  zpos = z - num_0_5*d3;
  // periodic boundary conditions:
  pacman_coordinate(&xpos, L1);
  pacman_coordinate(&ypos, L2);
  pacman_coordinate(&zpos, L3);

  cell_index1[0] = static_cast<ULONG>(xpos/d1); // indices of the cell of the particle
  cell_index1[1] = static_cast<ULONG>(ypos/d2);
  cell_index1[2] = static_cast<ULONG>(zpos/d3);

  cell_index1[0] = (cell_index1[0] + N1) % N1; // periodic boundary conditions
  cell_index1[1] = (cell_index1[1] + N2) % N2;
  cell_index1[2] = (cell_index1[2] + N3) % N3;

  cell_index2[0] = (cell_index1[0] + 1) % N1;
  cell_index2[1] = (cell_index1[1] + 1) % N2;
  cell_index2[2] = (cell_index1[2] + 1) % N3;
}

// EGP: CIC weights, used in the CIC density estimator and CIC interpolator.
// tx is just 1-dx. Note that these arrays must have length 3, and cell_index1
// as well. The latter can be calculated using getCICcells.
void getCICweights(real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec x, real_prec y, real_prec z, ULONG *cell_index1, real_prec *dx, real_prec *tx)
{
  real_prec xpos, ypos, zpos, xc, yc, zc;
  
  xpos = x - num_0_5*d1;
  ypos = y - num_0_5*d2;
  zpos = z - num_0_5*d3;
// periodic boundary conditions:
  pacman_coordinate(&xpos, L1);
  pacman_coordinate(&ypos, L2);
  pacman_coordinate(&zpos, L3);

  xc = static_cast<real_prec>(cell_index1[0]); 
  yc = static_cast<real_prec>(cell_index1[1]);
  zc = static_cast<real_prec>(cell_index1[2]);

  dx[0] = xpos/d1 - xc;
  dx[1] = ypos/d2 - yc;
  dx[2] = zpos/d3 - zc;

  tx[0] = num_1 - dx[0];
  tx[1] = num_1 - dx[1];
  tx[2] = num_1 - dx[2];
}

// Single coordinate interpolator without all the array mess.
real_prec interpolate_CIC(ULONG N1, ULONG N2, ULONG N3, real_prec L1,
                          real_prec L2, real_prec L3, real_prec d1,
                          real_prec d2, real_prec d3, real_prec xp,
                          real_prec yp, real_prec zp, real_prec *input) {
  real_prec dx[3], tx[3];
  ULONG i[3], ii[3];
  getCICcells(N1, N2, N3, L1, L2, L3, d1, d2, d3, xp, yp, zp, i, ii);
  getCICweights(L1, L2, L3, d1, d2, d3, xp, yp, zp, i, dx, tx);

  #define FIELD(i, j, k) input[k[2]+N3*(j[1]+N2*i[0])]
  real_prec output = FIELD(i, i, i)    * tx[0]*tx[1]*tx[2] +
                     FIELD(ii, i, i)   * dx[0]*tx[1]*tx[2] +
                     FIELD(i, ii, i)   * tx[0]*dx[1]*tx[2] +
                     FIELD(i, i, ii)   * tx[0]*tx[1]*dx[2] +
                     FIELD(ii, ii, i)  * dx[0]*dx[1]*tx[2] +
                     FIELD(ii, i, ii)  * dx[0]*tx[1]*dx[2] +
                     FIELD(i, ii, ii)  * tx[0]*dx[1]*dx[2] +
                     FIELD(ii, ii, ii) * dx[0]*dx[1]*dx[2];
  #undef FIELD

  return output;
}


// EGP: this CIC interpolator was added because getDensity_CIC was used to do
// this, which obviously doesn't work.
void interpolate_CIC(ULONG N1, ULONG N2, ULONG N3, real_prec L1, real_prec L2,
                     real_prec L3, real_prec d1, real_prec d2, real_prec d3,
                     real_prec *xp, real_prec *yp, real_prec *zp,
                     real_prec *field_grid, ULONG N_OBJ,
                     real_prec *interpolation) {
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif  // MULTITHREAD
  for (ULONG n = 0; n < N_OBJ; n++) {
    interpolation[n] = interpolate_CIC(N1, N2, N3, L1, L2, L3, d1, d2, d3,
                                       xp[n], yp[n], zp[n], field_grid);
  }
}


////////////////////////////////////////////////////////////////
// TSC: interpolation
////////////////////////////////////////////////////////////////
// Interpolating with TSC is in a sense even easier than CIC,
// since the 27 nearest cells are trivially found; they are the
// cells around the cell that contains the particle. 
// Kernel values taken from p. 8 of:
// http://www-hpcc.astro.washington.edu/simulations/DARK_MATTER/adap.ps
////////////////////////////////////////////////////////////////

// Single coordinate interpolator without all the array mess.
real_prec interpolate_TSC(ULONG N1_ul, ULONG N2_ul, ULONG N3_ul, real_prec d1, real_prec d2, real_prec d3, real_prec xp,
                          real_prec yp, real_prec zp, real_prec *field) {
  unsigned N1 = static_cast<unsigned>(N1_ul);
  unsigned N2 = static_cast<unsigned>(N2_ul);
  unsigned N3 = static_cast<unsigned>(N3_ul);

  real_prec output = 0;
  // particle coordinates in "kernel units" (divided by d1/2/3)
  real_prec xp_kern = xp/d1;
  real_prec yp_kern = yp/d2;
  real_prec zp_kern = zp/d3;
  // indices of the cell of the particle
  unsigned ix_x_cell = static_cast<unsigned>(xp_kern);
  unsigned ix_y_cell = static_cast<unsigned>(yp_kern);
  unsigned ix_z_cell = static_cast<unsigned>(zp_kern);
  // cell_center coordinates (kernel units)
  real_prec x_cell = static_cast<real_prec>(ix_x_cell) + 0.5;
  real_prec y_cell = static_cast<real_prec>(ix_y_cell) + 0.5;
  real_prec z_cell = static_cast<real_prec>(ix_z_cell) + 0.5;
  // distance of particle from cell center
  real_prec dx_p_cell = xp_kern - x_cell;
  real_prec dy_p_cell = yp_kern - y_cell;
  real_prec dz_p_cell = zp_kern - z_cell;

  // weights per dimension
  real_prec wx[3], wy[3], wz[3];
  wx[1] = 0.75 - dx_p_cell*dx_p_cell;
  wy[1] = 0.75 - dy_p_cell*dy_p_cell;
  wz[1] = 0.75 - dz_p_cell*dz_p_cell;
  wx[0] = 0.5 * gsl_pow_2(1.5 - abs(dx_p_cell+1));
  wy[0] = 0.5 * gsl_pow_2(1.5 - abs(dy_p_cell+1));
  wz[0] = 0.5 * gsl_pow_2(1.5 - abs(dz_p_cell+1));
  wx[2] = 0.5 * gsl_pow_2(1.5 - abs(dz_p_cell-1));
  wy[2] = 0.5 * gsl_pow_2(1.5 - abs(dz_p_cell-1));
  wz[2] = 0.5 * gsl_pow_2(1.5 - abs(dz_p_cell-1));

  // indices per dimension, periodic box-ized
  unsigned ix_x[3], ix_y[3], ix_z[3];
  ix_x[1] = ix_x_cell;
  ix_y[1] = ix_y_cell;
  ix_z[1] = ix_z_cell;
  ix_x[0] = (ix_x_cell - 1 + N1) % N1;
  ix_y[0] = (ix_y_cell - 1 + N2) % N2;
  ix_z[0] = (ix_z_cell - 1 + N3) % N3;
  ix_x[2] = (ix_x_cell + 1) % N1;
  ix_y[2] = (ix_y_cell + 1) % N2;
  ix_z[2] = (ix_z_cell + 1) % N3;

  // add the 27 contributions
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k) {
        unsigned ix_f = (ix_x[i]*N2 + ix_y[j])*N3 + ix_z[k];
        output  +=    wx[i]    *   wy[j]     *   wz[k] * field[ix_f] ;
      }

  return output;
}


void interpolate_TSC(ULONG N1, ULONG N2, ULONG N3, real_prec d1, real_prec d2, real_prec d3, real_prec *xp, real_prec *yp,
                     real_prec *zp, real_prec *field_grid, ULONG N_OBJ, real_prec *interpolation) {
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif  // MULTITHREAD
  for (ULONG n = 0; n < N_OBJ; n++) {
    interpolation[n] = interpolate_TSC(N1, N2, N3, d1, d2, d3, xp[n], yp[n], zp[n], field_grid);
  }
}


// multi-field variant
// Single coordinate interpolator without all the array mess.
void interpolate_TSC_multi(ULONG N1_ul, ULONG N2_ul, ULONG N3_ul, real_prec d1, real_prec d2, real_prec d3, real_prec xp,
                           real_prec yp, real_prec zp, real_prec **fields, int N_fields, real_prec *out) {
  unsigned N1 = static_cast<unsigned>(N1_ul);
  unsigned N2 = static_cast<unsigned>(N2_ul);
  unsigned N3 = static_cast<unsigned>(N3_ul);

  for (int ix_fields = 0; ix_fields < N_fields; ++ix_fields) {
    out[ix_fields] = 0;
  }
  // particle coordinates in "kernel units" (divided by d1/2/3)
  real_prec xp_kern = xp/d1;
  real_prec yp_kern = yp/d2;
  real_prec zp_kern = zp/d3;
  // indices of the cell of the particle
  unsigned ix_x_cell = static_cast<unsigned>(xp_kern);
  unsigned ix_y_cell = static_cast<unsigned>(yp_kern);
  unsigned ix_z_cell = static_cast<unsigned>(zp_kern);
  // cell_center coordinates (kernel units)
  real_prec x_cell = static_cast<real_prec>(ix_x_cell) + 0.5;
  real_prec y_cell = static_cast<real_prec>(ix_y_cell) + 0.5;
  real_prec z_cell = static_cast<real_prec>(ix_z_cell) + 0.5;
  // distance of particle from cell center
  real_prec dx_p_cell = xp_kern - x_cell;
  real_prec dy_p_cell = yp_kern - y_cell;
  real_prec dz_p_cell = zp_kern - z_cell;

  // weights per dimension
  real_prec wx[3], wy[3], wz[3];
  wx[1] = 0.75 - dx_p_cell*dx_p_cell;
  wy[1] = 0.75 - dy_p_cell*dy_p_cell;
  wz[1] = 0.75 - dz_p_cell*dz_p_cell;
  wx[0] = 0.5 * gsl_pow_2(1.5 - abs(dx_p_cell+1));
  wy[0] = 0.5 * gsl_pow_2(1.5 - abs(dy_p_cell+1));
  wz[0] = 0.5 * gsl_pow_2(1.5 - abs(dz_p_cell+1));
  wx[2] = 0.5 * gsl_pow_2(1.5 - abs(dz_p_cell-1));
  wy[2] = 0.5 * gsl_pow_2(1.5 - abs(dz_p_cell-1));
  wz[2] = 0.5 * gsl_pow_2(1.5 - abs(dz_p_cell-1));

  // indices per dimension, periodic box-ized
  unsigned ix_x[3], ix_y[3], ix_z[3];
  ix_x[1] = ix_x_cell;
  ix_y[1] = ix_y_cell;
  ix_z[1] = ix_z_cell;
  ix_x[0] = (ix_x_cell - 1 + N1) % N1;
  ix_y[0] = (ix_y_cell - 1 + N2) % N2;
  ix_z[0] = (ix_z_cell - 1 + N3) % N3;
  ix_x[2] = (ix_x_cell + 1) % N1;
  ix_y[2] = (ix_y_cell + 1) % N2;
  ix_z[2] = (ix_z_cell + 1) % N3;

  // add the 27 contributions
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k) {
        unsigned ix_f =         (ix_x[i]*N2 + ix_y[j])*N3 + ix_z[k];
        // ...for each field
        for (int ix_fields = 0; ix_fields < N_fields; ++ix_fields) {
          real_prec *field = fields[ix_fields];
          out[ix_fields]  +=    wx[i]    *   wy[j]     *   wz[k] * field[ix_f];
        }
      }
}


void interpolate_TSC_multi(ULONG N1, ULONG N2, ULONG N3, real_prec d1, real_prec d2, real_prec d3, real_prec *xp,
                           real_prec *yp, real_prec *zp, ULONG N_OBJ, real_prec **field_grids,
                           real_prec **interpolations, int N_fields) {
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif  // MULTITHREAD
  for (ULONG n = 0; n < N_OBJ; n++) {
    real_prec out[N_fields];
    // interpolations[n] = interpolate_TSC(N1, N2, N3, L1, L2, L3, d1, d2, d3,
    //                                    xp[n], yp[n], zp[n], field_grid);
    interpolate_TSC_multi(N1, N2, N3, d1, d2, d3, xp[n], yp[n], zp[n], field_grids, N_fields, out);
    for (int ix_fields = 0; ix_fields < N_fields; ++ix_fields) {
      interpolations[ix_fields][n] = out[ix_fields];
    }
  }
}






////////////////////////////////////////////////////////////////
// k-vector calculators
////////////////////////////////////////////////////////////////

real_prec k_squared(unsigned int i, unsigned int j, unsigned int k, real_prec L1, real_prec L2, real_prec L3,
                    unsigned int N1, unsigned int N2, unsigned int N3)
{
  real_prec k2=0.;
  real_prec  kfac=static_cast<real_prec>(2.*M_PI);

  real_prec kx=0.;
  real_prec ky=0.;
  real_prec kz=0.;

  if (i<=N1/2) kx = kfac/L1*static_cast<real_prec>(i);
  else kx= -kfac/L1*static_cast<real_prec>(N1-i);

  if (j<=N2/2) ky= kfac/L2*static_cast<real_prec>(j);
  else ky= -kfac/L2*static_cast<real_prec>(N2-j);

  if (k<=N3/2) kz = kfac/L3*static_cast<real_prec>(k);
  else kz = -kfac/L3*static_cast<real_prec>(N3-k);

  k2=kx*kx+ky*ky+kz*kz;

  return(k2);
}

// EGP: Note: calc_kx, y and z were exactly the same functions! Replaced with calc_ki
real_prec calc_ki(unsigned int i, real_prec Li, unsigned int Ni) {
  real_prec kfac=static_cast<real_prec>(2.*M_PI/Li);
  real_prec ki=0.;

  if (i<=Ni/2)
    ki = kfac*static_cast<real_prec>(i);
  else
    ki = -kfac*static_cast<real_prec>(Ni-i);

  return(ki);
}

real_prec calc_kx(unsigned int i, real_prec L1, unsigned int N1) {
  return(calc_ki(i, L1, N1));
}

real_prec calc_ky(unsigned int j, real_prec L2, unsigned int N2) {
  return(calc_ki(j, L2, N2));
}

real_prec calc_kz(unsigned int k, real_prec L3, unsigned int N3) {
  return(calc_ki(k, L3, N3));
}





////////////////////////////////////////////////////////////////
// Convolution and kernel functions
////////////////////////////////////////////////////////////////

void convolve(real_prec L1, real_prec L2, real_prec L3, unsigned N1, unsigned N2, unsigned N3, real_prec *in, real_prec *out,
              real_prec smol, bool zeropad, int filtertype)
{
  bool gauss=false;
  bool tophat=false;
  bool errfunc=false;

  switch (filtertype)
  {
    case 1:
      gauss=true;
      break;
    case 2:
      tophat=true;
      break;
    case 3:
      errfunc=true;
      break;
  }

  //EGP  ULONG N=N1*N2*N3;

  unsigned Nzp1=N1;
  unsigned Nzp2=N2;
  unsigned Nzp3=N3;

  real_prec Lzp1=L1;
  real_prec Lzp2=L2;
  real_prec Lzp3=L3;

  if (zeropad==true)
  {
    Lzp1=L1*num_2;
    Lzp2=L2*num_2;
    Lzp3=L3*num_2;

    Nzp1=N1*2;
    Nzp2=N2*2;
    Nzp3=N3*2;
  }

  ULONG Nzp=Nzp1*Nzp2*Nzp3;

  fftw_array<complex_prec>  AUX(Nzp), AUXb(Nzp);
  fftw_array<real_prec> wkernel(Nzp),dummy(Nzp);

  if (zeropad==true)
  {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
    for(ULONG i=0;i<Nzp;i++)
      dummy[i]=0.0;

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
    for(unsigned i=0;i<N1;i++)
      for(unsigned j=0;j<N2;j++)
        for(unsigned k=0;k<N3;k++)
        {
          unsigned k_zp=k+N3/2;
          unsigned j_zp=Nzp3*(j+N2/2);
          unsigned i_zp=Nzp3*Nzp2*(i+N1/2);

          unsigned k_eff=k;
          unsigned j_eff=N3*j;
          unsigned i_eff=N3*N2*i;

          dummy[i_zp+j_zp+k_zp]=in[i_eff+j_eff+k_eff];
        }
  }
  else
  {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
    for(ULONG i=0;i<Nzp;i++)
      dummy[i]=in[i];
  }

  FFT3dR2C (Nzp1,Nzp2,Nzp3,dummy,AUXb);

  //EGP  real_prec asmth=1.0;
  real_prec u;

#ifdef MULTITHREAD
#pragma omp parallel for private(u)
#endif // MULTITHREAD
  for (unsigned i=0;i<Nzp1;i++)
    for (unsigned j=0;j<Nzp2;j++)
      for (unsigned k=0;k<Nzp3;k++)
      {
        real_prec k2=k_squared(i,j,k,Lzp1,Lzp2,Lzp3,Nzp1,Nzp2,Nzp3);

        real_prec rS=smol;
        real_prec rS2=rS*rS;
        real_prec kcut=smol;//2.*M_PI/rS;
        real_prec sigma=static_cast<real_prec>(.3);

        if (tophat==true)
        {
          u = sqrt(k2);

          if (u>kcut)
            wkernel[k+(Nzp3)*(j+Nzp2*i)]=0.0;
          else
            wkernel[k+(Nzp3)*(j+Nzp2*i)]=1.0;
        }

        if (errfunc==true)
        {
          u = static_cast<real_prec>((sqrt(k2)-kcut)/(sqrt(2.)*sigma));
          real_prec fac = static_cast<real_prec>(erfc(u));
          wkernel[k+(Nzp3)*(j+Nzp2*i)]=fac;
        }

        if (gauss==true)
          wkernel[k+(Nzp3)*(j+Nzp2*i)]=static_cast<real_prec>(exp(-k2*rS2/2.));
      }


#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for(ULONG i=0;i<Nzp;i++)
  {
    re(AUX[i])=wkernel[i];
    im(AUX[i])=0.0;
  }

  {
    FFT3d (Nzp1,Nzp2,Nzp3, to_Rspace, AUX, AUX);

    real_prec wtot=0.;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:wtot)
#endif // MULTITHREAD
    for(ULONG i=0;i<Nzp;i++)
      wtot+=re(AUX[i]);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
    for(ULONG i=0;i<Nzp;i++)
    {
      wkernel[i]/=wtot;
    }
  }


#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for(ULONG i=0;i<Nzp;i++)
  {
    re(AUX[i])=wkernel[i];
    im(AUX[i])=0.0;
  }

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for(ULONG i=0;i<Nzp;i++)
  {
    re(AUXb[i])= re(AUXb[i])*re(AUX[i]);
    im(AUXb[i])= im(AUXb[i])*re(AUX[i]);
  }

  FFT3d (Nzp1,Nzp2,Nzp3, to_Rspace, AUXb, AUXb);
  if (zeropad==true)
  {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
    for(unsigned i=0;i<N1;i++)
      for(unsigned j=0;j<N2;j++)
        for(unsigned k=0;k<N3;k++)
        {
          unsigned k_zp=k+N3/2;
          unsigned j_zp=Nzp3*(j+N2/2);
          unsigned i_zp=Nzp3*Nzp2*(i+N1/2);

          unsigned k_eff=k;
          unsigned j_eff=N3*j;
          unsigned i_eff=N3*N2*i;

          out[i_eff+j_eff+k_eff]=re(AUXb[i_zp+j_zp+k_zp]);
        }
  }
  else
  {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
    for(ULONG i=0;i<Nzp;i++)
      out[i]=re(AUXb[i]);
  }
}


void kernelcomp(real_prec L1, real_prec L2, real_prec L3, unsigned N1, unsigned N2, unsigned N3, real_prec smol, int filtertype,
                struct DATA *data)
{
  bool gauss=false;
  bool errfunc=false;
  bool tophat=false;

  switch (filtertype)
  {
    case 1:
      gauss=true;
      break;
    case 2:
      tophat=true;
      break;
    case 3:
      errfunc=true;
      break;
  }

  ULONG N=N1*N2*N3;

  fftw_array<real_prec> out(N);
  fftw_array<complex_prec> AUX(N);

  //EGP  real_prec asmth=1.0;
  real_prec u;

#ifdef MULTITHREAD
#pragma omp parallel for private(u)
#endif // MULTITHREAD
  for (unsigned i=0;i<N1;i++)
    for (unsigned j=0;j<N2;j++)
      for (unsigned k=0;k<N3;k++)
      {
        ULONG ii=k+N3*(j+N2*i);

        real_prec k2=k_squared(i,j,k,L1,L2,L3,N1,N2,N3);

        real_prec rS=smol;
        real_prec rS2=rS*rS;
        real_prec kcut=smol;//2.*M_PI/rS;
        real_prec sigma=static_cast<real_prec>(.3);

        if (tophat==true)
        {
          u = sqrt(k2);

          if (u>kcut)
            re(AUX[ii])=0.0;
          else
            re(AUX[ii])=1.0;
        }

        if (errfunc==true)
        {
          u = static_cast<real_prec>((sqrt(k2)-kcut)/(sqrt(2.)*sigma));
          real_prec fac = static_cast<real_prec>(erfc(u));
          re(AUX[ii])=fac;
        }

        if (gauss==true)
          re(AUX[ii])=static_cast<real_prec>(exp(-k2*rS2/2.));


        im(AUX[ii])=0.0;
      }


#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for(ULONG i=0;i<N;i++)
  {
    out[i]=re(AUX[i]);
    im(AUX[i])=0.0;
  }

  FFT3d (N1,N2,N3, to_Rspace, AUX, AUX);

  real_prec wtot=0.;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:wtot)
#endif // MULTITHREAD
  for(ULONG i=0;i<N;i++)
    wtot+=re(AUX[i]);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for(ULONG i=0;i<N;i++)
  {
    out[i]/=wtot;
  }

  int bmax=100;
  char buffsl[bmax];
  sprintf(buffsl,"r%d",int(smol));
  string fname= data->numerical->dir + string("auxkernel")+buffsl;
  dump_scalar(out, N1, N2, N3, fname);
}


void convcomp(unsigned N1, unsigned N2, unsigned N3, real_prec *in, real_prec *out, real_prec smol, string dir)
{
  //bool gauss=false;
  //bool errfunc=false;
  //bool tophat=false;

  //switch (filtertype)
  //{
    //case 1:
      //gauss=true;
      //break;
    //case 2:
      //tophat=true;
      //break;
    //case 3:
      //errfunc=true;
      //break;
  //}


  ULONG N=N1*N2*N3;

  fftw_array<complex_prec> AUX(N);

  FFT3dR2C (N1,N2,N3,in,AUX);

  {
    int bmax=100;
    char buffsl[bmax];
    sprintf(buffsl,"r%d",int(smol));
    //EGP    char * fileN;
    string fname= dir + string("auxkernel")+buffsl;
    get_scalar(fname,out,N1,N2,N3);
  }

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (unsigned i=0;i<N1;i++)
    for (unsigned j=0;j<N2;j++)
      for (unsigned k=0;k<N3;k++)
      {
        ULONG kk=k+(N3)*(j+N2*i);
        ULONG ii=k+(N3)*(j+N2*i);

        re(AUX[kk])= re(AUX[kk])*out[ii];
        im(AUX[kk])= im(AUX[kk])*out[ii];
      }

  FFT3dC2R (N1,N2,N3,AUX,out);
}


void convcompb(unsigned N1, unsigned N2, unsigned N3, real_prec *in, real_prec *out)
{
  ULONG N=N1*N2*N3;

  fftw_array<complex_prec> AUX(N);

  FFT3dR2C (N1,N2,N3,in,AUX);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (unsigned i=0;i<N1;i++)
    for (unsigned j=0;j<N2;j++)
      for (unsigned k=0;k<N3;k++)
      {
        ULONG kk=k+(N3)*(j+N2*i);
        ULONG ii=k+(N3)*(j+N2*i);

        re(AUX[kk])= re(AUX[kk])*out[ii];
        im(AUX[kk])= im(AUX[kk])*out[ii];
      }

  FFT3dC2R (N1,N2,N3,AUX,out);
}






////////////////////////////////////////////////////////////////
// Gradient and other 3D derivative functions
////////////////////////////////////////////////////////////////

void gradfft(unsigned N1, unsigned N2, unsigned N3, real_prec L1, real_prec L2, real_prec L3, real_prec *in,
             real_prec *out, unsigned int dim)
{
  // ULONG N=N1*N2*N3;
  unsigned N3half = N3/2 + 1;
  ULONG Nhalf = N1 * N2 * N3half;
  // fftw_array<complex_prec> AUX(N),AUX2(N);
  fftw_array<complex_prec> AUX(Nhalf);

  // complexify_array(in, AUX2, N);
  // FFT3d(N1,N2,N3, true, AUX2, AUX);
  fftR2C(N1, N2, N3, in, AUX);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (unsigned i=0;i<N1;i++)
    for (unsigned j=0;j<N2;j++)
      // for (unsigned k=0;k<N3;k++)
      for (unsigned k = 0; k < N3half; ++k) {
        real_prec kl;
        switch (dim)
        {
          case 1:
            kl=calc_kx(i,L1,N1);
            break;
          case 2:
            kl=calc_ky(j,L2,N2);
            break;
          case 3:
            kl=calc_kz(k,L3,N3);
            break;
          default:
            throw BarcodeException("In gradfft: dim should be either 1, 2 or 3!");
            break;
        }

        ULONG ll = k + N3half * (j + N2*i);

        // EGP: switched minus sign between terms! Was wrong way round.
        // EGP: ... and switched it back; Fourier def. assumption was wrong!
        real_prec dummy = re(AUX[ll]);  // for in-place
        re(AUX[ll]) = -kl * im(AUX[ll]);
        im(AUX[ll]) =  kl * dummy;
        
        // From http://math.mit.edu/~stevenj/fft-deriv.pdf: Nyquist components
        // should be zero, in every dimension, otherwise the result is complex.
        if ((i==N1/2) || (j==N2/2) || (k==N3/2))
        {
          re(AUX[ll]) = 0.;
          im(AUX[ll]) = 0.;
        }
      }
  // FFT3d(N1,N2,N3, false, AUX2, AUX);
  // real_part_array(AUX, out, N);
  fftC2R(N1, N2, N3, AUX, out);
}


void gradfindif(unsigned N1, real_prec L1, real_prec *in, real_prec *out, unsigned int dim)
{
  if (N1 > INT_MAX) {
    std::string msg = "Box dimensions must not be larger than INT_MAX in gradfindif!";
    throw std::logic_error(msg);
    // otherwise the signed integer indices below will overflow
  }
  real_prec fac=static_cast<real_prec>(N1/(2.*L1));

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for(unsigned x = 0; x < N1; x++)
    for(unsigned y = 0; y < N1; y++)
      for(unsigned z = 0; z < N1; z++)
      {
        unsigned xr, xrr;
        int xl, xll;
        unsigned yr, yrr;
        int yl, yll;
        unsigned zr, zrr;
        int zl, zll;

        xrr = xr = x;
        yrr = yr = y;
        zrr = zr = z;

        xll = xl = static_cast<int>(x);
        yll = yl = static_cast<int>(y);
        zll = zl = static_cast<int>(z);

        int *il, *ill;
        unsigned *ir, *irr, *ii;

        switch (dim) {
          case 1: {
            ii = &x; il = &xl; ill = &xll; ir = &xr; irr = &xrr;
            break;
          }
          case 2: {
            ii = &y; il = &yl; ill = &yll; ir = &yr; irr = &yrr;
            break;
          }
          case 3: {
            ii = &z; il = &zl; ill = &zll; ir = &zr; irr = &zrr;
            break;
          }
          default: {
            throw BarcodeException("dim must be 1, 2 or 3 in gradfindif!");
          }
        }

        *ir = *ii + 1;
        *il = static_cast<int>(*ii) - 1;
        *irr = *ii + 2;
        *ill = static_cast<int>(*ii) - 2;
        if(*ir >= N1)
          *ir -= N1;
        if(*irr >= N1)
          *irr -= N1;
        if(*il < 0)
          *il += N1;
        if(*ill < 0)
          *ill += N1;

        ULONG ix = z+N1*(y+N1*x);
        ULONG l = static_cast<unsigned>(zl)+N1*(static_cast<unsigned>(yl)+N1*static_cast<unsigned>(xl));
        ULONG r = zr+N1*(yr+N1*xr);
        ULONG ll = static_cast<unsigned>(zll)+N1*(static_cast<unsigned>(yll)+N1*static_cast<unsigned>(xll));
        ULONG rr = zrr+N1*(yrr+N1*xrr);

        out[ix] = -static_cast<real_prec>(fac*((4.0 / 3)*(in[l] - in[r]) - (1.0 / 6) * (in[ll] - in[rr])));
      }
}

// EGP: added this simultaneous (a) one component fourier space gradient and (b) inverse laplacian calculator, in Fourier space:
void grad_inv_lap_FS(unsigned N1, unsigned N2, unsigned N3, real_prec L1, real_prec L2, real_prec L3, complex_prec *in,
                     complex_prec *out, unsigned int index, bool rfft)
{
  unsigned kz_max = N3;
  if (rfft)
    kz_max = N3/2 + 1;

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif // MULTITHREAD
  for (unsigned i=0;i<N1;i++)
    for (unsigned j=0;j<N2;j++)
      for (unsigned k = 0; k < kz_max; ++k) {
        ULONG ii=k + kz_max * (j + N2*i);

        real_prec kx=calc_kx(i,L1,N1);
        real_prec ky=calc_ky(j,L2,N2);
        real_prec kz=calc_kz(k,L3,N3);

        real_prec kmod = kx*kx + ky*ky + kz*kz;
        real_prec fac_kmod = 0.;
        if (kmod > 0)
          fac_kmod = 1/kmod;

        real_prec ki_over_kmod=0.;
        switch (index)
        {
          case 1:
            ki_over_kmod=kx*fac_kmod;
            break;
          case 2:
            ki_over_kmod=ky*fac_kmod;
            break;
          case 3:
            ki_over_kmod=kz*fac_kmod;
            break;
          default:
            throw BarcodeException("In grad_inv_lap_FS: index must be either 1, 2 or 3!");
            break;
        }
        real_prec dummy = re(in[ii]); // to make in-place possible (in == out)
        re(out[ii]) =  ki_over_kmod*im(in[ii]);
        im(out[ii]) = -ki_over_kmod*dummy;

        // From http://math.mit.edu/~stevenj/fft-deriv.pdf: Nyquist components
        // should be zero, in every dimension, otherwise the result is complex.
        // This function is not actually treated in that document, but I would
        // say it also counts as an odd-ordered derivative.
        if ((i==N1/2) || (j==N2/2) || (k==N3/2))
        {
          re(out[ii]) = 0.;
          im(out[ii]) = 0.;
        }
      }
}





////////////////////////////////////////////////////////////////
// Random functions
////////////////////////////////////////////////////////////////

real_prec GR_NUM(gsl_rng * SEED, real_prec sigma ,int GR_METHOD)
{
  real_prec val=0.;
  switch (GR_METHOD)
  {
    case 0:
      {
        val=static_cast<real_prec>(gsl_ran_gaussian(SEED, sigma));
      }
      break;

    case 1:
      {
        val=static_cast<real_prec>(gsl_ran_gaussian_ziggurat(SEED, sigma));
      }
      break;

    case 2:
      {
        val=static_cast<real_prec>(gsl_ran_gaussian_ratio_method(SEED, sigma));
      }
      break;

  }
  return(val);
}


void create_GARFIELD(unsigned N1, unsigned N2, unsigned N3, real_prec L1, real_prec L2, real_prec L3,real_prec *delta,real_prec * Power,gsl_rng * seed)
{
  /// FOLLOWING HUGO MARTEL
  //
  // New version, EGP: uses a "resolution independent random grid",
  // which makes sure that, when given the same seed, the random grid
  // that is generated for gridsize N is also generated for gridsize
  // 2*N, plus the higher k terms (smaller scale structure).
  // Even though there are no more random calls in the for loop, there
  // is the possibility that one array location is accessed simultaneously
  // by multiple threads, so no multithreading here! Unless you want to
  // put atomics around all the calls...
  

  ULONG N=N1*N2*N3;
  ULONG  Nhalf=N1*N2*(N3/2+1);

  real_prec Vol = L1*L2*L3;

  fftw_array<complex_prec> GRF_array(N),GRF_array2(Nhalf);

  // Build a resolution independent random grid
  vector< complex<real_prec> > random_grid = resolution_independent_random_grid_FS<real_prec>( N1, seed, false );
  complex_prec *random_grid_array = reinterpret_cast<complex_prec *>(random_grid.data());
  copyArray(random_grid_array, GRF_array, N);

// Factor for conversion from power spectrum to discrete fourier transform
// amplitude (see Martel 2005):
  real_prec ps2dft_amp = 1.;
#ifdef FOURIER_DEF_1
  ps2dft_amp = 1./Vol;
#endif // FOURIER_DEF_1
#ifdef FOURIER_DEF_2
  ps2dft_amp = static_cast<real_prec>(N*N)/Vol;
#endif // FOURIER_DEF_2

  // TODO:
  // We can shorten the code considerably by just multiplying the whole array
  // by sigma first in a separate (multithreaded) for-loop (regardless of all
  // the special cases in the complicated loop below) and only after that 
  // doing all the special cases below, but we can then leave out all the 
  // *= sigma parts. Actually, I think what is done here is just two things:
  // 1. making the array hermitian, and 2. a couple of special cases that need
  // to be non-imaginary. The former we can do in a separate function (which is
  // already in the code, I believe) and the latter we can just do with a few
  // direct assignments, instead of having to go through all the if's at every
  // step of the loop.
  
  //EGP  real_prec kr=0.;
  real_prec sigma=0.;
//#ifdef MULTITHREAD
//#pragma omp parallel for
//#endif // MULTITHREAD
  for (unsigned i=0 ; i<=N1/2;i++)
    for (unsigned j=0 ; j<=N2/2;j++)
      for (unsigned k=0 ; k<=N3/2;k++)
      {
        sigma = sqrt(ps2dft_amp * Power[k+N3*(j+N2*i)]/num_2);

        if( (i>0 && i<N1/2) && (j>0 && j<N2/2) && (k>0 && k<N3/2))
        {
          unsigned ii=N1-i; unsigned jj=N2-j; unsigned kk=N3-k;

          re(GRF_array[k+N3*(j+N2*i)]) *= sigma;
          im(GRF_array[k+N3*(j+N2*i)]) *= sigma;

          re(GRF_array[kk+N3*(jj+N2*ii)])=re(GRF_array[k+N3*(j+N2*i)]);
          im(GRF_array[kk+N3*(jj+N2*ii)])=-im(GRF_array[k+N3*(j+N2*i)]);

          ///**********************************************************

          re(GRF_array[k+N3*(j+N2*ii)]) *= sigma;
          im(GRF_array[k+N3*(j+N2*ii)]) *= sigma;

          re(GRF_array[kk+N3*(jj+N2*i)])=re(GRF_array[k+N3*(j+N2*ii)]);
          im(GRF_array[kk+N3*(jj+N2*i)])=-im(GRF_array[k+N3*(j+N2*ii)]);

          ///**********************************************************

          re(GRF_array[k+N3*(jj+N2*i)]) *= sigma;
          im(GRF_array[k+N3*(jj+N2*i)]) *= sigma;

          re(GRF_array[kk+N3*(j+N2*ii)])=re(GRF_array[k+N3*(jj+N2*i)]);
          im(GRF_array[kk+N3*(j+N2*ii)])=-im(GRF_array[k+N3*(jj+N2*i)]);

          ///**********************************************************

          re(GRF_array[kk+N3*(j+N2*i)]) *= sigma;
          im(GRF_array[kk+N3*(j+N2*i)]) *= sigma;

          re(GRF_array[k+N3*(jj+N2*ii)])=re(GRF_array[kk+N3*(j+N2*i)]);
          im(GRF_array[k+N3*(jj+N2*ii)])=-im(GRF_array[kk+N3*(j+N2*i)]);

        }

        if( (i>0 && i<N1/2) && (j>0 && j<N2/2) && (k==N3/2))
        {
          unsigned ii=N1-i; unsigned jj=N2-j;

          re(GRF_array[N3/2+N3*(j+N2*i)]) *= sigma;
          im(GRF_array[N3/2+N3*(j+N2*i)]) *= sigma;

          re(GRF_array[N3/2+N3*(jj+N2*ii)])=re(GRF_array[N3/2+N3*(j+N2*i)]);
          im(GRF_array[N3/2+N3*(jj+N2*ii)])=-im(GRF_array[N3/2+N3*(j+N2*i)]);

          ///**********************************************************

          re(GRF_array[N3/2+N3*(j+N2*ii)]) *= sigma;
          im(GRF_array[N3/2+N3*(j+N2*ii)]) *= sigma;

          re(GRF_array[N3/2+N3*(jj+N2*i)])=re(GRF_array[N3/2+N3*(j+N2*ii)]);
          im(GRF_array[N3/2+N3*(jj+N2*i)])=-im(GRF_array[N3/2+N3*(j+N2*ii)]);

        }

        if( (i>0 && i<N1/2) && (j==N2/2) && (k>0 && k<N3/2))
        {
          unsigned ii=N1-i; unsigned kk=N3-k;

          re(GRF_array[k+N3*(N2/2+N2*i)]) *= sigma;
          im(GRF_array[k+N3*(N2/2+N2*i)]) *= sigma;

          re(GRF_array[kk+N3*(N2/2+N2*ii)])=re(GRF_array[k+N3*(N2/2+N2*i)]);
          im(GRF_array[kk+N3*(N2/2+N2*ii)])=-im(GRF_array[k+N3*(N2/2+N2*i)]);

          ///**********************************************************

          re(GRF_array[k+N3*(N2/2+N2*ii)]) *= sigma;
          im(GRF_array[k+N3*(N2/2+N2*ii)]) *= sigma;

          re(GRF_array[kk+N3*(N2/2+N2*i)])=re(GRF_array[k+N3*(N2/2+N2*ii)]);
          im(GRF_array[kk+N3*(N2/2+N2*i)])=-im(GRF_array[k+N3*(N2/2+N2*ii)]);

        }

        if( (i==N1/2) && (j>0 && j<N2/2) && (k>0 && k<N3/2))
        {
          unsigned jj=N2-j; unsigned kk=N3-k;

          re(GRF_array[k+N3*(j+N2*N1/2)]) *= sigma;
          im(GRF_array[k+N3*(j+N2*N1/2)]) *= sigma;

          re(GRF_array[kk+N3*(jj+N2*N1/2)])=re(GRF_array[k+N3*(j+N2*N1/2)]);
          im(GRF_array[kk+N3*(jj+N2*N1/2)])=-im(GRF_array[k+N3*(j+N2*N1/2)]);

          ///**********************************************************

          re(GRF_array[k+N3*(jj+N2*N1/2)]) *= sigma;
          im(GRF_array[k+N3*(jj+N2*N1/2)]) *= sigma;

          re(GRF_array[kk+N3*(j+N2*N1/2)])=re(GRF_array[k+N3*(jj+N2*N1/2)]);
          im(GRF_array[kk+N3*(j+N2*N1/2)])=-im(GRF_array[k+N3*(jj+N2*N1/2)]);

        }

        if( (i==N1/2) && (j==N2/2) && (k>0 && k<N3/2))
        {
          unsigned kk=N3-k;

          re(GRF_array[k+N3*(N2/2+N2*N1/2)]) *= sigma;
          im(GRF_array[k+N3*(N2/2+N2*N1/2)]) *= sigma;

          re(GRF_array[kk+N3*(N2/2+N2*N1/2)])=re(GRF_array[k+N3*(N2/2+N2*N1/2)]);
          im(GRF_array[kk+N3*(N2/2+N2*N1/2)])=-im(GRF_array[k+N3*(N2/2+N2*N1/2)]);

          ///**********************************************************
        }

        if( (i>0 && i<N1/2) && (j==N2/2) && (k==N3/2))
        {
          unsigned ii=N1-i;

          re(GRF_array[N3/2+N3*(N2/2+N2*i)]) *= sigma;
          im(GRF_array[N3/2+N3*(N2/2+N2*i)]) *= sigma;

          re(GRF_array[N3/2+N3*(N2/2+N2*ii)])=re(GRF_array[N3/2+N3*(N2/2+N2*i)]);
          im(GRF_array[N3/2+N3*(N2/2+N2*ii)])=-im(GRF_array[N3/2+N3*(N2/2+N2*i)]);

          ///**********************************************************
        }

        if( (i==N1/2) && (j>0 && j<N2/2) && (k==N3/2))
        {
          unsigned jj=N2-j;

          re(GRF_array[N3/2+N3*(j+N2*N1/2)]) *= sigma;
          im(GRF_array[N3/2+N3*(j+N2*N1/2)]) *= sigma;

          re(GRF_array[N3/2+N3*(jj+N2*N1/2)])=re(GRF_array[N3/2+N3*(j+N2*N1/2)]);
          im(GRF_array[N3/2+N3*(jj+N2*N1/2)])=-im(GRF_array[N3/2+N3*(j+N2*N1/2)]);

          ///**********************************************************
        }

        if( (i==N1/2) && (j==N2/2) && (k==N3/2))
        {
          re(GRF_array[N3/2+N3*(N2/2+N2*N1/2)]) *= sqrt(num_2)*sigma;
          im(GRF_array[N3/2+N3*(N2/2+N2*N1/2)])=0.;

          ///**********************************************************
        }

        if( (i>0 && i<N1/2) && (j>0 && j<N2/2) && k==0)
        {
          unsigned ii=N1-i; unsigned jj=N2-j;

          re(GRF_array[0+N3*(j+N2*i)]) *= sigma;
          im(GRF_array[0+N3*(j+N2*i)]) *= sigma;

          re(GRF_array[0+N3*(jj+N2*ii)])=re(GRF_array[k+N3*(j+N2*i)]);
          im(GRF_array[0+N3*(jj+N2*ii)])=-im(GRF_array[k+N3*(j+N2*i)]);

          ///**********************************************************

          re(GRF_array[0+N3*(j+N2*ii)]) *= sigma;
          im(GRF_array[0+N3*(j+N2*ii)]) *= sigma;

          re(GRF_array[0+N3*(jj+N2*i)])=re(GRF_array[0+N3*(j+N2*ii)]);
          im(GRF_array[0+N3*(jj+N2*i)])=-im(GRF_array[0+N3*(j+N2*ii)]);
        }

        if( (i>0 && i<N1/2) && (j==0) && (k>0 && k<N3/2))
        {
          unsigned ii=N1-i; unsigned kk=N3-k;

          re(GRF_array[k+N3*(0+N2*i)]) *= sigma;
          im(GRF_array[k+N3*(0+N2*i)]) *= sigma;

          re(GRF_array[kk+N3*(0+N2*ii)])=re(GRF_array[k+N3*(0+N2*i)]);
          im(GRF_array[kk+N3*(0+N2*ii)])=-im(GRF_array[k+N3*(0+N2*i)]);

          ///**********************************************************

          re(GRF_array[k+N3*(0+N2*ii)]) *= sigma;
          im(GRF_array[k+N3*(0+N2*ii)]) *= sigma;

          re(GRF_array[kk+N3*(0+N2*i)])=re(GRF_array[k+N3*(0+N2*ii)]);
          im(GRF_array[kk+N3*(0+N2*i)])=-im(GRF_array[k+N3*(0+N2*ii)]);
        }

        if( (i==0) && (j>0 && j<N2/2) && (k>0 && k<N3/2))
        {
          unsigned jj=N2-j; unsigned kk=N3-k;

          re(GRF_array[k+N3*(j+N2*0)]) *= sigma;
          im(GRF_array[k+N3*(j+N2*0)]) *= sigma;

          re(GRF_array[kk+N3*(jj+N2*0)])=re(GRF_array[k+N3*(j+N2*0)]);
          im(GRF_array[kk+N3*(jj+N2*0)])=-im(GRF_array[k+N3*(j+N2*0)]);

          ///**********************************************************

          re(GRF_array[k+N3*(jj+N2*0)]) *= sigma;
          im(GRF_array[k+N3*(jj+N2*0)]) *= sigma;

          re(GRF_array[kk+N3*(j+N2*0)])=re(GRF_array[k+N3*(jj+N2*0)]);
          im(GRF_array[kk+N3*(j+N2*0)])=-im(GRF_array[k+N3*(jj+N2*0)]);
        }

        if( (i>0 && i<N1/2) && (j==0) && (k==0))
        {
          unsigned ii=N1-i;

          re(GRF_array[0+N3*(0+N2*i)]) *= sigma;
          im(GRF_array[0+N3*(0+N2*i)]) *= sigma;

          re(GRF_array[0+N3*(0+N2*ii)])=re(GRF_array[0+N3*(0+N2*i)]);
          im(GRF_array[0+N3*(0+N2*ii)])=-im(GRF_array[0+N3*(0+N2*i)]);
          ///**********************************************************

        }

        if( (i==0) && (j==0) && (k>0 && k<N3/2))
        {
          unsigned kk=N3-k;

          re(GRF_array[k+N3*(0+N2*0)]) *= sigma;
          im(GRF_array[k+N3*(0+N2*0)]) *= sigma;

          re(GRF_array[kk+N3*(0+N2*0)])=re(GRF_array[k+N3*(0+N2*0)]);
          im(GRF_array[kk+N3*(0+N2*0)])=-im(GRF_array[k+N3*(0+N2*0)]);
          ///**********************************************************

        }

        if( (i==0) && (j>0 && j<N2/2) && (k==0))
        {
          unsigned jj=N2-j;

          re(GRF_array[0+N3*(j+N2*0)]) *= sigma;
          im(GRF_array[0+N3*(j+N2*0)]) *= sigma;

          re(GRF_array[0+N3*(jj+N2*0)])=re(GRF_array[0+N3*(j+N2*0)]);
          im(GRF_array[0+N3*(jj+N2*0)])=-im(GRF_array[0+N3*(j+N2*0)]);
          ///**********************************************************

        }

        if( (i==0) && (j==0) && (k==0))
        {
          re(GRF_array[0])=0.;
          im(GRF_array[0])=0.;
          ///**********************************************************

        }

        if( (i>0 && i<N1/2) && (j==N2/2) && (k==0))
        {
          unsigned ii=N1-i;

          re(GRF_array[0+N3*(N2/2+N2*i)]) *= sigma;
          im(GRF_array[0+N3*(N2/2+N2*i)]) *= sigma;

          re(GRF_array[0+N3*(N2/2+N2*ii)])=re(GRF_array[0+N3*(N2/2+N2*i)]);
          im(GRF_array[0+N3*(N2/2+N2*ii)])=-im(GRF_array[0+N3*(N2/2+N2*i)]);

          ///**********************************************************

        }

        if( (i>0 && i<N1/2) && (j==0) && (k==N3/2))
        {
          unsigned ii=N1-i;

          re(GRF_array[N3/2+N3*(0+N2*i)]) *= sigma;
          im(GRF_array[N3/2+N3*(0+N2*i)]) *= sigma;

          re(GRF_array[N3/2+N3*(0+N2*ii)])=re(GRF_array[N3/2+N3*(0+N2*i)]);
          im(GRF_array[N3/2+N3*(0+N2*ii)])=-im(GRF_array[N3/2+N3*(0+N2*i)]);

          ///**********************************************************

        }

        if( (i==N1/2) && (j>0 && j<N2/2) && (k==0))
        {
          unsigned jj=N2-j;

          re(GRF_array[0+N3*(j+N2*N1/2)]) *= sigma;
          im(GRF_array[0+N3*(j+N2*N1/2)]) *= sigma;

          re(GRF_array[0+N3*(jj+N2*N1/2)])=re(GRF_array[0+N3*(j+N2*N1/2)]);
          im(GRF_array[0+N3*(jj+N2*N1/2)])=-im(GRF_array[0+N3*(j+N2*N1/2)]);

          ///**********************************************************

        }


        if( (i==0) && (j>0 && j<N2/2) && (k==N3/2))
        {
          unsigned jj=N2-j;

          re(GRF_array[N3/2+N3*(j+N2*0)]) *= sigma;
          im(GRF_array[N3/2+N3*(j+N2*0)]) *= sigma;

          re(GRF_array[N3/2+N3*(jj+N2*0)])=re(GRF_array[N3/2+N3*(j+N2*0)]);
          im(GRF_array[N3/2+N3*(jj+N2*0)])=-im(GRF_array[N3/2+N3*(j+N2*0)]);

          ///**********************************************************
        }

        if( (i==N1/2) && (j==0) && (k>0 && k<N3/2))
        {
          unsigned kk=N3-k;

          re(GRF_array[k+N3*(0+N2*N1/2)]) *= sigma;
          im(GRF_array[k+N3*(0+N2*N1/2)]) *= sigma;

          re(GRF_array[kk+N3*(0+N2*N1/2)])=re(GRF_array[k+N3*(0+N2*N1/2)]);
          im(GRF_array[kk+N3*(0+N2*N1/2)])=-im(GRF_array[k+N3*(0+N2*N1/2)]);

          ///**********************************************************

        }

        if( (i==0) && (j==N2/2) && (k>0 && k<N3/2))
        {
          unsigned kk=N3-k;

          re(GRF_array[k+N3*(N2/2+N2*0)]) *= sigma;
          im(GRF_array[k+N3*(N2/2+N2*0)]) *= sigma;

          re(GRF_array[kk+N3*(N2/2+N2*0)])=re(GRF_array[k+N3*(N2/2+N2*0)]);
          im(GRF_array[kk+N3*(N2/2+N2*0)])=-im(GRF_array[k+N3*(N2/2+N2*0)]);

          ///**********************************************************

        }

        if( (i==0) && (j==0) && (k==N3/2))
        {
          re(GRF_array[N3/2+N3*(0+N2*0)]) *= sqrt(num_2)*sigma;
          im(GRF_array[N3/2+N3*(0+N2*0)])=0.;
          ///**********************************************************

        }

        if( (i==0) && (j==N2/2) && (k==0))
        {
          re(GRF_array[0+N3*(N2/2+N2*0)]) *= sqrt(num_2)*sigma;
          im(GRF_array[0+N3*(N2/2+N2*0)])=0.;
          ///**********************************************************

        }

        if( (i==N1/2) && (j==0) && (k==0))
        {
          re(GRF_array[0+N3*(0+N2*N1/2)]) *= sqrt(num_2)*sigma;
          im(GRF_array[0+N3*(0+N2*N1/2)])=0.;
          ///**********************************************************

        }

        if( (i==N1/2) && (j==N2/2) && (k==0))
        {
          re(GRF_array[0+N3*(N2/2+N2*N1/2)]) *= sqrt(num_2)*sigma;
          im(GRF_array[0+N3*(N2/2+N2*N1/2)])=0.;
          ///**********************************************************

        }

        if( (i==N1/2) && (j==0) && (k==N3/2))
        {
          re(GRF_array[N3/2+N3*(0+N2*N1/2)]) *= sqrt(num_2)*sigma;
          im(GRF_array[N3/2+N3*(0+N2*N1/2)])=0.;
          ///**********************************************************

        }

        if( (i==0) && (j==N2/2) && (k==N3/2))
        {
          re(GRF_array[N3/2+N3*(N2/2+N2*0)]) *= sqrt(num_2)*sigma;
          im(GRF_array[N3/2+N3*(N2/2+N2*0)])=0.;
          ///**********************************************************
        }

      }

  /*
     FFT3d (N1,N2,N3, to_Rspace, GRF_array, GRF_array);

     for (ULONG i=0 ; i<N;i++)
     delta[i] = re(GRF_array[i]);
     */

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for (unsigned i=0 ; i<N1;i++)
    for (unsigned j=0 ; j<N2;j++)
      for (unsigned k=0 ; k<=N3/2;k++)
      {
        ULONG ihalf=k+(N3/2+1)*(j+N2*i);
        ULONG iind=k+N3*(j+N2*i);

        re(GRF_array2[ihalf])=re(GRF_array[iind]);
        im(GRF_array2[ihalf])=im(GRF_array[iind]);
      }

  fftC2R(N1,N2,N3,GRF_array2,delta);

}
