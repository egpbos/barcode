/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include <fftw3.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>

#include <fstream>  //NOLINT
#include <iomanip>
#include <cassert>
#include <ncurses.h>

#include "struct_main.h"
#include "struct_hamil.h"
#include "fftw_array.h"
#include "math_funcs.h"
#include "IOfunctions.h"
#include "Lag2Eul.h"
#include "field_statistics.h"
#include "convenience.h"

#include "HMC_mass.h"

using namespace std;


void likeli_force_power(struct HAMIL_DATA *hd, real_prec *signal,
  real_prec *likeli_power, real_prec *kmode, struct DATA * data) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  fftw_array<real_prec> likeli_force(n->N);

  likelihood_grad_log_like(hd, signal, likeli_force);

  measure_spectrum(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, likeli_force,
    kmode, likeli_power, n->N_bin);
  dump_measured_spec(kmode, likeli_power,
    data->numerical->dir + string("forcespec.dat"), n->N_bin);
}


void Hamiltonian_mass_likeli_force(struct HAMIL_DATA *hd, real_prec *signal,
  struct DATA *data, real_prec *out) {
  struct HAMIL_NUMERICAL *n = hd->numerical;

  real_prec likeli_power[n->N_bin];
  real_prec kmode[n->N_bin];
  likeli_force_power(hd, signal, likeli_power, kmode, data);

  real_prec kmax = sqrt(k_squared(n->N1/2, n->N2/2, n->N3/2, n->L1, n->L2,
    n->L3, n->N1, n->N2, n->N3));
  real_prec dk = kmax/static_cast<real_prec>(n->N_bin);

  real_prec NORM = num_1;  // [>static_cast<real_prec>(N);// care
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif  // MULTITHREAD
  for (unsigned i = 0; i < n->N1; i++)
    for (unsigned j = 0; j < n->N2; j++)
      for (unsigned k = 0; k < n->N3; k++) {
        ULONG l = k+n->N3*(j+n->N2*i);

        real_prec kr = sqrt(k_squared(i, j, k, n->L1, n->L2, n->L3, n->N1,
          n->N2, n->N3));
        ULONG nbin = static_cast<ULONG>(kr/dk);

        if (kr>0.)
          out[l] = likeli_power[nbin]*NORM;
        else
          out[l] = 0.0;
      }
}


real_prec Hamiltonian_mass_mean_likeli_force(struct HAMIL_DATA *hd,
  real_prec *signal, struct DATA *data) {
  struct HAMIL_NUMERICAL *n = hd->numerical;

  real_prec likeli_power[n->N_bin];
  real_prec kmode[n->N_bin];
  likeli_force_power(hd, signal, likeli_power, kmode, data);

  real_prec kmax = sqrt(k_squared(n->N1/2, n->N2/2, n->N3/2, n->L1, n->L2,
    n->L3, n->N1, n->N2, n->N3));
  real_prec dk = kmax/static_cast<real_prec>(n->N_bin);

  real_prec force_mean_1D = 0.0;
  real_prec k_volume_1D = 0.0;
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:force_mean_1D)
  #endif  // MULTITHREAD
  for (ULONG i = 0; i < n->N_bin; i++)
    force_mean_1D += 4.*M_PI * kmode[i] * kmode[i] * dk * likeli_power[i];
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:k_volume_1D)
  #endif  // MULTITHREAD
  for (ULONG i = 0; i < n->N_bin; i++)
    k_volume_1D += 4.*M_PI * kmode[i] * kmode[i] * dk;

  force_mean_1D /= k_volume_1D;

  return(force_mean_1D);
}


real_prec inv_ps(struct HAMIL_DATA *hd, ULONG ix) {
  real_prec invP;
  if (hd->signal_PS[ix] > 0.0)
    invP = 1./hd->signal_PS[ix];
  else
    invP = 0.;
  return (invP);
}


void likeli_force_mass(struct HAMIL_DATA *hd, real_prec *signal,
  struct DATA *data, real_prec fac_conv) {
  // invP + (mean) likelihood-force mass
  struct HAMIL_NUMERICAL *n = hd->numerical;

  Hamiltonian_mass_likeli_force(hd, signal, data, hd->mass_f);

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif  // MULTITHREAD
  for (ULONG i = 0; i < n->N; ++i) {
    real_prec invP = inv_ps(hd, i);

    hd->mass_f[i] = fac_conv * ( 2*invP + sqrt(invP * hd->mass_f[i]) );
  }
}


void mean_likeli_force_mass(struct HAMIL_DATA *hd, real_prec *signal,
  struct DATA *data, real_prec fac_conv) {
  // invP + (mean) likelihood-force mass
  struct HAMIL_NUMERICAL *n = hd->numerical;

  real_prec force_mean_1D = Hamiltonian_mass_mean_likeli_force(hd, signal, data);

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif  // MULTITHREAD
  for (ULONG i = 0; i < n->N; ++i) {
    real_prec invP = inv_ps(hd, i);

    hd->mass_f[i] = fac_conv * ( 2*invP + sqrt(invP * force_mean_1D) );
  }
}


void inverse_power_spectrum_mass(struct HAMIL_DATA *hd, real_prec fac_conv) {
  // invP + (mean) likelihood-force mass
  struct HAMIL_NUMERICAL *n = hd->numerical;

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif  // MULTITHREAD
  for (ULONG i = 0; i < n->N; ++i)
    hd->mass_f[i] = fac_conv * inv_ps(hd, i);
}





// 1st order likelihood force expansion (like Jasche+13) //
void Wprime_il(struct HAMIL_DATA *hd, ULONG l, real_prec *out_x,
  real_prec *out_y, real_prec *out_z) {
  struct HAMIL_NUMERICAL *n = hd->numerical;

  // cell center for cell l
  ULONG N3 = static_cast<ULONG>(n->N3);
  ULONG N2 = static_cast<ULONG>(n->N2);
  int ix_x = static_cast<int>((l / N3) / N2);
  int ix_y = static_cast<int>((l / N3) % N2);
  int ix_z = static_cast<int>(l % N3);
  real_prec xl = (static_cast<real_prec>(ix_x) + 0.5) * n->d1;
  real_prec yl = (static_cast<real_prec>(ix_y) + 0.5) * n->d2;
  real_prec zl = (static_cast<real_prec>(ix_z) + 0.5) * n->d3;

  real_prec h = n->particle_kernel_h;  // SPH scale length

  real_prec norm = 1. / (M_PI * gsl_pow_5(h));

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif  // MULTITHREAD
  for (ULONG i = 0; i < n->N; ++i) {
    real_prec dx = hd->posx[i] - xl;
    real_prec dy = hd->posy[i] - yl;
    real_prec dz = hd->posz[i] - zl;
    pacman_difference(&dx, n->L1);
    pacman_difference(&dy, n->L2);
    pacman_difference(&dz, n->L3);

    real_prec r = sqrt(dx * dx + dy * dy + dz * dz);
    real_prec q = r / h;

    if (q >= 2) {  // check q >= 2 first, because it will be true most cases
      out_x[i] = 0;
      out_y[i] = 0;
      out_z[i] = 0;
    } else if (q >= 1) {  // after that q >= 1 will be most prevalent
      real_prec common = norm * (3 - 0.75 * q - 3./q);
      out_x[i] = dx * common;
      out_y[i] = dy * common;
      out_z[i] = dz * common;
    } else {  // don't need to check for q > 0, sqrt(sq) in r keeps it positive
      real_prec common = norm * (2.25 * q - 3);
      out_x[i] = dx * common;
      out_y[i] = dy * common;
      out_z[i] = dz * common;
    }
  }
}


void likeli_force_1st_order_diagonal_mass(struct HAMIL_DATA *hd, real_prec *signal, struct DATA *data) {
  // C_ih term, diagonal only (C_ii) //
  // In this version, C_ii is the whole mass, so don't need extra C_ii array.

  // 0. prepare data structures
  struct HAMIL_NUMERICAL *n = hd->numerical;
  fftw_array<real_prec> Wprime_il_x(n->N), Wprime_il_y(n->N), Wprime_il_z(n->N);
  fftw_array<complex_prec> Wprime_il_A_f(n->Nhalf), C_ii_l_f(n->Nhalf);

  fillZero(hd->mass_r, n->N);

  // also have to run Lag2Eul, because pos arrays are needed in Wprime_il
  {
    unsigned facL=1;
    bool reggrid=true;
    gsl_rng *seed = nullptr; // empty: reggrid is true anyway
    real_prec kernel_scale = SPH_kernel_scale(hd);
    if (hd->rsd_model)
      Lag2Eul_rsd_zeldovich(signal, hd->deltaX, hd->posx, hd->posy, hd->posz, n->N1, n->N2, n->N3, n->L1, n->L2, n->L3,
                            n->d1, n->d2, n->d3, n->min1, n->min2, n->min3, hd->D1, hd->ascale, hd->OM, hd->OL, n->mk,
                            facL,
                            reggrid, seed, kernel_scale, n->xobs, n->yobs, n->zobs, n->planepar, n->periodic,
                            n->R2Cplan, n->C2Rplan);
    else
      Lag2Eul(signal, hd->deltaX, hd->posx, hd->posy, hd->posz, n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1, n->d2,
              n->d3, n->min1, n->min2, n->min3, hd->D1, hd->D2, hd->ascale, hd->OM, hd->OL, data->numerical->sfmodel, n->mk, n->kth, facL,
              reggrid, seed, "", kernel_scale, n->R2Cplan, n->C2Rplan);
  }

  bool rfft = true;  // use real FFT's to cut calculation time in half
  // N.B.: the normalisation is (-1)(-1)m^2. -m^2 comes from C_ih itself,
  // the other minus comes from the fact that M_ih = S^-1_ih - C_ih, so there's
  // an extra minus there!
  real_prec minus_minus_SPH_mass_sq = gsl_pow_2(hd->rho_c * n->vol /
    static_cast<real_prec>(n->N));

  // boost::progress_display show_progress(n->N);

  for (ULONG l = 0; l < n->N; ++l) {  // 1.
    if (hd->window[l] > 0) {  // 2.
      Wprime_il(hd, l, Wprime_il_x, Wprime_il_y, Wprime_il_z);  // 4.

      fftR2C(n->N1, n->N2, n->N3, Wprime_il_x, Wprime_il_A_f);  // 5. x
      grad_inv_lap_FS(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, Wprime_il_A_f,
        C_ii_l_f, 1, rfft);  // 6. calculate -ik_x/k^2 ^W'_il_x, store in C_ii_l

      fftR2C(n->N1, n->N2, n->N3, Wprime_il_y, Wprime_il_A_f);  // 5. y
      grad_inv_lap_FS(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, Wprime_il_A_f,
        Wprime_il_A_f, 2, rfft);  // 6. calculate -ik_y/k^2 ^W'_il_y, in-place
      add_to_array(Wprime_il_A_f, C_ii_l_f, n->Nhalf);  // 7. sum to temp array

      fftR2C(n->N1, n->N2, n->N3, Wprime_il_z, Wprime_il_A_f);  // 5. z
      grad_inv_lap_FS(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, Wprime_il_A_f,
        Wprime_il_A_f, 3, rfft);  // 6. calculate -ik_z/k^2 ^W'_il_z, in-place
      add_to_array(Wprime_il_A_f, C_ii_l_f, n->Nhalf);  // 7. sum to temp array

      // 8. IFFT temporary array from 7, reusing Wprime_il_x for the result
      // ("div_inv_lap_Wprime_il" if we had to name it):
      fftC2R(n->N1, n->N2, n->N3, C_ii_l_f, Wprime_il_x);

      // 9. take the square of all elements, 10. multiply by w_l/sigma_l^2,
      // 11. sum elements to the appropriate C_ii, all in one loop-step:
      #ifdef MULTITHREAD
      #pragma omp parallel for
      #endif  // MULTITHREAD
      for (ULONG i = 0; i < n->N; ++i)
        // C_ii[i] += hd->window[i] * gsl_pow_2(Wprime_il_A[i] / hd->noise[i]);
        hd->mass_r[i] += hd->window[i] * gsl_pow_2(Wprime_il_x[i] / hd->noise[i]);
    }
    // ++show_progress;
    double progress_percent = 100. * static_cast<double>(l) / static_cast<double>(n->N);
    wprintw(data->curses->status, "computing Hamiltonian mass... %5.1d\%", progress_percent);
    wrefresh(data->curses->status);
  }
  // 12.
  multiply_factor_array(minus_minus_SPH_mass_sq, hd->mass_r, hd->mass_r, n->N);
}







// Momenta: mass for Gaussian likelihood! //
void Hamiltonian_mass(struct HAMIL_DATA *hd, real_prec *signal, struct DATA *data) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  // FT[M] = 1/P+re(FT[w*Nmean])

  wprintw(data->curses->status, "computing Hamiltonian mass...");
  wrefresh(data->curses->status);

  real_prec fac_conv = 1;

  switch (n->mass_type) {
    case 0:  // all ones (R)
      fill_one(hd->mass_r, n->N);
      break;
    case 1:  // inverse power spectrum (FS)
      inverse_power_spectrum_mass(hd, fac_conv);
      break;
    case 2:  // inverse power spectrum + likelihood force (FS)
      likeli_force_mass(hd, signal, data, fac_conv);
      break;
    case 3:  // inverse power spectrum + *mean* likelihood force (Wang+12) (FS)
      mean_likeli_force_mass(hd, signal, data, fac_conv);
      break;
    case 4:  // mass == P(k) (FS)
      copyArray(hd->signal_PS, hd->mass_f, n->N);
      break;
    case 5:  // inverse power spectrum (FS) + 1st order likelihood force
             // expansion (Jasche+13) (R)
      inverse_power_spectrum_mass(hd, fac_conv);
      likeli_force_1st_order_diagonal_mass(hd, signal, data);
      break;
    case 6:  // 1st order likelihood force expansion (Jasche+13) (R)
      likeli_force_1st_order_diagonal_mass(hd, signal, data);
      break;
    case 60:  // type 0 until burn-in, type 6 afterwards
      // N.B.: n->s_eps_total is quite random, but it should not be a lot more.
      // The amount of steps s is also rejected steps for the epsilon stepping,
      // so probably if you count by n->iGibbs, you'll be a long way into the
      // burn-in if you take e.g. a multiple of s_eps_total.
      if (n->iGibbs < n->s_eps_total) {
        fill_one(hd->mass_r, n->N);
      } else {
        likeli_force_1st_order_diagonal_mass(hd, signal, data);
      }
      break;
  }

  // Testing:
  if (n->mass_fs)
    multiply_factor_array(n->mass_factor, hd->mass_f, hd->mass_f, n->N);

#ifdef DEBUG
  dump_mass_spec(hd, data);
#endif  // DEBUG
}









// helper function //
void dump_mass_spec(struct HAMIL_DATA *hd, struct DATA *data) {
  struct HAMIL_NUMERICAL *n = hd->numerical;

  int bmax = 100;
  char buffer1[bmax];
  sprintf(buffer1, "mass_spec.dat");
  string outputFileName = data->numerical->dir + string(buffer1);
  ofstream outStream(outputFileName.data());
  assert(outStream.is_open());
  real_prec NORM = 1.0;

  real_prec kr_old = 0.;
  for (unsigned i = 0; i < n->N1; i++)
    for (unsigned j = 0; j < n->N2; j++)
      for (unsigned k = 0; k < n->N3; k++) {
        real_prec k2 = k_squared(i, j, k, n->L1, n->L2, n->L3, n->N1, n->N2,
          n->N3);
        real_prec kr = sqrt(k2);

        if (i == 0 && j == 0 && k == 0)
          hd->mass_f[k+n->N3*(j+n->N2*i)]=1.;

        if (kr > kr_old) {
          kr_old = kr;
          outStream << kr << "   " << NORM*hd->mass_f[k+n->N3*(j+n->N2*i)]
            << endl;
        }
      }
  outStream.close();
}


// ABOUT fac_conv IN Hamiltonian_mass:
//   I was experimenting with this factor, but didn't really find a suitable
//   conclusion. Below is what I previously wrote down there. In the function
//   let's stick with fac_conv = 1 for now.
// EGP: a convolution factor is needed here, because in reality what we have
// is that the "operators" M and P must be each others' inverses. This means
// that their (continuous) convolution must give a Dirac delta. However,
// since we want a discrete version of this mass, we need to use the
// convolution factor V or V/N (see Martel 2005).
// N.B.: not sure if it must also be in front of the Wang/force-masses...
// real_prec fac_conv;
// #ifdef FOURIER_DEF_1
// fac_conv = 1./n->vol;
// #endif
// #ifdef FOURIER_DEF_2
// fac_conv = static_cast<real_prec>(n->N)/n->vol;
// #endif
// EGP:
// Possibly, another factor is needed though. And maybe the above isn't
// necessary at all... However, one can imagine that the V and N factors
// from the power spectrum need to be canceled first before adding new
// ones... I've got it on paper.
//#ifdef FOURIER_DEF_2
//fac_conv = n->vol/static_cast<real_prec>(n->N*n->N*n->N);
//#endif
// EGP: definitely, this should not be in front the Wang/force-mass terms...
// anyway... this term includes the N/V factor. Another possibility would be:
//#ifdef FOURIER_DEF_2
//fac_conv = gsl_pow_2(n->vol) / static_cast<real_prec>(n->N*n->N*n->N*n->N);
//#endif
// EGP: doesn't seem like this is going to work. Let's try the old one Paco used
// (which was just N, but now the GRF function is changed, so we compensate for
// that with V/N^2):
//#ifdef FOURIER_DEF_2
//fac_conv = static_cast<real_prec>(n->N)*n->vol/static_cast<real_prec>(n->N*n->N); // yeah, that's really just V/N, but for clarity of the derivation, let's keep it like this.
//#endif
// EGP: okay, let's just test how the mass behaves when we vary V and N in input.par.
// fac_conv = 1.;

