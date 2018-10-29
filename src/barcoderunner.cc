/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include "struct_main.h"
#include "fftw_array.h"

#include <cmath>
#include <iomanip>
#include <cstring>
#include <cassert>
#include <sstream>
#include <algorithm>  // std::max_element

#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>

#include "transf.h"
#include "field_statistics.h"
#include "Lag2Eul.h"
#include "IOfunctions.h"
#include "IOfunctionsGen.h"
#include "math_funcs.h" // power_mean
#include "convolution.hpp" // kernelcomp, convcomp
#include "random.hpp" // create_GARFIELD

#include "protocol.h"
#include "sample_maker.h"

#include "convenience.h"

#include "HMC_models.h" // SPH_kernel_scale
#include "hmc/likelihood/lognormal_independent.hpp" // lognormal_likelihood_f_delta_x_i_calc

using namespace std;

void setup_random_test(struct DATA *data, real_prec *delta_lag, real_prec *delta_eul, unsigned int facL,
                       real_prec *posx, real_prec *posy, real_prec *posz, gsl_rng *gBaseRand)
{
  struct NUMERICAL *n = data->numerical;
  struct OBSERVATIONAL *o = data->observational;
  struct COSMOLOGY *c = data->cosmology;

  fftw_array<real_prec> kmode(n->N_bin), power(n->N_bin);

  // lagrangian
  create_GARFIELD(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, delta_lag, o->Power, gBaseRand);

  //dump_scalar(delta_lag, n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, 0, n->dir + string("deltaLAGtest"));
  measure_spectrum(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, delta_lag, kmode, power, n->N_bin);
  dump_measured_spec(kmode, power, n->dir + string("specLAGtest.dat"), n->N_bin);

  // eulerian
  bool reggrid=true;
  real_prec kernel_scale = SPH_kernel_scale<struct DATA>(data);
  //real_prec kernel_scale;
  //if (n->mk == 3)
    //kernel_scale = n->particle_kernel_h;
  //else
    //kernel_scale = 0;

  if (n->random_test_rsd)
  {
    // cout << " **** random RSD test! **** " << endl;
    Lag2Eul_rsd_zeldovich(delta_lag, delta_eul, posx, posy, posz, n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1,
                          n->d2, n->d3, n->xllc, n->yllc, n->zllc, c->D1, c->ascale, c->omega_m, c->omega_q, n->mk,
                          facL,
                          reggrid, gBaseRand, kernel_scale, n->xobs, n->yobs, n->zobs, n->planepar, n->periodic,
                          n->R2Cplan, n->C2Rplan);
  }
  else
  {
    // cout << " **** normal random test! **** " << endl;
    Lag2Eul(delta_lag, delta_eul, posx, posy, posz, n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1, n->d2, n->d3,
            n->xllc, n->yllc, n->zllc, c->D1, c->D2, c->ascale, c->omega_m, c->omega_q, n->sfmodel, n->mk, n->slength, facL,
            reggrid, gBaseRand, n->dir, kernel_scale, n->R2Cplan, n->C2Rplan);
  }

  //dump_scalar(delta_eul, n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, 0, n->dir + string("deltaEULtest"));
  dump_deltas(data, delta_lag, delta_eul, string("test"));
  //measure_and_dump_spectrum(delta_eul, N1, N2, N3, L1, L2, L3, N_bin, n->dir + string("specEULtest.dat")); // build this function
  measure_spectrum(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, delta_eul, kmode, power, n->N_bin);
  dump_measured_spec(kmode, power,  n->dir + string("specEULtest.dat"), n->N_bin);

  // window
  switch (n->window_type)
  {
    case 1: // fill with ones
      fill_one(o->window, n->N);
      break;
    case 10: // fill half with ones, half with zeros (for testing)
      fill_one(o->window, n->N);
      fillZero(o->window, n->N/2);
      break;
    case 23: // put zeros where delta_eul > 3, ones otherwise
      for (ULONG i = 0; i < n->N; ++i)
      {
        if (delta_eul[i] > 3)
          o->window[i] = 1;
        else
          o->window[i] = 0;
      }
      break;
    default:
      stringstream message;
      message << "in barcoderunner: window_type = " << n->window_type << " is not a valid choice!";
      throw runtime_error(message.str());
  }
  dump_scalar(o->window, n->N1, n->N2, n->N3, n->dir + string("win"));

  // nobs and noise_sf
  switch (n->data_model)
  {
    case 0: // linear data model 
#ifdef MULTITHREAD_RNG
#pragma omp parallel for
#endif // MULTITHREAD_RNG
      for(ULONG i=0;i<n->N;i++)
      {
        real_prec Lambda = o->rho_c * (num_1+delta_eul[i]);
        if (o->window[i]>0.)
        { 
          switch (n->likelihood)
          {
            case 0: // Poissonian likelihood
              o->nobs[i]=static_cast<real_prec>(gsl_ran_poisson(gBaseRand,Lambda));
              break;
            case 1: // Gaussian likelihood
              {
                real_prec sigma_obs = o->sigma_min + o->sigma_fac*Lambda;
                o->noise_sf[i] = sigma_obs;

                real_prec nobs = Lambda + static_cast<real_prec>(GR_NUM(gBaseRand, sigma_obs, 0));
                // Keep above 0 (or not):
                if (!n->negative_obs && nobs < 0) {
                  nobs = 0;
                }
                o->nobs[i] = nobs;
              }
              break;
            case 3: // GRF
              {
                // sigma quadratic in terms of lambda (test)
                real_prec sigma_obs = o->sigma_min + o->sigma_fac * gsl_pow_2(delta_lag[i]);
                o->noise_sf[i] = sigma_obs;

                o->nobs[i] = delta_lag[i] + static_cast<real_prec>(GR_NUM(gBaseRand, sigma_obs, 0));
              }
              break;
            default:
              throw runtime_error("in barcoderunner: linear data model was chosen (additive error), but incompatible likelihood!");
          }
        }
        else
          o->nobs[i] = static_cast<real_prec>(0.0);
      }
      break;
    case 1: // log-normal data model
#ifdef MULTITHREAD_RNG
#pragma omp parallel for
#endif // MULTITHREAD_RNG 
      for (ULONG i = 0; i < n->N; i++)
      {
        real_prec delta_i = delta_eul[i];
        real_prec Lambda = lognormal_likelihood_f_delta_x_i_calc(o->rho_c, o->delta_min, delta_i);

        if (o->window[i]>0.)
        { 
          real_prec sigma_obs = o->sigma_fac;
          o->noise_sf[i] = sigma_obs;

          real_prec nobs = Lambda + static_cast<real_prec>(GR_NUM(gBaseRand, sigma_obs, 0));
          o->nobs[i] = nobs;
        }
        else
          o->nobs[i] = static_cast<real_prec>(log(gsl_pow_2(o->rho_c * (1 + o->delta_min)))); // huh, why squared?
      }
      break;
    default:
      stringstream message;
      message << "in barcoderunner: data_model = " << n->data_model << " is not a valid choice!";
      throw runtime_error(message.str());
  }

  // Gaussian (RSD) and GRF likelihood can't have zero noise terms, or you'll get NaNs.
  if ((n->likelihood == 1) || (n->likelihood == 3))
    for (ULONG i = 0; i < n->N; i++)
      if (o->window[i] > 0)
        if (o->noise_sf[i] == 0.) {
          stringstream message;
          message << "in barcoderunner(): noise = 0 found! Index " << i;
          throw runtime_error(message.str());
        }

  dump_scalar(o->nobs, n->N1, n->N2, n->N3, n->dir + string("nobs"));
  dump_scalar(o->noise_sf, n->N1, n->N2, n->N3, n->dir + string("sigma"));

  measure_spectrum(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3,o->nobs,kmode,power, n->N_bin);
  dump_measured_spec(kmode,power, n->dir + string("spec_nobs.dat"), n->N_bin);
}

void make_initial_guess(struct DATA *data, gsl_rng *gBaseRand)
{
  struct NUMERICAL *n = data->numerical;
  struct OBSERVATIONAL *o = data->observational;
  // Setting the initial guess in Lagrangian (signal) space 
  switch (n->initial_guess)
  {
    case 0: // zero initial guess 
      fillZero(o->signal, n->N);
      break;
    case 1: // initial guess from file 
      get_scalar(n->dir + n->initial_guess_file, o->signal, n->N1, n->N2, n->N3);
      break;
    case 2: // GRF initial guess 
      create_GARFIELD(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, o->signal, o->Power, gBaseRand);
      break;
    case 3: // smoothed GRF initial guess 
      {
        create_GARFIELD(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, o->signal, o->Power, gBaseRand);
        int filtertype = n->initial_guess_smoothing_type;
        real_prec smoothing_length = n->initial_guess_smoothing_scale;
        kernelcomp(n->L1, n->L2, n->L3, n->N1, n->N2, n->N3, smoothing_length, filtertype, data);
        convcomp(n->N1, n->N2, n->N3, o->signal, o->signal, smoothing_length, n->dir);
      }
      break;
    case 4: // zero plus some random noise
      {
        fillZero(o->signal, n->N);
#ifdef MULTITHREAD_RNG
#pragma omp parallel for
#endif // MULTITHREAD_RNG
        for (ULONG i = 0; i < n->N; ++i)
          o->signal[i] += static_cast<real_prec>(gsl_ran_gaussian(gBaseRand, 1.e-1));
      }
      break;
    default:
      stringstream message;
      message << "In barcoderunner: invalid choice of initial_guess (" << n->initial_guess << ")!";
      throw runtime_error(message.str());
  }
}

unsigned int initial_iteration_number(struct DATA *data)
{
  struct NUMERICAL *n = data->numerical;

  unsigned resnum = 0;
#ifdef RESTART_FILE
  string fname= n->dir + string("restart");
  ifstream inStream;
  inStream.open(fname.data());

  if (inStream.is_open() == true )
  {
    cout<<endl;
    cout<<"attention: restart command!"<<endl;

    string fname2= n->dir + string("restart.prt");
    ifstream inStream2;
    inStream2.open(fname2.data());
    assert(inStream2.is_open());

    inStream2 >> resnum;

    inStream2.close();
#else // RESTART_FILE, i.e. the following is for when RESTART_FILE is not defined:
  if (n->start_at > 0)
  {
    resnum = n->start_at;
#endif // RESTART_FILE
    cout<<"restart file: "<<resnum<<endl;
  }

  return(resnum);
}


void load_initial_fields(struct DATA *data, ULONG resnum, gsl_rng *gBaseRand)
{
  struct NUMERICAL *n = data->numerical;
  struct OBSERVATIONAL *o = data->observational;

  if (resnum > 0) // restart from existing snapshots
  {
    {
      string fname_res= n->dir + string("deltaLAG");
      int bmax_res=100;
      char buffer_res[bmax_res];
      sprintf(buffer_res,"_%lu",resnum);
      string FileName_res=string(fname_res)+buffer_res;
      char * file_res;
      file_res = new char[FileName_res.length() + 1];
      strcpy(file_res, FileName_res.c_str());
      get_scalar(file_res, o->signal, n->N1, n->N2, n->N3);
      delete[] file_res;
    }

    get_scalar( n->dir + string("win"), o->window, n->N1, n->N2, n->N3);
    get_scalar( n->dir + string("nobs"), o->nobs, n->N1, n->N2, n->N3);
    get_scalar( n->dir + string("sigma"),  o->noise_sf,  n->N1,  n->N2,  n->N3);
  }
  else
  {
    if (n->random_test)
    {
      fftw_array<real_prec> delta_lag(n->N), delta_eul(n->N);
      unsigned facL=1; // if you change this, you have to change the size of the pos arrays as well
      fftw_array<real_prec> posx(n->N),posy(n->N),posz(n->N);
      setup_random_test(data, delta_lag, delta_eul, facL, posx, posy, posz, gBaseRand);
    }
    else
    {
      get_scalar( n->dir + string("win"),  o->window, n->N1, n->N2, n->N3);
      get_scalar( n->dir + string("nobs"),  o->nobs, n->N1, n->N2, n->N3);
      get_scalar( n->dir + string("sigma"),  o->noise_sf, n->N1, n->N2, n->N3);
    }

    make_initial_guess(data, gBaseRand); // in lagrangian (signal) space
    // dump initial guess:
    dump_scalar(o->signal, n->N1, n->N2, n->N3, n->dir + string("initial_guess"));
    fftw_array<real_prec> kmode(n->N_bin), power(n->N_bin);
    measure_spectrum(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, o->signal,kmode,power,n->N_bin);
    dump_measured_spec(kmode,power, n->dir + string("spec_initial_guess.dat"), n->N_bin);
  }

  // The signal in Eulerian space doesn't need to be set, the
  // first time it will be used it will be filled with Lag2Eul first.
  // But, let's set it to zero, just for fun.
  fillZero( o->signalX, n->N);
  // EGP: that was fun!

  // Do remember to do the above when you need it! For instance, I forgot to do
  // Lag2Eul inside the likeli_force_1st_order_diagonal_mass function and before
  // that function no Lag2Eul calls have been made yet and posx,y,z were set
  // initially to zero, so that gave crap results.
  // I considered filling the arrays correctly here already, but actually the
  // pos arrays are not passed along to the HAMIL_DATA struct.
}


void initialize_performance_log(struct DATA *data, ULONG resnum)
{
  if (resnum == 0)
    data->numerical->performance_log.open("performance_log.txt");
  else
    data->numerical->performance_log.open("performance_log.txt", ofstream::out | ofstream::app); // "out" for writing (must be specified), "app" for appending instead of overwriting (necessary at restart)

  assert(data->numerical->performance_log.is_open());
  if (resnum == 0)
  {
    data->numerical->performance_log << "accepted\tepsilon\tNeps\tdH\tdK\tdE\tdprior\tdlikeli\t";
    data->numerical->performance_log << "psi_prior_i\tpsi_prior_f\tpsi_likeli_i\tpsi_likeli_f\tH_kin_i\tH_kin_f" << endl;
  }
}



void barcoderunner(struct DATA *data,gsl_rng * gBaseRand)
{
  struct NUMERICAL *n = data->numerical;
  struct OBSERVATIONAL *o = data->observational;

  // prepare smoothing kernels and transfer functions
  {
    real_prec kthsc=n->slength;//attention!!
    int ftype=1;//attention!!

    kernelcomp(n->L1, n->L2, n->L3, n->N1, n->N2, n->N3, n->slength, ftype, data);

    if (kthsc!=n->slength)
      kernelcomp(n->L1, n->L2, n->L3, n->N1, n->N2, n->N3, kthsc, ftype, data);

#ifdef TRANSF
    int sftype;
    switch (n->sfmodel)
    {
      case 0:
        {
          cout<<"does not apply"<<endl;
        }
        break;
      case 1:
        {
          sftype=1;
          transflpt(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1, n->d2, n->d3,n->slength, sftype, n->dir);
        }
        break;
      case 2: // same as case 3
      case 3:
        {
          sftype=1;
          transflpt(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1, n->d2, n->d3,n->slength, sftype, n->dir);
          sftype=2;
          transflpt(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1, n->d2, n->d3,n->slength, sftype, n->dir);
        }
        break;
    }
#endif
  }
  // END prepare smoothing kernels and transfer functions

  unsigned resnum = initial_iteration_number(data); // restart number (0 if not restarting)
  unsigned iistart = resnum + 1;

  // read out possibly pre-existing log to set rejections and epsilon back to
  // what they were before restarting
  ifstream plog_read ("performance_log.txt");
  auto plog_line_count = count(istreambuf_iterator<char>(plog_read),
                               istreambuf_iterator<char>(), '\n');
  plog_read.seekg(0, plog_read.beg);  // move "cursor" back to start of stream

  if (plog_line_count > 0 && resnum > 0)
  {
    cout << "Performance log found with " << plog_line_count << " lines. Reading..." << endl;
    unsigned accepted = 0; // counter for accepted steps
    unsigned acc; // dummy
    string line;
    std::vector<std::string> plog_row;
    ULONG ix_leapfrog;
    real_prec epsilon_convert;
    while ( std::getline(plog_read, line) && (accepted < resnum) )
    {
      plog_row = split(line, '\t');
      std::istringstream acc_ss(plog_row[0]);
      acc_ss >> acc;
      if (acc == 0)
        ++n->rejections;
      else
      {
        ++accepted;
      }

      // update epsilon circular buffer
      ix_leapfrog = n->rejections + accepted - 1;
      std::istringstream eps_ss(plog_row[1]);
      eps_ss >> epsilon_convert;
      n->epsilon_N_a[ix_leapfrog % n->N_a_eps_update] = epsilon_convert;
    }

    std::cout << "accepted/rejected steps reset to " << accepted
              << "/" << n->rejections << std::endl;

    // reset eps_fac to what it was before restart
    switch (n->eps_fac_update_type) {
      case 0:
      break;
      case 1:
      {
        ULONG total_steps_before_restart = resnum + n->rejections;
        ULONG updates = total_steps_before_restart / n->s_eps_total;
        for (ULONG i = 0; i < updates; ++i)
        {
          n->eps_fac = power_mean(n->eps_fac, n->eps_fac_target, n->eps_fac_power);
          cout << "  updating eps_fac to " << n->eps_fac << endl;
        }
      }
      break;
      case 2:
      case 3:
      {
        // Not actually exact, but probably close enough. Will give a slightly
        // smaller epsilon initially, leading to a slightly higher initial
        // acceptance rate than before the restart. Already a lot better than
        // wasting hours due to ~0 initial acceptance rate.
        n->eps_fac = *std::max_element(n->epsilon_N_a.begin(),
                                       n->epsilon_N_a.end());
        std::cout << "Maximum epsilon (eps_fac) reset to " << n->eps_fac << std::endl;
      }
      break;
    }

  }
  plog_read.close();

  load_initial_fields(data, resnum, gBaseRand);

  INIT_PROTOCOL_CONV(data);
  INIT_PROTOCOL_SPEC(data);

  fftw_array<real_prec> kmode(n->N_bin), power(n->N_bin);

  initialize_performance_log(data, resnum);

  // start MCMC iterations 
  for (ULONG ii=iistart;ii<=n->N_Gibbs;ii++)
  {
    // if (n->N_Gibbs>0)
      // cout<<"\n"<<">>> MCMC-sampling iterations: "<<ii<<endl;

    // set MCMC iteration counter

    n->iGibbs=ii;

    string fnamefe= n->dir + string("fastexit");
    ifstream inStreamfe;
    inStreamfe.open(fnamefe.data());

    if (inStreamfe.is_open())
    {
      stringstream message;
      message << "attention: fast exit command!";
      throw runtime_error(message.str());
    }
    inStreamfe.close();

    bool write_output = 0==ii%n->outnum || ii <= 10;
    bool write_ps_output = 0==ii%n->outnum_ps || ii <= 10;

#ifdef RESTART_FILE
    if (( write_output || ii==1) && ii>1)
      PROTOCOL_RESTART(ii-1, data);
#endif // RESTART_FILE

    sample_maker(data,gBaseRand);

    if (write_output)
    {
      int bmax = 100;
      char buffer[bmax];  	    
      sprintf(buffer, "_%lu", ii);
      dump_deltas(data, o->signal, o->signalX, buffer);
    }

    if (write_ps_output)
    {
      measure_spectrum(n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, o->signal, kmode, power, n->N_bin);
      dump_ps_it(data,kmode,power);
    }

    // ULONG total_steps = n->iGibbs + n->rejections;

    // real_prec accepted_steps_f = static_cast<real_prec>(n->iGibbs);
    // real_prec total_steps_f = static_cast<real_prec>(total_steps);
    // real_prec acceptance_rate = accepted_steps_f/total_steps_f;

    // cout << "***** Current acceptance rate: " << acceptance_rate << endl;
  }
}
