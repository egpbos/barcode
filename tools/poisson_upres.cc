/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <iostream>

#include "define_opt.h"
#include "struct_main.h"
#include "fftw_array.h"
#include <cmath>
#include <sstream>
#include <gsl/gsl_randist.h>
#include "IOfunctionsGen.h"

#include "massFunctions.h" // getDensity

using namespace std;

vector< vector<real_prec> > discrete_poisson_sample(real_prec *lambda, unsigned N1, unsigned N2, unsigned N3, real_prec L1, real_prec L2, real_prec L3, gsl_rng *gBaseRand)
{
  real_prec d1 = L1/static_cast<real_prec>(N1);
  real_prec d2 = L2/static_cast<real_prec>(N2);
  real_prec d3 = L3/static_cast<real_prec>(N3);

  ULONG N = static_cast<ULONG>(N1) * static_cast<ULONG>(N2) * static_cast<ULONG>(N3);

  fftw_array<ULONG> number(N);
  ULONG Npart = 0;

  cout << "... drawing particle numbers per cell" << endl;
  for (ULONG i = 0; i < N; ++i)
  {
    number[i] = static_cast<ULONG>( floor(static_cast<real_prec>( gsl_ran_poisson(gBaseRand, lambda[i]) )) );
    Npart += number[i];
  }

  vector< vector<real_prec> > pos_pois(Npart, vector<real_prec>(3));
  ULONG part_ix = 0;

  cout << "... creating particle positions" << endl;
  for (unsigned i=0;i<N1;i++)
    for (unsigned j=0;j<N2;j++)
      for (unsigned k=0;k<N3;k++)
      {
        ULONG lc=k+N3*(j+N2*i);                     
        for(ULONG kk=0; kk < number[lc]; kk++)
        {
          real_prec rx = static_cast<real_prec>(gsl_rng_uniform(gBaseRand)) * d1;
          real_prec ry = static_cast<real_prec>(gsl_rng_uniform(gBaseRand)) * d2;
          real_prec rz = static_cast<real_prec>(gsl_rng_uniform(gBaseRand)) * d3;

          pos_pois[part_ix][0] = d1 * (static_cast<real_prec>(i)) + rx;
          pos_pois[part_ix][1] = d2 * (static_cast<real_prec>(j)) + ry;
          pos_pois[part_ix][2] = d3 * (static_cast<real_prec>(k)) + rz;

          ++part_ix;
        }
      }

  return (pos_pois);
}

void load_arguments(int argc, char *argv[], string &fname_in, unsigned &N1, real_prec &L1, unsigned &N1_out, real_prec &Nbar, ULONG &seed, string &fname_out)
{
  stringstream N1_arg, L1_arg, N1_out_arg, Nbar_arg, seed_arg;

  if (argc >= 7)
  {
    fname_in = string(argv[1]);
    N1_arg << argv[2];
    N1_arg >> N1;
    L1_arg << argv[3];
    L1_arg >> L1;
    N1_out_arg << argv[4];
    N1_out_arg >> N1_out;
    Nbar_arg << argv[5];
    Nbar_arg >> Nbar;
    seed_arg << argv[6];
    seed_arg >> seed;

    if (argc >= 8)
      fname_out = string(argv[7]);
    else
      fname_out = fname_in + string("_poisCIC") + N1_out_arg.str() + string("_Nbar") + Nbar_arg.str();

    cout << "will output to file: " << fname_out << endl;
  }
  else
  {
    cerr << "Need 6 parameters (file in, N1, L1, N1_out, Nbar, seed)! " << endl
         << "N.B.: filenames must be given without extension (which must be "
            ".dat)." << endl
         << "Optional: fname_out (default: "
            "fname_in+'_poisCIC[N1_out]_Nbar[Nbar]')." << endl << endl
         << "Nbar is the number of poisson particles per cell. For a decent "
            "representation of the previous density field, an Nbar of "
            "8^log2(N1_out/N1) seems to work well."
         << endl;
    exit(1);
  }
}

int main(int argc, char *argv[])
{
  unsigned N1, N1_out;
  ULONG seed;
  real_prec L1, Nbar;
  string fname_in, fname_out;
  load_arguments(argc, argv, fname_in, N1, L1, N1_out, Nbar, seed, fname_out);

  ULONG N = N1*N1*N1;
  ULONG N_out = N1_out*N1_out*N1_out;

  //real_prec d1 = L1 / static_cast<real_prec>(N1);
  real_prec d1_out = L1 / static_cast<real_prec>(N1_out);

  fftw_array<real_prec> input(N);

  // get data from input file
  get_scalar(fname_in, input, N1, N1, N1);
  // convert input from delta to lambda
  for (ULONG i = 0; i < N; ++i)
    input[i] = Nbar * (1 + input[i]);

  // initialize RNG
  gsl_rng *gBaseRand;
  gBaseRand = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (gBaseRand, seed);

  // do poisson sample; get positions
  vector< vector<real_prec> > pos_pois = discrete_poisson_sample(input, N1, N1, N1, L1, L1, L1, gBaseRand);
  // convert to usable format
  ULONG Npart = pos_pois.size();
  fftw_array<real_prec> posx(Npart), posy(Npart), posz(Npart), pmass(Npart);
  for (ULONG i = 0; i < Npart; ++i)
  {
    posx[i] = pos_pois[i][0];
    posy[i] = pos_pois[i][1];
    posz[i] = pos_pois[i][2];
    //pmass[i] = 1.; // particle mass (not really useful, but necessary for getDensity)
  }

  // calculate new density (CIC))
  cout << "... calculate new density from particle sample" << endl;
  fftw_array<real_prec> output(N_out);
  bool pmass_used = false;
  getDensity_CIC(N1_out, N1_out, N1_out, L1, L1, L1, d1_out, d1_out, d1_out, 0, 0, 0, posx, posy, posz, pmass, Npart, output, pmass_used);

  // output results
  quick_dump_scalar(output, N1_out, fname_out, 0, false);

  return(0);
}
