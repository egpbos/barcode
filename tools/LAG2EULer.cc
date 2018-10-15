/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include "define_opt.h"
#include "struct_main.h"
#include "fftw_array.h"
#include <sstream>
#include "math_funcs.h" // create_GARFIELD
#include "IOfunctionsGen.h"
#include "Lag2Eul.h"
#include "init_par.h" // INIT_COSMOLOGY

using namespace std;

void load_arguments(int argc, char *argv[], string &fname_in, unsigned &N1, real_prec &L1, string &fname_out, real_prec &ascale)
{
  stringstream N1_arg, L1_arg, ascale_arg;

  if (argc >= 5)
  {
    fname_in = string(argv[1]);
    N1_arg << argv[2];
    N1_arg >> N1;
    L1_arg << argv[3];
    L1_arg >> L1;
    fname_out = string(argv[4]);
    if (argc >= 6)
    {
      ascale_arg << argv[5];
      ascale_arg >> ascale;
    }
    else
    {
      ascale = 1.;
    }
  }
  else
  {
    cerr << "Need 4 parameters (file in, N1, L1, file out)! N.B.: filenames must be given without extension (which must be .dat). Optional: a_scale (default 1)." << endl;
    exit(1);
  }
}

int main(int argc, char *argv[])
{
  unsigned N1;
  real_prec L1, ascale;
  string fname_in, fname_out;
  load_arguments(argc, argv, fname_in, N1, L1, fname_out, ascale);
 
  ULONG N = N1*N1*N1;

  real_prec d1 = L1/static_cast<real_prec>(N1);
  real_prec min1 = 0.;

  fftw_array<real_prec> delta_lag(N), delta_eul(N), posx(N), posy(N), posz(N);

  // load input file
  get_scalar(fname_in, delta_lag, N1, N1, N1);

  struct COSMOLOGY *c = new COSMOLOGY;
  c->ascale = ascale;
  INIT_COSMOLOGY(c, string("LAG2EULer"));

  // other
  int mk = 1; // CIC -> don't need kernel_scale
  real_prec kth = 4.; // whatever

  // Lag2Eul
  unsigned facL = 1;
  bool reggrid = true;
  gsl_rng *gBaseRand = nullptr; // rng -> not necessary
  real_prec kernel_scale = 0.; // CIC -> don't need kernel_scale

  int sfmodel = 1;
  std::cout << "N.B.: sfmodel is fixed to 1 (Zel'dovich)!" << std::endl;

  Lag2Eul(delta_lag, delta_eul, posx, posy, posz, N1, N1, N1, L1, L1, L1, d1, d1, d1, min1, min1, min1, c->D1, c->D2,
          c->ascale, c->omega_m, c->omega_q, sfmodel, mk, kth, facL,
          reggrid, gBaseRand, "", kernel_scale, nullptr, nullptr);

  quick_dump_scalar(delta_eul, N1, fname_out, 0, false);

  return(0);
}
