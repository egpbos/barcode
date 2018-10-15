/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include <string>  // string
#include <sstream>  // stringstream

#include "define_opt.h"  // real_prec, ULONG
#include "fftw_array.h"  // fftw_array
#include "IOfunctionsGen.h"  // get_scalar, dump_scalar
#include "massFunctions.h"  // getDensity_SPH
#include "convenience.h" // fill_one

using std::cout;
using std::cerr;
using std::endl;
using std::string;

void load_arguments(int argc, char *argv[], string &fname_inx,
                    string &fname_iny, string &fname_inz, unsigned &N1,
                    real_prec &L1, ULONG &N_part, string &fname_out) {
  stringstream N1_arg, L1_arg, N_part_arg;

  int N_required = 6;

  if (argc - 1 >= N_required) {
    fname_inx = string(argv[1]);
    fname_iny = string(argv[2]);
    fname_inz = string(argv[3]);
    N1_arg << argv[4];
    N1_arg >> N1;
    L1_arg << argv[5];
    L1_arg >> L1;
    N_part_arg << argv[6];
    N_part_arg >> N_part;

    if (argc - 1 >= N_required + 1)
      fname_out = string(argv[7]);
    else
      fname_out = string("density");
  } else {
    cerr << "Need " << N_required << " parameters (input files x, y, z, N1, L1, N_part)!" << endl
         << "N.B.: filenames must be given without extension (which must be"
            " .dat)." << endl
         << "Optional: fname_out (default: 'density')." << endl;
    exit(1);
  }
}

int main(int argc, char *argv[]) {
  unsigned N1;
  ULONG N_part;
  real_prec L1;
  string fn_inx, fn_iny, fn_inz, fn_out;
  load_arguments(argc, argv, fn_inx, fn_iny, fn_inz, N1, L1, N_part, fn_out);

  real_prec d1 = L1/static_cast<real_prec>(N1);

  ULONG N = N1*N1*N1;

  fftw_array<real_prec> density(N), x(N_part), y(N_part), z(N_part);

  // get data from input file
  read_array(fn_inx, x, N_part);
  read_array(fn_iny, y, N_part);
  read_array(fn_inz, z, N_part);

  // do density estimation
  real_prec min1 = 0.;
  fftw_array<real_prec> particle_mass(N_part);
  fill_one(particle_mass, N_part);
  bool weightmass = true;
  real_prec kernel_scale = d1;  // default in barcode

  getDensity_SPH(N1, N1, N1, L1, L1, L1, d1, d1, d1, min1, min1, min1, x, y, z,
                 particle_mass, N_part, density, weightmass, kernel_scale);

  // output results
  quick_dump_scalar(density, N1, fn_out, 0, false);

  return(0);
}
