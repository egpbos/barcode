/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include <string>  // string
#include <sstream>  // stringstream

#include "../define_opt.h"  // real_prec, ULONG
#include "../fftw_array.h"  // fftw_array
#include "../IOfunctionsGen.h"  // get_scalar, dump_scalar
#include "../IOfunctions.h"  // dump_measured_spec
#include "../field_statistics.h" // measure_spectrum

using std::cout;
using std::cerr;
using std::endl;
using std::string;

void load_arguments(int argc, char *argv[], string &fname_in, unsigned &N1,
                    real_prec &L1, ULONG &N_bin, string &fname_out) {
  stringstream N1_arg, L1_arg, N_bin_arg;

  if (argc >= 5) {
    fname_in = string(argv[1]);
    N1_arg << argv[2];
    N1_arg >> N1;
    L1_arg << argv[3];
    L1_arg >> L1;
    N_bin_arg << argv[4];
    N_bin_arg >> N_bin;

    if (argc >= 6)
      fname_out = string(argv[5]);
    else
      fname_out = fname_in + string("_pow");
  } else {
    cerr << "Need 4 parameters (file in, N1, L1, N_bin)!" << endl
         << "N.B.: filenames must be given without extension (which must be"
            " .dat)." << endl
         << "Optional: fname_out (default: fname_in+'_pow')." << endl;
    exit(1);
  }
}

int main(int argc, char *argv[]) {
  unsigned N1;
  ULONG N_bin;
  real_prec L1;
  string fname_in, fname_out;
  load_arguments(argc, argv, fname_in, N1, L1, N_bin, fname_out);

  ULONG N = N1*N1*N1;

  fftw_array<real_prec> power(N_bin), kmode(N_bin);
  fftw_array<real_prec> input(N);

  // get data from input file
  get_scalar(fname_in, input, N1, N1, N1);

  // do measurement
  measure_spectrum(N1, N1, N1, L1, L1, L1, input, kmode, power, N_bin);

  // output results
  dump_measured_spec(kmode, power, fname_out, N_bin);

  return(0);
}
