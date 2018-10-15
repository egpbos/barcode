/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include <cmath>
#include <string>
#include <sstream>
#include <vector>

#include "define_opt.h"
#include "struct_main.h"
#include "fftw_array.h"
#include "IOfunctionsGen.h"
#include "math_funcs.h"  // interpolate_CIC

// using namespace std;  // forbidden by Google C++ Style Guide
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;

void load_arguments(int argc, char *argv[], string &fname_in, unsigned &N1,
                    real_prec &L1, unsigned &N1_out, string &fname_out) {
  stringstream N1_arg, L1_arg, N1_out_arg, Nbar_arg, seed_arg;

  if (argc >= 5) {
    fname_in = string(argv[1]);
    N1_arg << argv[2];
    N1_arg >> N1;
    L1_arg << argv[3];
    L1_arg >> L1;
    N1_out_arg << argv[4];
    N1_out_arg >> N1_out;

    if (argc >= 6) {
      fname_out = string(argv[5]);
    } else {
      fname_out = fname_in + string("_interpCIC") + N1_out_arg.str();
    }

    cout << "will output to file: " << fname_out << endl;
  } else {
    cerr << "Need 4 parameters (file in, N1, L1, N1_out)! " << endl
         << "N.B.: filenames must be given without extension (which must be "
            ".dat)." << endl
         << "Optional: fname_out (default: fname_in+'_interpCIC[N1_out]')."
         << endl;
    exit(1);
  }
}

void interp_field(real_prec *input, unsigned N1, unsigned N2, unsigned N3, real_prec L1,
                  real_prec L2, real_prec L3, unsigned N1_out, unsigned N2_out,
                  unsigned N3_out, real_prec *output) {
  real_prec dx = L1 / static_cast<real_prec>(N1);
  real_prec dy = L2 / static_cast<real_prec>(N2);
  real_prec dz = L3 / static_cast<real_prec>(N3);
  real_prec dx_out = L1 / static_cast<real_prec>(N1_out);
  real_prec dy_out = L2 / static_cast<real_prec>(N2_out);
  real_prec dz_out = L3 / static_cast<real_prec>(N3_out);

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif  // MULTITHREAD
  for (unsigned i = 0; i < N1_out; ++i)
    for (unsigned j = 0; j < N2_out; ++j)
      for (unsigned k = 0; k < N3_out; ++k) {
        ULONG ix = k + N3_out*(j+N2_out*i);

        // +0.5, assuming we associate grid values with values at center of grid
        // cells
        real_prec posx = dx_out * (0.5 + static_cast<real_prec>(i));
        real_prec posy = dy_out * (0.5 + static_cast<real_prec>(j));
        real_prec posz = dz_out * (0.5 + static_cast<real_prec>(k));

        output[ix] = interpolate_CIC(N1, N2, N3, L1, L2, L3, dx, dy, dz,
                                     posx, posy, posz, input);
      }
}


int main(int argc, char *argv[]) {
  unsigned N1, N1_out;
  real_prec L1;
  string fname_in, fname_out;
  load_arguments(argc, argv, fname_in, N1, L1, N1_out, fname_out);

  ULONG N = N1*N1*N1;
  ULONG N_out = N1_out*N1_out*N1_out;

  fftw_array<real_prec> input(N);

  // get data from input file
  get_scalar(fname_in, input, N1, N1, N1);

  // interpolate
  cout << "... interpolate" << endl;
  fftw_array<real_prec> result(N_out);
  interp_field(input, N1, N1, N1, L1, L1, L1, N1_out, N1_out, N1_out, result);

  // output results
  quick_dump_scalar(result, N1_out, fname_out, 0, false);

  return(0);
}

