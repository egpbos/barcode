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

#include "IOfunctionsGen.h"

#include <iomanip>
#include <cassert>
#include <sstream>
#include <netinet/in.h>
#include <fstream>

#include <gsl/gsl_randist.h>

// for dump_deltas:
#include "Lag2Eul.h"
#include "HMC_models.h"

#include <chrono> // for timestamp in quick_dump_scalar

using namespace std;

ULONG count_lines(const string& fname)
{
  ifstream inStream;
  inStream.open(fname.data());
  assert(inStream.is_open());
  string line;
  ULONG linenum = 0;
  while (getline (inStream, line))
  {
    linenum++;     
  }
  inStream.close();

  return(linenum);
}


void change_rm2cm (real_prec *A,float *B,int N1,int N2,int N3)
{
  /*
     Change Row major to column major
     */
  for(int i=0;i<N1;i++)
    for(int j=0;j<N2;j++)
      for(int k=0;k<N3;k++)
      {		
        B[i+N1*(j+N2*k)]=float(A[k+N3*(j+N2*i)]);
        //B[k+N3*(j+N2*i)]=A[k+N3*(j+N2*i)];
      }
}


void calc_BoundingBox(real_prec *BBox, real_prec L1, real_prec L2, real_prec L3)
{
  // the length is calculated from Voxel center to Voxel center
  BBox[0]=0.;
  BBox[1]=L1;

  BBox[2]=0.;
  BBox[3]=L2;

  BBox[4]=0.;
  BBox[5]=L3;
}



void change_cm2rm (float *A,real_prec *B,int N1,int N2,int N3)
{
  /*
     Change Row major to column major
     */
  int i_eff_cm=0;
  int j_eff_cm=0;
  int k_eff_cm=0;

  int i_eff_rm=0;
  int j_eff_rm=0;
  int k_eff_rm=0;

  for(int i=0;i<N1;i++)
    for(int j=0;j<N2;j++)
      for(int k=0;k<N3;k++)
      {
        i_eff_rm=N3*N2*i;
        j_eff_rm=N3*j;
        k_eff_rm=k;

        i_eff_cm=i;
        j_eff_cm=N1*j;
        k_eff_cm=N2*N1*k;

        B[i_eff_rm+j_eff_rm+k_eff_rm]=float (A[i_eff_cm+j_eff_cm+k_eff_cm]);
      }
}



// EGP: handy function for debugging, only for cubical boxes and only when not using amira mesh!
// Note: N1 is the linear size.
void quick_dump_scalar(real_prec *A_rm, unsigned int N1, const string& fname, unsigned int sample_number,
                       bool prepend_timestamp)
{
  stringstream filename;
  if (prepend_timestamp)
  {
    chrono::milliseconds ms = chrono::duration_cast< chrono::milliseconds >(
      chrono::high_resolution_clock::now().time_since_epoch()
      );
    filename << ms.count() << "_";
  }
  filename << fname;
  if (sample_number)
    filename << "_" << sample_number;

  dump_scalar(A_rm, N1, N1, N1, filename.str());
}

// Dump all delta's
void dump_deltas(struct DATA *data, real_prec *deltaLAG, real_prec *deltaS, const string& fname_append)
{
  // TODO: add appendage to signify exact Eulerian model (Zel'dovich, 2LPT, ALPT)

  struct NUMERICAL *n = data->numerical;
  struct COSMOLOGY *c = data->cosmology;

  string fname = n->dir + string("deltaLAG") + fname_append;
  dump_scalar(deltaLAG, n->N1, n->N2, n->N3, fname);

  if (!n->rsd_model)
  {
    // without rsd_model, deltaS (S = q + Psi) is simply in Eulerian space
    fname= n->dir + string("deltaEUL") + fname_append;
    dump_scalar(deltaS, n->N1, n->N2, n->N3, fname);
  } 
  else
  {
    // The deltaS is now the redshift space density field!
    fname = n->dir + string("deltaRSS") + fname_append;
    dump_scalar(deltaS, n->N1, n->N2, n->N3, fname);

    // Also output the real Eulerian density field
    unsigned facL = 1;
    bool reggrid = true;
    gsl_rng *seed = nullptr; // empty: reggrid is true anyway
    real_prec kernel_scale = SPH_kernel_scale<struct DATA>(data);
    fftw_array<real_prec> deltaEUL(n->N), posx(n->N), posy(n->N), posz(n->N);
    Lag2Eul(deltaLAG, deltaEUL, posx, posy, posz, n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1, n->d2, n->d3,
            n->xllc, n->yllc, n->zllc, c->D1, c->D2, c->ascale, c->omega_m, c->omega_q, n->sfmodel, n->mk, n->slength, facL,
            reggrid, seed, "", kernel_scale, n->R2Cplan, n->C2Rplan);

    fname = n->dir + string("deltaEUL") + fname_append;
    dump_scalar(deltaEUL, n->N1, n->N2, n->N3, fname);
  }
}

int get_scalar(const string& FNAME, real_prec *OUT, unsigned int N1, unsigned int N2, unsigned int N3) {
  read_array(FNAME, OUT, N1, N2, N3);
  return 0;
}

void read_array(const string& fname, real_prec *out, unsigned int N1, unsigned int N2, unsigned int N3)
{
  ULONG N=N1*N2*N3;
  read_array(fname, out, N);
}


std::string add_extension_if_missing(const std::string& fn, std::string ext = ".dat") {
  if (fn.rfind('.') == std::string::npos) {
    return fn + ext;
  } else {
    return fn;
  }
}


void read_array(const string& FNAME, real_prec *out, ULONG N) {
  string fname = add_extension_if_missing(FNAME);

  std::ifstream inStream(fname.c_str(), std::ios::binary);
  if(inStream.is_open()) {
    inStream.read(reinterpret_cast<char *>(out), static_cast<std::streamsize>(N * sizeof(real_prec)));
  } else {
    throw BarcodeException(std::string("In read_array: error opening file ") + fname);
  }
}


void dump_scalar(real_prec *A_rm, unsigned int N1, unsigned int N2, unsigned int N3, const string& fname) {
  write_array(fname, A_rm, N1, N2, N3);
}

void write_array(const string& fname, real_prec *A_rm, unsigned int N1, unsigned int N2, unsigned int N3)
{
  ULONG N=N1*N2*N3;
  write_array(fname, A_rm, N);
}

void write_array(const string& fname, real_prec *A_rm, ULONG N) {
  string FNAME = add_extension_if_missing(fname);

#ifdef DEBUG
  cout<<"... writing file "<<FNAME<<endl;
#endif // DEBUG

  std::ofstream outStream(FNAME.c_str(), std::ios::binary);
  if (outStream.is_open()) {
    outStream.write(reinterpret_cast<const char *>(A_rm), static_cast<std::streamsize>(N * sizeof(real_prec)));
  } else {
    throw BarcodeException(std::string("In write_array: error opening file ") + FNAME);
  }
}


void dump_signal_it(ULONG iGibbs, unsigned int N1, unsigned int N2, unsigned int N3, real_prec *signal, const string& filnam) {
  ULONG N=N1*N2*N3;

  int bmax=100;
  char buffer1[bmax];
  sprintf(buffer1,"_%lu.dat", iGibbs); // EGP: %d -> %lu

  string FileName=filnam+buffer1;
  write_array(FileName, signal, N);
}
