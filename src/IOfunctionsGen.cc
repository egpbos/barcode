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

#include <iomanip>
#include <cassert>
#include <sstream>
#include <netinet/in.h>

#include <gsl/gsl_randist.h>

#include "../planck/bstream.h"
#include "../planck/paramfile.h"

// for dump_deltas:
#include "Lag2Eul.h"
#include "HMC_models.h"

#include <chrono> // for timestamp in quick_dump_scalar

using namespace std;

ULONG count_lines(string fname)
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

int get_scalar(string FNAME, real_prec *OUT, unsigned int N1, unsigned int N2, unsigned int N3) {
  ULONG N=N1*N2*N3;	
  fftw_array<real_prec> dummy(N);

  string fname=FNAME+string(".dat");
  // cout<<"... reading file "<<FNAME<<endl;
  bifstream inStream(fname.data());
  assert(inStream.is_open());
  inStream.get(dummy.data,N);
  inStream.close();

  for(ULONG i=0;i<N;i++)
    OUT[i]=dummy[i];

  return 0;
}


void dump_scalar(real_prec *A_rm, unsigned int N1, unsigned int N2, unsigned int N3, string fname) {
  ULONG N=N1*N2*N3;

  fftw_array<real_prec> dummy(N);

  for(ULONG i=0;i<N;i++)
    dummy[i]=A_rm[i];

  string FNAME=fname+string(".dat");
#ifdef DEBUG
  cout<<"... writing file "<<FNAME<<endl;
#endif // DEBUG
  bofstream outStream(FNAME.data());
  assert(outStream.is_open());
  outStream.put(dummy.data,N);
  outStream.close();
}

// EGP: handy function for debugging, only for cubical boxes and only when not using amira mesh!
// Note: N1 is the linear size.
void quick_dump_scalar(real_prec *A_rm, unsigned int N1, string fname, unsigned int sample_number = 0,
                       bool prepend_timestamp = true)
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
void dump_deltas(struct DATA *data, real_prec *deltaLAG, real_prec *deltaS, string fname_append)
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

// EGP: another version, which only takes one length argument
void read_array(string FNAME, real_prec *out, ULONG N)
{
  fftw_array<real_prec> dummy(N);

  string fname=FNAME+string(".dat");
  // cout<<"... reading file "<<FNAME<<endl;
  
  bifstream inStream(fname.data());
  assert(inStream.is_open());
  inStream.get(dummy.data,N);
  inStream.close();

  for(ULONG i=0;i<N;i++)
    out[i]=dummy[i];
}

void read_array(string fname, real_prec *out, unsigned int N1, unsigned int N2, unsigned int N3)
{
  ULONG N=N1*N2*N3;
  read_array(fname, out, N);
}

void write_array(string fname, real_prec *A_rm, ULONG N)
{
  fftw_array<real_prec> dummy(N);

  for(ULONG i=0;i<N;i++)
    dummy[i]=A_rm[i];

  string FNAME=fname+string(".dat");
#ifdef DEBUG
  cout<<"... writing file "<<FNAME<<endl;
#endif // DEBUG
  bofstream outStream(FNAME.data());
  assert(outStream.is_open());
  outStream.put(dummy.data,N);
  outStream.close();
}

void write_array(string fname, real_prec *A_rm, unsigned int N1, unsigned int N2, unsigned int N3)
{
  ULONG N=N1*N2*N3;
  write_array(fname, A_rm, N);
}


void dump_signal_it(ULONG iGibbs, unsigned int N1, unsigned int N2, unsigned int N3, real_prec *signal, string filnam) {
  ULONG N=N1*N2*N3;
  fftw_array<real_prec> dummy(N);

  for(ULONG i=0;i<N;i++)
    dummy[i]=signal[i];

  int bmax=100;
  char buffer1[bmax];  	    
  sprintf(buffer1,"_%lu.dat",iGibbs); // EGP: %d -> %lu

  string FileName=filnam+buffer1;
#ifdef DEBUG
  cout<<"... writing file "<<FileName<<endl;
#endif // DEBUG
  bofstream outStream(FileName.data());
  assert(outStream.is_open());
  outStream.put(dummy.data,N);
  outStream.close();
}
