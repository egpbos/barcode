/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include "struct_main.h"

#include <iomanip>
#include <cassert>

using namespace std;

void dump_measured_spec(real_prec *x, real_prec *y, const std::string& outputFileName,
  ULONG N_bin) {
  ofstream outStream(outputFileName.data());
  assert(outStream.is_open());

#ifdef DEBUG
  cout << " >>> dump power-spectrum in : " << outputFileName << " ...."
    << endl << endl;
#endif  // DEBUG

  for (ULONG i = 0; i < N_bin; i++)
    if (y[i] > 0. && x[i] > 0.)
     outStream << x[i] << "   " << y[i] << endl;

  outStream.close();
}


void dump_ps_it(struct DATA *data,real_prec *x, const real_prec *y)
{
  //int N1=data->numerical->N1;
  //int N2=data->numerical->N2;
  //int N3=data->numerical->N3;

  //ULONG N=N1*N2*N3;

  //real_prec L1=data->numerical->L1;
  //real_prec L2=data->numerical->L2;
  //real_prec L3=data->numerical->L3;
  
  ULONG N_bin=data->numerical->N_bin;
  
  int bmax=100;
  char buffer1[bmax];  	    
  sprintf(buffer1,"it%lu.dat",data->numerical->iGibbs); // EGP: %d -> %lu
	    
  string outputFileName= data->numerical->dir + string("powSpec")+buffer1;	
  ofstream outStream(outputFileName.data());
  
  assert(outStream.is_open());
    
  auto NORM=static_cast<real_prec>(1.0); // output normalisation
  // EGP: this is now done in measure_spec already (EqSolvers.cc)

//#ifdef	  FOURIER_DEF_1 
  //NORM=static_cast<real_prec>(L1*L2*L3/4./M_PI); // output normalisation
//#endif  
//#ifdef	  FOURIER_DEF_2 
  //NORM=static_cast<real_prec>(L1*L2*L3/real_prec(N)/real_prec(N)/4./M_PI); // output normalisation
//#endif

#ifdef DEBUG
  cout << " >>> dump power-spectrum in : "<<outputFileName <<"...."<<endl<<endl;
#endif // DEBUG
  
  for(ULONG i=0;i<N_bin;i++)
    {
      if(y[i]>0. && x[i]>0.)
	{
	  outStream << x[i] <<"   " << NORM*y[i]<<endl;
	  
	}
    }
  outStream.close();
}
