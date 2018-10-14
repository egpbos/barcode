/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include <math.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cassert>

#include "planck/paramfile.h"

#include "struct_main.h"

using namespace std;

void INIT_PROTOCOL_CONV(struct DATA *data)
{
  string outputFileName= data->numerical->dir + string("convergence.prt");
  
  ofstream outStream(outputFileName.data());
  assert(outStream.is_open());
  
  outStream.close();
}

void UPDATE_PROTOCOL_CONV(ULONG it, real_prec res, struct DATA *data)
{
  string outputFileName= data->numerical->dir + string("convergence.prt");
  
  ofstream outStream(outputFileName.data(), ios::app );
  assert(outStream.is_open());
  
  outStream << it <<" "<< res  <<endl;
  outStream.close();		
}

void PROTOCOL_RESTART(ULONG it, struct DATA *data)
{
  string outputFileName= data->numerical->dir + string("restart.prt");
		
  ofstream outStream(outputFileName.data());
  assert(outStream.is_open());

  outStream << it <<endl;
  outStream.close();	       
}

void INIT_PROTOCOL_SPEC(struct DATA *data)
{
  unsigned N1=data->numerical->N1;
  unsigned N2=data->numerical->N2;
  unsigned N3=data->numerical->N3;
  
  real_prec L1=data->numerical->L1;
  real_prec L2=data->numerical->L2;
  real_prec L3=data->numerical->L3;

  ULONG N_bin=data->numerical->N_bin;

  string outputFileName= data->numerical->dir + string("spec_protocol.prt");
  
  ofstream outStream(outputFileName.data());
  assert(outStream.is_open());

  string codename=data->numerical->codename;
  
  // cout << " >>> "<<codename<<": initialising power-spectrum protocol...."<<endl<<endl;
  // cout<<"\n >>> "<<codename<<": create power-spectrum protocol file: "<<outputFileName<<" ...\n";
  wprintw(data->curses->status, "initialising power-spectrum protocol...");
  wrefresh(data->curses->status);
  
  outStream <<" "<<codename<<"_POWER_SPECTRUM_PROTOCOL "<<endl;
  outStream << "-----------------------------"<<endl;
  outStream << N1  <<endl;
  outStream << N2  <<endl;
  outStream << N3  <<endl;	
  outStream << "-----------------------------"<<endl;
  outStream << L1  <<endl;
  outStream << L2  <<endl;
  outStream << L3  <<endl;	
  outStream << "-----------------------------"<<endl;
  outStream << N_bin  <<endl;
  outStream << "-----------------------------"<<endl;
  outStream.close();  
}

void UPDATE_PROTOCOL_SPEC(struct DATA *data, string specname)
{
//EGP  int N1=data->numerical->N1;
//EGP  int N2=data->numerical->N2;
//EGP  int N3=data->numerical->N3;
  
//EGP  real_prec L1=data->numerical->L1;
//EGP  real_prec L2=data->numerical->L2;
//EGP  real_prec L3=data->numerical->L3;
  
//EGP  ULONG N_bin=data->numerical->N_bin;
  
  string outputFileName= data->numerical->dir + string("spec_protocol.prt");
  
  ofstream outStream(outputFileName.data(), ios::app );
  assert(outStream.is_open());
  
  string codename=data->numerical->codename;

  cout << " >>> "<<codename<<": update power-spectrum protocol...."<<endl<<endl;
  outStream << specname  <<endl;
  outStream.close();  
}
