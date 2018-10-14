/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#pragma once
#include <string>
#include "define_opt.h"
#include "struct_main.h"

using namespace std;

ULONG count_lines(string fname);

void change_rm2cm (real_prec *A,float *B,int N1,int N2,int N3);

void change_cm2rm (float *A,real_prec *B,int N1,int N2,int N3);

void calc_BoundingBox(real_prec *BBox, real_prec L1, real_prec L2, real_prec L3);

void dump_scalar(real_prec *A_rm, unsigned int N1, unsigned int N2, unsigned int N3, string fname);
void quick_dump_scalar(real_prec *A_rm, unsigned int N1, string fname, unsigned int sample_number = 0,
                       bool prepend_timestamp = true);
void dump_deltas(struct DATA *data, real_prec *deltaLAG, real_prec *deltaS, string fname_append);

int get_scalar(string FNAME, real_prec *OUT, unsigned int N1, unsigned int N2, unsigned int N3);

int read_amiramesh(char* filename, float* data, unsigned int* dims, float* bnds, float& time, bool readdata);

void read_array(string FNAME, real_prec *out, ULONG N);
void read_array(string fname, real_prec *out, unsigned int N1, unsigned int N2, unsigned int N3);

void write_array(string fname, real_prec *A_rm, ULONG N);
void write_array(string fname, real_prec *A_rm, unsigned int N1, unsigned int N2, unsigned int N3);

void dump_signal_it(ULONG iGibbs, unsigned int N1, unsigned int N2, unsigned int N3, real_prec *signal, string filnam);
