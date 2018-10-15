/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#include "struct_main.h"

#include <gsl/gsl_randist.h>

#include "call_hamil.h"


void security_recursion(struct DATA *data, gsl_rng * seed)
{
  if(data->numerical->INV_SUCCESS==0)
    {
      call_hamil(data, seed);
      
      security_recursion(data,seed);
    }
}


void sample_maker(struct DATA *data, gsl_rng * seed)
{    
  if(data->numerical->INV_SUCCESS==0)    
    security_recursion(data, seed);
  
  data->numerical->INV_SUCCESS=0;	     
}

