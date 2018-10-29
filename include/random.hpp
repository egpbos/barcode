/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


#pragma once
#include <gsl/gsl_randist.h> // includes gsl_rng.h
#include <vector>
#include <complex>

template <class out_type>
std::complex<out_type> complex_gaussian_random_number(gsl_rng *rng)
{
  // _ugaussian is the same as _gaussian, but with sigma = 1.0
  auto val_real = static_cast<out_type>(gsl_ran_ugaussian(rng)); // real part
  auto val_imag = static_cast<out_type>(gsl_ran_ugaussian(rng)); // imaginary part
  std::complex<out_type> val(val_real, val_imag);
  return(val);
}

// Can only do a "cubical" grid, i.e. there is only one linear gridsize.
// The grid is generated for use in Fourier space. Its size is determined by the
// half_size parameter; if false, it will give a full size cube, otherwise a
// "half + 1" sized grid for use in a real-FFT. Note that this function returns
// a complex vector.
template <class out_type>
std::vector<std::complex<out_type> > resolution_independent_random_grid_FS(unsigned gridsize, gsl_rng *rng, bool half_size = true)
{
  unsigned i, j, k; // just some for-loop indices

  // Build output vector
  ULONG output_size;
  if (half_size)
    output_size = gridsize * gridsize * (gridsize/2+1);
  else
    output_size = gridsize * gridsize * gridsize;
  std::vector<std::complex<out_type> > out_grid(output_size);

  // The strides are also dependent on grid size:
  unsigned j_stride;
  if (half_size)
    j_stride = (gridsize/2+1);
  else
    j_stride = gridsize;
  unsigned i_stride = gridsize * j_stride;
  // Note: i_stride does not necessarily go with i in this function, but
  // rather with the i that one would have in an outer for-loop with index
  // i of the form:
  // for i=[0-N) for j=[0-N) for k=[0-N) (or k=[0-(N/2+1)))
  // Same goes for the j_stride. k_stride == 1 in this way.

  // N.B.: don't make the for loops multi-threaded!
  
  // Fill out:
  for (i = 0; i < gridsize/2; i++) // "layer" count; start in the corners, and move out one cubical layer at a time 
  {
    for (k = 0; k < i+1; k++) // first do the two "walls"
    {
      for (j = 0; j < i; j++) // "slim" side
      {
        out_grid[i * i_stride + j * j_stride + k] = complex_gaussian_random_number<out_type>(rng); // corner 1
        out_grid[(gridsize - 1 - i) * i_stride + j * j_stride + k] = complex_gaussian_random_number<out_type>(rng); // corner 2
        out_grid[i * i_stride + (gridsize - 1 - j) * j_stride + k] = complex_gaussian_random_number<out_type>(rng); // corner 3
        out_grid[(gridsize - 1 - i) * i_stride + (gridsize - 1 - j) * j_stride + k] = complex_gaussian_random_number<out_type>(rng); // corner 4
        if (not half_size)
        {
          out_grid[i * i_stride + j * j_stride + (gridsize - 1 - k)] = complex_gaussian_random_number<out_type>(rng); // corner 5
          out_grid[(gridsize - 1 - i) * i_stride + j * j_stride + (gridsize - 1 - k)] = complex_gaussian_random_number<out_type>(rng); // corner 6
          out_grid[i * i_stride + (gridsize - 1 - j) * j_stride + (gridsize - 1 - k)] = complex_gaussian_random_number<out_type>(rng); // corner 7
          out_grid[(gridsize - 1 - i) * i_stride + (gridsize - 1 - j) * j_stride + (gridsize - 1 - k)] = complex_gaussian_random_number<out_type>(rng); // corner 8
        }
      }
      for (j = 0; j < i + 1; j++) // "broad" side
      {
        out_grid[j * i_stride + i * j_stride + k] = complex_gaussian_random_number<out_type>(rng); // corner 1
        out_grid[(gridsize - 1 - j) * i_stride + i * j_stride + k] = complex_gaussian_random_number<out_type>(rng); // corner 2
        out_grid[j * i_stride + (gridsize - 1 - i) * j_stride + k] = complex_gaussian_random_number<out_type>(rng); // corner 3
        out_grid[(gridsize - 1 - j) * i_stride + (gridsize - 1 - i) * j_stride + k] = complex_gaussian_random_number<out_type>(rng); // corner 4
        if (not half_size)
        {
          out_grid[j * i_stride + i * j_stride + (gridsize - 1 - k)] = complex_gaussian_random_number<out_type>(rng); // corner 5
          out_grid[(gridsize - 1 - j) * i_stride + i * j_stride + (gridsize - 1 - k)] = complex_gaussian_random_number<out_type>(rng); // corner 6
          out_grid[j * i_stride + (gridsize - 1 - i) * j_stride + (gridsize - 1 - k)] = complex_gaussian_random_number<out_type>(rng); // corner 7
          out_grid[(gridsize - 1 - j) * i_stride + (gridsize - 1 - i) * j_stride + (gridsize - 1 - k)] = complex_gaussian_random_number<out_type>(rng); // corner 8
        }
      }
    }
    for (j = 0; j < i; j++) // then the "roof" (what remains on the top)
      for (k = 0; k < i; k++)
      {
        out_grid[j * i_stride + k * j_stride + i] = complex_gaussian_random_number<out_type>(rng); // corner 1
        out_grid[(gridsize - 1 - j) * i_stride + k * j_stride + i] = complex_gaussian_random_number<out_type>(rng); // corner 2
        out_grid[j * i_stride + (gridsize - 1 - k) * j_stride + i] = complex_gaussian_random_number<out_type>(rng); // corner 3
        out_grid[(gridsize - 1 - j) * i_stride + (gridsize - 1 - k) * j_stride + i] = complex_gaussian_random_number<out_type>(rng); // corner 4
        if (not half_size)
        {
          out_grid[j * i_stride + k * j_stride + (gridsize - 1 - i)] = complex_gaussian_random_number<out_type>(rng); // corner 5
          out_grid[(gridsize - 1 - j) * i_stride + k * j_stride + (gridsize - 1 - i)] = complex_gaussian_random_number<out_type>(rng); // corner 6
          out_grid[j * i_stride + (gridsize - 1 - k) * j_stride + (gridsize - 1 - i)] = complex_gaussian_random_number<out_type>(rng); // corner 7
          out_grid[(gridsize - 1 - j) * i_stride + (gridsize - 1 - k) * j_stride + (gridsize - 1 - i)] = complex_gaussian_random_number<out_type>(rng); // corner 8
        }
      }
  }

  if (half_size) // Finally, fill out the nyquist plane:
    for (i = 0; i < gridsize; i++) 
      for (j = 0; j < gridsize; j++) 
        out_grid[i_stride*i + j_stride*j + gridsize/2] = complex_gaussian_random_number<out_type>(rng);

  return out_grid;
}
