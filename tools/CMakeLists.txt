#
# Barcode
# Copyright E.G.P. Bos and F.S. Kitaura
#
# Distributed under the terms of the MIT License.
# The full license is in the file LICENSE, distributed with this software.
#

include(AddDefaultBinary)

## density
ADD_DEFAULT_BINARY(density density.cc)

## convert primordial density fluctuation field ("lagrangian") to z=0 density field ("eulerian")
ADD_DEFAULT_BINARY(LAG2EULer LAG2EULer.cc)

## correlation function calculator
ADD_DEFAULT_BINARY(corr_fct corr_fct.cc)

## powspec
ADD_DEFAULT_BINARY(powspec powspec.cc)

## 2D correlation function generators (with or without interpolating resolution upscaling)
ADD_DEFAULT_BINARY(2D_corr_fct 2D_corr_fct.cc)
ADD_DEFAULT_BINARY(2D_corr_fct_interp 2D_corr_fct_interp.cc)
target_link_libraries(2D_corr_fct_interp PRIVATE ${FFTW_DOUBLE_LIB})
if (MULTITHREAD_FFTW)
    target_link_libraries(2D_corr_fct_interp PRIVATE ${FFTW_DOUBLE_OPENMP_LIB})
endif()

## 2D power spectrum calculator
ADD_DEFAULT_BINARY(2D_powspec 2D_powspec.cc)

## tools for upscaling resolution (with poisson sampling or interpolation)
ADD_DEFAULT_BINARY(interp_upres interp_upres.cc)
ADD_DEFAULT_BINARY(poisson_upres poisson_upres.cc)
