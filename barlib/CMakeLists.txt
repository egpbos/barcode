#
# Barcode
# Copyright E.G.P. Bos and F.S. Kitaura
#
# Distributed under the terms of the MIT License.
# The full license is in the file LICENSE, distributed with this software.
#

add_library(barlib
        # CORE
        src/barcoderunner.cc
        src/protocol.cc         # initialization of ... not sure what exactly. TODO: REMOVE, ONLY ONE THAT'S USED IS RESTART_PROTOCOL
        src/sample_maker.cc
        # CONFIG FILE
        src/ini_reader.cpp
        src/init_par.cc
        # HMC
        src/call_hamil.cc
        src/HMC.cc
        # split off of old HMC.cc:
        src/hmc/leapfrog/time_step.cpp
        # --- end split off HMC.cc
        src/HMC_help.cc
        src/HMC_mass.cc
        src/HMC_momenta.cc
        src/HMC_models.cc
        # split off of old HMC_models.cc:
        src/hmc/prior/gaussian.cpp
        src/hmc/likelihood/gaussian_random_field.cpp
        src/hmc/likelihood/poissonian.cpp
        src/hmc/likelihood/gaussian_independent.cpp
        src/hmc/likelihood/lognormal_independent.cpp
        src/HMC_models_testing.cpp
        src/SPH_kernel.cpp
        # --- end split off HMC_models.cc
        # ARRAY, MATH AND STRING MANIPULATIONS:
        src/convenience.cc
        src/math_funcs.cc
        # split off of old math_funcs.cc:
        src/random.cpp  # included includes in random.hpp
        src/interpolate_grid.cpp
        src/pacman.cpp
        src/scale_space.cpp
        src/convolution.cpp
        src/gradient.cpp
        # COSMOLOGY
        src/cosmo.cc
        src/calc_power.cc       # calculate theoretical cosmological power spectra
        # INTERFACE
        src/curses_funcs.cc
        # DEBUGGING
        src/debug.cc
        # LAG2EUL
        src/Lag2Eul.cc
        src/disp_part.cc
        src/EqSolvers.cc
        # FFT
        src/fftwrapper.cc
        # IO
        src/IOfunctions.cc
        src/IOfunctionsGen.cc
        # DENSITY
        src/massFunctions.cc
        src/transf.cpp
        src/rankorder.cc         # aren't these also just transfer functions? TODO: maybe merge with transf.cpp
        src/field_statistics.cpp # field statistics; power spectrum, etc
        # REDSHIFT SPACE
        src/rsd.cc)

target_compile_features(barlib PUBLIC cxx_std_11)
TARGET_ADD_OPTIONS(barlib)

target_include_directories(barlib PUBLIC include ${GSL_INCLUDE_DIRS} ${FFTW_INCLUDE_DIRS})

target_link_libraries(barlib PRIVATE ${GSL_LIBRARIES} m ncurses)
TARGET_ADD_FFTW(barlib)
if (MULTITHREAD)
    target_link_libraries(barlib PRIVATE OpenMP::OpenMP_CXX)
endif()

#--------------------------------------- installation
install(TARGETS barlib
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
