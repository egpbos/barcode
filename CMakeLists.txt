#
# Barcode
# Copyright E.G.P. Bos and F.S. Kitaura
#
# Distributed under the terms of the MIT License.
# The full license is in the file LICENSE, distributed with this software.
#

cmake_minimum_required(VERSION 3.9.0)
# 3.8 for target_compile_features(target PRIVATE/PUBLIC cxx_std_11/17), 3.9.0 for FindOpenMP, all other features are in 3.0 at least (didn't check 2.x)


#--------------------------------------- download and setup FindFFTW
configure_file(cmake/downloadFindFFTW.cmake.in findFFTW-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
if(result)
    message(FATAL_ERROR "CMake step for findFFTW failed: ${result}")
    else()
    message("CMake step for findFFTW completed (${result}).")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
if(result)
    message(FATAL_ERROR "Build step for findFFTW failed: ${result}")
endif()

set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)

#--------------------------------------- CMake modules

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/" "${findFFTW_DIR}")

project(barcode)

#--------------------------------------- warnings (for all binaries)
add_compile_options(-Wundef -Wshadow -Wwrite-strings -Wredundant-decls -Woverloaded-virtual -Wcast-qual -Wcast-align -Wpointer-arith -Wconversion -Wold-style-cast -Wall -Wextra)


#--------------------------------------- OPTIONS

option(INTEL "Use intel library and compiler" OFF)

#--------------------------------------- Basic operation mode of code
option(DEBUG "Debugging mode" OFF)
option(MULTITHREAD "Enables OpenMP directives on most for-loops, except once with RNGs." ON)
option(MULTITHREAD_FFTW "Use multi-threaded version of the FFTW3 library." ON)
option(MULTITHREAD_RNG "Enable OpenMP directives on for-loops containing random number generators. This causes an extra degree of randomness, due to the order in which cores access the RNG, which will differ every run, because the scheduler controls this and other running processes will have an (unpredictable) influence on this." OFF)
option(NAN_DETECTION "Turns on not-a-number detection. This causes floating point exceptions if a NaN is created somewhere in the code. If you don't do this, the NaNs may freely propagate and your code will keep running, producing nonsense results." ON)
option(RESTART_FILE "Using this option, if you want to restart your simulation, you have to create an empty file called \"restart\" in the output directory. This will then read in the restart.prt file that is created every step. If you don't use this option, you can restart at a step of your choice by simply calling the code with an extra command line option specifying the number you want to start from, e.g. \"barcode.x 203\". The disadvantage of using RESTART_FILE is that it can cause big I/O lags on computers with network storage if you're running multiple instances on the same computer at once." OFF)

#--------------------------------------- Hamiltonian stuff
option(MASKING "[NOT TESTED!] Turns handling of masking on." OFF)

#--------------------------------------- Modelled physics
# Note that these options should be replaced by a sfmodel parameter in
# the input.par file, like in the patchy code! This is work in progress.
option(TRANSF "[NOT TESTED!] Activate transfer function convolution for Zel'dovich and 2LPT models." OFF)
option(TRANSFSC "[NOT TESTED!] Convolve the SC part of Psi with a transfer function as well. This seems to be the same transfer function as the one used for the Zel'dovich+TF model. This option gives a \"ALPT+TF\" model." OFF)
option(BARYONS "Power spectrum with baryons (Eisenstein & Hu 1998). When set to OFF, the power spectrum will be a power law, incl. apodisation at Nyquist scale. Both options only apply when not reading the power spectrum from a file, which can be set in the input file." ON)

#--------------------------------------- Rank ordering
option(RANKORD "NOT IMPLEMENTED. Switch on rank ordering." OFF)
option(RANDORDTYPE "NOT IMPLEMENTED. Do rank ordering split up according to t-web type. If not, just order by density, regardless of type." OFF)
option(LOGRO "NOT USED. Do a logarithmic rank ordering, i.e. transform density to log-density before ordering." OFF)

#--------------------------------------- Fourier transform convention
option(FOURIER_DEF_1 "Defines the type of (discrete) Fourier transforms. With DEF_1, the FT itself has the 1/N factor in front. With DEF_2, the inverse FT has the 1/N factor. Default: FOURIER_DEF_2. Note: code was only tested with FOURIER_DEF_2!" OFF)
option(FOURIER_DEF_2 "Defines the type of (discrete) Fourier transforms. With DEF_1, the FT itself has the 1/N factor in front. With DEF_2, the inverse FT has the 1/N factor. Default: FOURIER_DEF_2. Note: code was only tested with FOURIER_DEF_2!" ON)
option(FOURIER_DEF_2_20151021 "[NOT TESTED!] Experimental update to DEF_2. Can only be used simultaneously with FOURIER_DEF_2, not with DEF_1 and not by itself." OFF)

if ( FOURIER_DEF_1 AND FOURIER_DEF_2 )
    message(FATAL_ERROR "Cannot use FOURIER_DEF_1 and FOURIER_DEF_2 at same time!")
endif()

if ( ( FOURIER_DEF_1 AND FOURIER_DEF_2_20151021 ) OR (FOURIER_DEF_2_20151021 AND NOT FOURIER_DEF_2))
    message(FATAL_ERROR "Cannot use FOURIER_DEF_2_20151021 with FOURIER_DEF_1 or by itself! Only together with DEF_2.")
endif()

#--------------------------------------- Single/Double Precision
option(SINGLE_PREC "Choose single or double precision calculations and input/output." OFF)
option(DOUBLE_PREC "Choose single or double precision calculations and input/output." ON)

if ( SINGLE_PREC AND DOUBLE_PREC )
    message(FATAL_ERROR "Cannot use SINGLE_PREC and DOUBLE_PREC at same time!")
endif()

#--------------------------------------- Testing and Debugging options
option(PLANCK_SILENT "Makes the Planck LevelS package functions silent (used for parsing configuration input file). Turn off to make them more verbose." ON)

#--------------------------------------- Things for special behaviour
option(SCSMOO "[NOT TESTED!] Smooth the spherical collapse part of the displacement field with a smoothing length kthsc2. This variable is currently not implemented in the code. kthsc is implemented and is equal to kth, which equals ksmooth from input.par. Something similar should probably be defined for this." OFF)

#--------------------------------------- NEEDS ATTENTION

# CONSIDER EITHER TAKING THESE OPTIONS OUT OR ACTUALLY GIVING THEM CODE-WIDE IMPACT:
option(GFFT "Use FFT based derivative (NOT REALLY USED, currently only used in non-Zel'dovich Lag2Eul via calc_m2v_mem and in EigenValuesTweb)." OFF)
option(GFINDIFF "Use finite difference derivative (NOT REALLY USED, currently only used in non-Zel'dovich Lag2Eul via calc_m2v_mem and in EigenValuesTweb)." ON)

if ( GFFT AND GFINDIFF )
    message(FATAL_ERROR "Cannot use GFFT and GFINDIFF at same time!")
endif()

if ( NOT (GFINDIFF OR GFFT) )
    message(FATAL_ERROR "Must choose one of GFFT and GFINDIFF!")
endif()

set(BARCODE_SOURCE_FILES
# CORE
        src/barcoderunner.cc
        src/protocol.cc         # initialization of ... not sure what exactly. TODO: REMOVE, ONLY ONE THAT'S USED IS RESTART_PROTOCOL
        src/sample_maker.cc
# CONFIG FILE
        planck/cxxutils.cc  # PARSING
        src/init_par.cc
# HMC
        src/call_hamil.cc
        src/HMC.cc
        src/HMC_help.cc
        src/HMC_mass.cc
        src/HMC_models.cc
        src/HMC_momenta.cc
# ARRAY, MATH AND STRING MANIPULATIONS:
        src/convenience.cc
        src/math_funcs.cc
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

include_directories(include)

add_library(barlib ${BARCODE_SOURCE_FILES})
target_compile_features(barlib PUBLIC cxx_std_11)

add_executable(barcode src/main.cc)
target_link_libraries(barcode PRIVATE barlib)

#--------------------------------------- dependencies
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
target_link_libraries(barcode PRIVATE ${GSL_LIBRARIES})
target_link_libraries(barlib PRIVATE ${GSL_LIBRARIES})

find_package(OpenMP COMPONENTS CXX)
if (NOT OPENMP_CXX_FOUND)
    message("Disabling MULTITHREAD and MULITHREAD_FFTW options, since OpenMP support was not detected.")
    set(MULTITHREAD OFF)
    set(MULTITHREAD_FFTW OFF)
endif()

find_package(FFTW)

include_directories(${FFTW_INCLUDE_DIRS})
if (DOUBLE_PREC)
    target_link_libraries(barcode PRIVATE ${FFTW_DOUBLE_LIB})
#    target_link_libraries(barlib ${FFTW_DOUBLE_LIB})
    if (MULTITHREAD_FFTW)
        target_link_libraries(barcode PRIVATE ${FFTW_DOUBLE_OPENMP_LIB})
#        target_link_libraries(barlib ${FFTW_DOUBLE_OPENMP_LIB})
    endif()
elseif(SINGLE_PREC)
    target_link_libraries(barcode PRIVATE ${FFTW_FLOAT_LIB})
#    target_link_libraries(barlib ${FFTW_FLOAT_LIB})
    if (MULTITHREAD_FFTW)
        target_link_libraries(barcode PRIVATE ${FFTW_FLOAT_OPENMP_LIB})
#        target_link_libraries(barlib ${FFTW_FLOAT_OPENMP_LIB})
    endif()
endif()


# convert all options to preprocessor flag compilation options
if (INTEL)
    add_definitions(-DINTEL)
endif()

if (DEBUG OR (CMAKE_BUILD_TYPE MATCHES Debug) OR (CMAKE_BUILD_TYPE MATCHES RelWithDebInfo))
    add_definitions(-DDEBUG)
endif()

if (MULTITHREAD)
    add_definitions(-DMULTITHREAD)
endif()

if (MULTITHREAD_FFTW)
    add_definitions(-DMULTITHREAD_FFTW)
endif()

if (MULTITHREAD_RNG)
    add_definitions(-DMULTITHREAD_RNG)
endif()

if (NAN_DETECTION)
    add_definitions(-DNAN_DETECTION)
endif()

if (RESTART_FILE)
    add_definitions(-DRESTART_FILE)
endif()

if (MASKING)
    add_definitions(-DMASKING)
endif()

if (ZELD)
    add_definitions(-DZELD)
endif()

if (TRANSF)
    add_definitions(-DTRANSF)
endif()

if (TRANSFSC)
    add_definitions(-DTRANSFSC)
endif()

if (BARYONS)
    add_definitions(-DBARYONS)
endif()

if (RANKORD)
    add_definitions(-DRANKORD)
endif()

if (RANKORDTYPE)
    add_definitions(-DRANKORDTYPE)
endif()

if (LOGRO)
    add_definitions(-DLOGRO)
endif()

if (FOURIER_DEF_1)
    add_definitions(-DFOURIER_DEF_1)
endif()

if (FOURIER_DEF_2)
    add_definitions(-DFOURIER_DEF_2)
endif()

if (FOURIER_DEF_2_20151021)
    add_definitions(-DFOURIER_DEF_2_20151021)
endif()

if (SINGLE_PREC)
    add_definitions(-DSINGLE_PREC)
endif()

if (DOUBLE_PREC)
    add_definitions(-DDOUBLE_PREC)
endif()

if (PLANCK_SILENT)
    add_definitions(-DPLANCK_SILENT)
endif()

if (SCSMOO)
    add_definitions(-DSCSMOO)
endif()

if (GFFT)
    add_definitions(-DGFFT)
endif()

if (GFINDIFF)
    add_definitions(-DGFINDIFF)
endif()

# some additional things that must be added based on the options
if (MULTITHREAD)
    target_link_libraries(barcode PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(barlib PRIVATE OpenMP::OpenMP_CXX)
endif()

# final things that are always used
target_link_libraries(barcode PRIVATE m ncurses)
target_link_libraries(barlib PRIVATE m ncurses)

# make a link for the most recent build
# if Mercurial is present, use it to add the revision number to the output executable
# N.B.: this causes the barcode binary output file to not be present, so every build will always do a new linking step!
#find_package(Hg)
#if(HG_FOUND)
#    get_filename_component(OUTPUT_FN $<TARGET_FILE_NAME:barcode> NAME)
#    add_custom_command(TARGET barcode
#            POST_BUILD
#            COMMAND bash ${CMAKE_SOURCE_DIR}/cmake/hg_rev_filename.sh ${PROJECT_BINARY_DIR} ${OUTPUT_FN}
#            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#            COMMENT "Renamed output binary to include Mercurial revision information and created symbolic link barcode_latest_build to new path.")
#endif(HG_FOUND)


##################
# OTHER BINARIES #
##################
include(AddDefaultBinary)

#############
### TOOLS ###
#############

## density
ADD_DEFAULT_BINARY(density tools/density.cc)

## convert primordial density fluctuation field ("lagrangian") to z=0 density field ("eulerian")
ADD_DEFAULT_BINARY(LAG2EULer tools/LAG2EULer.cc)

## correlation function calculator
ADD_DEFAULT_BINARY(corr_fct tools/corr_fct.cc)

## powspec
ADD_DEFAULT_BINARY(powspec tools/powspec.cc)

## 2D correlation function generators (with or without interpolating resolution upscaling)
ADD_DEFAULT_BINARY(2D_corr_fct tools/2D_corr_fct.cc)
ADD_DEFAULT_BINARY(2D_corr_fct_interp tools/2D_corr_fct_interp.cc)
target_link_libraries(2D_corr_fct_interp PRIVATE ${FFTW_DOUBLE_LIB})
if (MULTITHREAD_FFTW)
    target_link_libraries(2D_corr_fct_interp PRIVATE ${FFTW_DOUBLE_OPENMP_LIB})
endif()

## 2D power spectrum calculator
ADD_DEFAULT_BINARY(2D_powspec tools/2D_powspec.cc)

## tools for upscaling resolution (with poisson sampling or interpolation)
ADD_DEFAULT_BINARY(interp_upres tools/interp_upres.cc)
ADD_DEFAULT_BINARY(poisson_upres tools/poisson_upres.cc)
