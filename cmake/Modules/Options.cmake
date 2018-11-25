#
# Barcode
# Copyright E.G.P. Bos and F.S. Kitaura
#
# Distributed under the terms of the MIT License.
# The full license is in the file LICENSE, distributed with this software.
#

option(INTEL "Use intel library and compiler" OFF)

#--------------------------------------- Basic operation mode of code
option(DEBUG "Debugging mode" OFF)
option(MULTITHREAD "Enables OpenMP directives on most for-loops, except once with RNGs." ON)
option(MULTITHREAD_FFTW "Use multi-threaded version of the FFTW3 library." ON)
option(MULTITHREAD_RNG "Enable OpenMP directives on for-loops containing random number generators. This causes an extra degree of randomness, due to the order in which cores access the RNG, which will differ every run, because the scheduler controls this and other running processes will have an (unpredictable) influence on this." OFF)

if (NOT OPENMP_CXX_FOUND)
    message("Disabling MULTITHREAD option, since OpenMP support was not detected.")
    set(MULTITHREAD OFF)
endif()

if (MULTITHREAD_FFTW)
    if     (DOUBLE_PREC AND NOT FFTW_DOUBLE_OPENMP_LIB_FOUND AND NOT FFTW_DOUBLE_THREADS_LIB_FOUND)
        message("Disabling MULTITHREAD_FFTW option, since the double precision OpenMP and pthreads FFTW libraries were not detected.")
        set(MULTITHREAD_FFTW OFF)
    elseif (SINGLE_PREC AND NOT FFTW_FLOAT_OPENMP_LIB_FOUND AND NOT FFTW_FLOAT_THREADS_LIB_FOUND)
        message("Disabling MULTITHREAD_FFTW option, since the single precision OpenMP and pthreads FFTW libraries were not detected.")
        set(MULTITHREAD_FFTW OFF)
    endif()
endif()


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

#--------------------------------------- Things for special behaviour
option(SCSMOO "[NOT TESTED!] Smooth the spherical collapse part of the displacement field with a smoothing length kthsc2. This variable is currently not implemented in the code. kthsc is implemented and is equal to kth, which equals ksmooth from input.par. Something similar should probably be defined for this." OFF)

OPTION(COVERAGE "Enable coverage compile flags (gcc only!)" OFF)

#--------------------------------------- NEEDS ATTENTION

if (COVERAGE)
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        message(FATAL_ERROR "Coverage can only be used with GNU and Clang compilers!")
    endif()

    message("Adding coverage compile flags")
    #    set_target_properties(${TARGET} PROPERTIES COMPILE_FLAGS -g -O0 --coverage)
    #    add_compile_options(-g -O0 --coverage)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")     # debug, no optimisation
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage") # enabling coverage
endif()

# CONSIDER EITHER TAKING THESE OPTIONS OUT OR ACTUALLY GIVING THEM CODE-WIDE IMPACT:
option(GFFT "Use FFT based derivative (NOT REALLY USED, currently only used in non-Zel'dovich Lag2Eul via calc_m2v_mem and in EigenValuesTweb)." OFF)
option(GFINDIFF "Use finite difference derivative (NOT REALLY USED, currently only used in non-Zel'dovich Lag2Eul via calc_m2v_mem and in EigenValuesTweb)." ON)

if ( GFFT AND GFINDIFF )
    message(FATAL_ERROR "Cannot use GFFT and GFINDIFF at same time!")
endif()

if ( NOT (GFINDIFF OR GFFT) )
    message(FATAL_ERROR "Must choose one of GFFT and GFINDIFF!")
endif()



# convert all options to preprocessor flag compilation options
function(TARGET_ADD_OPTIONS target)
    if (INTEL)
        target_compile_definitions(${target} PUBLIC INTEL)
    endif()

    if (DEBUG OR (CMAKE_BUILD_TYPE MATCHES Debug) OR (CMAKE_BUILD_TYPE MATCHES RelWithDebInfo))
        target_compile_definitions(${target} PUBLIC DEBUG)
    endif()

    if (MULTITHREAD)
        target_compile_definitions(${target} PUBLIC MULTITHREAD)
    endif()

    if (MULTITHREAD_FFTW)
        target_compile_definitions(${target} PUBLIC MULTITHREAD_FFTW)
    endif()

    if (MULTITHREAD_RNG)
        target_compile_definitions(${target} PUBLIC MULTITHREAD_RNG)
    endif()

    if (OPENMP_CXX_FOUND)
        target_compile_definitions(${target} PUBLIC WITH_OPENMP)
    endif()

    if (NAN_DETECTION)
        target_compile_definitions(${target} PUBLIC NAN_DETECTION)
    endif()

    if (RESTART_FILE)
        target_compile_definitions(${target} PUBLIC RESTART_FILE)
    endif()

    if (MASKING)
        target_compile_definitions(${target} PUBLIC MASKING)
    endif()

    if (ZELD)
        target_compile_definitions(${target} PUBLIC ZELD)
    endif()

    if (TRANSF)
        target_compile_definitions(${target} PUBLIC TRANSF)
    endif()

    if (TRANSFSC)
        target_compile_definitions(${target} PUBLIC TRANSFSC)
    endif()

    if (BARYONS)
        target_compile_definitions(${target} PUBLIC BARYONS)
    endif()

    if (RANKORD)
        target_compile_definitions(${target} PUBLIC RANKORD)
    endif()

    if (RANKORDTYPE)
        target_compile_definitions(${target} PUBLIC RANKORDTYPE)
    endif()

    if (LOGRO)
        target_compile_definitions(${target} PUBLIC LOGRO)
    endif()

    if (FOURIER_DEF_1)
        target_compile_definitions(${target} PUBLIC FOURIER_DEF_1)
    endif()

    if (FOURIER_DEF_2)
        target_compile_definitions(${target} PUBLIC FOURIER_DEF_2)
    endif()

    if (FOURIER_DEF_2_20151021)
        target_compile_definitions(${target} PUBLIC FOURIER_DEF_2_20151021)
    endif()

    if (SINGLE_PREC)
        target_compile_definitions(${target} PUBLIC SINGLE_PREC)
    endif()

    if (DOUBLE_PREC)
        target_compile_definitions(${target} PUBLIC DOUBLE_PREC)
    endif()

    if (SCSMOO)
        target_compile_definitions(${target} PUBLIC SCSMOO)
    endif()

    if (GFFT)
        target_compile_definitions(${target} PUBLIC GFFT)
    endif()

    if (GFINDIFF)
        target_compile_definitions(${target} PUBLIC GFINDIFF)
    endif()
endfunction(TARGET_ADD_OPTIONS)


function(TARGET_ADD_FFTW target)
    if (DOUBLE_PREC)
        target_link_libraries(${target} PRIVATE ${FFTW_DOUBLE_LIB})
    elseif(SINGLE_PREC)
        target_link_libraries(${target} PRIVATE ${FFTW_FLOAT_LIB})
    endif()

    if (MULTITHREAD_FFTW)
        if     (DOUBLE_PREC AND OPENMP_CXX_FOUND AND FFTW_DOUBLE_OPENMP_LIB_FOUND)
            target_link_libraries(${target} PRIVATE ${FFTW_DOUBLE_OPENMP_LIB})
        elseif (DOUBLE_PREC AND FFTW_DOUBLE_THREADS_LIB_FOUND) # pthreads version as fallback
            target_link_libraries(${target} PRIVATE ${FFTW_DOUBLE_THREADS_LIB})
        elseif (SINGLE_PREC AND OPENMP_CXX_FOUND AND FFTW_FLOAT_OPENMP_LIB_FOUND)
            target_link_libraries(${target} PRIVATE ${FFTW_FLOAT_OPENMP_LIB})
        elseif (SINGLE_PREC AND FFTW_FLOAT_THREADS_LIB_FOUND)
            target_link_libraries(${target} PRIVATE ${FFTW_FLOAT_THREADS_LIB})
        else()
            message("Something went wrong, shouldn't have reached this point!")
        endif()
    endif()
endfunction(TARGET_ADD_FFTW)
