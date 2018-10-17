#
# Barcode
# Copyright E.G.P. Bos and F.S. Kitaura
#
# Distributed under the terms of the MIT License.
# The full license is in the file LICENSE, distributed with this software.
#

function(ADD_DEFAULT_BINARY binary_name source)
    add_executable(${binary_name} ${source})

    # optional argument: library dependencies, barlib by default
    if (ARG3)
        set(LIBS ${ARG3})
    else()
        set(LIBS barlib)
    endif(ARG3)

    target_link_libraries(${binary_name} PRIVATE ${LIBS} ${GSL_LIBRARIES} m ncurses)

    if (DOUBLE_PREC)
        target_link_libraries(${binary_name} PRIVATE ${FFTW_DOUBLE_LIB})
        if (MULTITHREAD_FFTW)
            target_link_libraries(${binary_name} PRIVATE ${FFTW_DOUBLE_OPENMP_LIB})
        endif()
    elseif(SINGLE_PREC)
        target_link_libraries(${binary_name} PRIVATE ${FFTW_FLOAT_LIB})
        if (MULTITHREAD_FFTW)
            target_link_libraries(${binary_name} PRIVATE ${FFTW_FLOAT_OPENMP_LIB})
        endif()
    endif()

    if (MULTITHREAD)
        target_link_libraries(${binary_name} PRIVATE OpenMP::OpenMP_CXX)
    endif()
endfunction(ADD_DEFAULT_BINARY)
