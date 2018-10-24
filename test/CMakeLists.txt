# Prepare "Catch" library for other executables
set(CATCH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/ext/catch2)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

# catch main, speeds up compilation of individual tests (https://github.com/catchorg/Catch2/blob/master/docs/slow-compiles.md)
add_library(catch_main catch_main.cpp)
target_link_libraries(catch_main PUBLIC Catch)
target_compile_features(catch_main PUBLIC cxx_std_11)

# test executables
function(add_default_test target_name source_file)
    add_executable(${target_name} ${source_file})
    target_link_libraries(${target_name} PRIVATE barlib catch_main)
endfunction(add_default_test)

add_default_test(io_array io_array.cpp)
#target_compile_features(io_array PRIVATE cxx_std_17)
add_default_test(parameter_input_file parameter_input_file.cpp)

if(COVERAGE)
    message("Adding coverage compile flags")
    #    set_target_properties(${TARGET} PROPERTIES COMPILE_FLAGS -g -O0 --coverage)
    #    add_compile_options(-g -O0 --coverage)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")     # debug, no optimisation
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage") # enabling coverage
endif()
