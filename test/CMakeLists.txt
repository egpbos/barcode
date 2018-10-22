# CTest (needed for catch_discover_tests?)
enable_testing()

# Prepare "Catch" library for other executables
set(CATCH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/ext/catch2)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

# for automatic CTest registration of Catch tests using catch_discover_tests()
include(${CATCH_INCLUDE_DIR}/Catch.cmake)

# catch main, speeds up compilation of individual tests (https://github.com/catchorg/Catch2/blob/master/docs/slow-compiles.md)
add_library(catch_main catch_main.cpp)
target_link_libraries(catch_main PUBLIC Catch)
target_compile_features(catch_main PUBLIC cxx_std_11)

# test executables
add_executable(io_array io_array.cpp)
target_link_libraries(io_array PRIVATE barlib catch_main)
catch_discover_tests(io_array WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
#target_compile_features(io_array PRIVATE cxx_std_17)
