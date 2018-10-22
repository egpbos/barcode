/*
 * Barcode
 * Copyright E.G.P. Bos
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <define_opt.h>  // real_prec
#include <fftw_array.h>
#include <IOfunctionsGen.h>

//#include <experimental/filesystem> // not included in macOS stl yet, arggh
#include <fstream> // istream to test for file existence
#include <cstdio>  // delete
#include <exception> // runtime_error
#include <numeric>  // iota
#include <algorithm> // for_each

#include <catch.hpp>

bool file_exists(const char *filename) {
  std::ifstream infile(filename);
  return infile.good();
}

void remove_file(const char *filename) {
  if(remove(filename) != 0) {
    throw std::runtime_error("for some reason the file could not be removed");
  }
}


SCENARIO("The array IO functions read_array and write_array should respectively read and produce binary files from arrays.") {
  GIVEN("An fftw_array<real_prec> of shape 2x2x2 containing numbers 1-8") {
    unsigned long N = 2 * 2 * 2;
    fftw_array<real_prec> array(N);
    std::iota(array + 0, array + N, 1);

    WHEN("write_array is called on it") {
      write_array("io_array_test", array, N);
      THEN("a file should be created with the given filename") {
//        REQUIRE(std::experimental::filesystem::exists("io_array_test.dat"));
        REQUIRE(file_exists("io_array_test.dat"));
      }
    }

    AND_WHEN("we read the file with read_array") {
      fftw_array<real_prec> in_array(N);
      read_array("io_array_test", in_array, N);
      THEN("the data should equal the original fftw_array's data") {
        for (std::size_t i = 0ul; i < 8; ++i) {
          REQUIRE(array[i] == in_array[i]);
        }
        remove_file("io_array_test.dat");
      }
    }
  }
}


// This is not an actual test, but just used to create a persistent
// test file that we add to the Git repo for testing whether new
// versions still load the data the same as older versions.
TEST_CASE("create persistent data file", "[!hide]") {
  unsigned long N = 2 * 2 * 2;
  fftw_array<real_prec> array(N);
  array[0] = 18012.18201;
  array[1] = 280.22;
  array[2] = 300021.850;
  array[3] = 3.14;
  array[4] = 2.;
  array[5] = 333888.;
  array[6] = 807520.20;
  array[7] = 170412.0;

  write_array("data/io_array", array, N);
}

// The above file is then loaded in this test:
SCENARIO("The function read_array must read binary files the same way as in previous versions.") {
  GIVEN("A persistent data file saved in the repo and an array of size 8") {
    unsigned long N = 2 * 2 * 2;
    fftw_array<real_prec> array(N);
    WHEN("the file is read into the array with read_array") {
      read_array("data/io_array", array, N);
      THEN("it should give the same numbers as in previous versions") {
        REQUIRE(array[0] == 18012.18201);
        REQUIRE(array[1] == 280.22);
        REQUIRE(array[2] == 300021.850);
        REQUIRE(array[3] == 3.14);
        REQUIRE(array[4] == 2.);
        REQUIRE(array[5] == 333888.);
        REQUIRE(array[6] == 807520.20);
        REQUIRE(array[7] == 170412.0);
      }
    }
  }
}
