/*
 * Barcode
 * Copyright E.G.P. Bos
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <string>
#include <ini_reader.hpp>
#include <catch.hpp>

SCENARIO("Reading a parameter input file.") {
  GIVEN("An filename for a parameter input file") {
    std::string filename("data/input.par");

    WHEN("we load the file") {
      parameter_inifile params(filename);
      THEN("we should be able to correctly load its contents") {
        REQUIRE(params.find<bool>("bool_true") == true);
        REQUIRE(params.find<bool>("bool_false") == false);
        REQUIRE(params.find<bool>("bool_comment") == false);

        REQUIRE(params.find<int>("zero") == 0);
        REQUIRE(params.find<unsigned int>("zero") == 0u);
        REQUIRE(params.find<long>("zero") == 0l);
        REQUIRE(params.find<unsigned long>("zero") == 0ul);
        REQUIRE(params.find<float>("zero") == 0.f);
        REQUIRE(params.find<double>("zero") == 0.);

        REQUIRE(params.find<int>("zero_comment") == 0);
        REQUIRE(params.find<unsigned int>("zero_comment") == 0u);
        REQUIRE(params.find<long>("zero_comment") == 0l);
        REQUIRE(params.find<unsigned long>("zero_comment") == 0ul);
        REQUIRE(params.find<float>("zero_comment") == 0.f);
        REQUIRE(params.find<double>("zero_comment") == 0.);

        REQUIRE(params.find<float>("float_one") == 1.f);
        REQUIRE(params.find<double>("float_one") == 1.);
        REQUIRE(params.find<float>("float_minus") == -1.2f);
        REQUIRE(params.find<double>("float_minus") == -1.2);

        REQUIRE(params.find<std::string>("string") == "hello_no_spaces_please");
        REQUIRE(params.find<std::string>("string_comment") == "stuff");
      }
    }
  }
}
