/*
 * Barcode
 * Copyright E.G.P. Bos
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#include <string>

#include "define_opt.h" // real_prec
#include <ini_reader.hpp>

#include <catch.hpp>


SCENARIO("Reading a parameter input file.") {
  GIVEN("An filename for a parameter input file") {
    std::string filename("data/input.par");

    WHEN("we load the file") {
      parameter_inifile params(filename);
      THEN("we should be able to correctly load its contents") {
        REQUIRE(params.find<bool>("correct_delta") == true);
        REQUIRE(params.find<int>("calc_h") == 2);
        REQUIRE(params.find<int>("particle_kernel") == 0);
        REQUIRE(params.find<real_prec>("particle_kernel_h_rel") == 1.);
        REQUIRE(params.find<int>("inputmode") == 0);
        REQUIRE(params.find<int>("seed") == 1);
        REQUIRE(params.find<bool>("random_test") == true);
        REQUIRE(params.find<bool>("random_test_rsd") == false);
        REQUIRE(params.find<int>("window_type") == 1);
        REQUIRE(params.find<int>("data_model") == 0);
        REQUIRE(params.find<bool>("negative_obs") == false);
        REQUIRE(params.find<int>("likelihood") == 1);
        REQUIRE(params.find<int>("prior") == 0);
        REQUIRE(params.find<int>("sfmodel") == 1);
        REQUIRE(params.find<bool>("rsd_model") == false);
        REQUIRE(params.find<real_prec>("sigma_min") == 1.0);
        REQUIRE(params.find<real_prec>("sigma_fac") == 0.0);
        REQUIRE(params.find<real_prec>("delta_min") == -0.999);
        REQUIRE(params.find<int>("initial_guess") == 0);
        REQUIRE(params.find<std::string>("initial_guess_file") == "deltaLAGtest");
        REQUIRE(params.find<int>("initial_guess_smoothing_type") == 1);
        REQUIRE(params.find<real_prec>("initial_guess_smoothing_scale") == 20.);
        REQUIRE(params.find<real_prec>("N_eps_fac") == 8.0);
        REQUIRE(params.find<int>("eps_fac_update_type") == 3);
        REQUIRE(params.find<real_prec>("eps_fac") == 0.0);
        REQUIRE(params.find<real_prec>("eps_fac_initial") == 0.5);
        REQUIRE(params.find<int>("eps_fac_power") == 2);
        REQUIRE(params.find<int>("N_a_eps_update") == 100);
        REQUIRE(params.find<real_prec>("acc_min") == 0.6);
        REQUIRE(params.find<real_prec>("acc_max") == 0.7);
        REQUIRE(params.find<int>("eps_down_smooth") == 5);
        REQUIRE(params.find<int>("eps_up_fac") == 1);
        REQUIRE(params.find<int>("mass_type") == 1);
        REQUIRE(params.find<int>("massnum_burn") == 0);
        REQUIRE(params.find<int>("massnum_post") == 0);
        REQUIRE(params.find<int>("outnum") == 10);
        REQUIRE(params.find<int>("outnum_ps") == 10);
        REQUIRE(params.find<std::string>("file") == "/path/to/input/file.dat");
        REQUIRE(params.find<std::string>("filec") == "file.dat");
        REQUIRE(params.find<bool>("readPS") == true);
        REQUIRE(params.find<std::string>("fnamePS") == "/path/to/power/spectrum/like/for/instance/WMAP7_CAMB.dat");
        REQUIRE(params.find<std::string>("dir") == "./data/");
        REQUIRE(params.find<real_prec>("slength") == 4.);
        REQUIRE(params.find<int>("Nx") == 64);
        REQUIRE(params.find<real_prec>("Lx") == 200.);
        REQUIRE(params.find<real_prec>("z") == .0);
        REQUIRE(params.find<int>("N_bin") == 200);
        REQUIRE(params.find<int>("N_Gibbs") == 10000);
        REQUIRE(params.find<int>("total_steps_lim") == 0);
        REQUIRE(params.find<int>("masskernel") == 3);
        REQUIRE(params.find<real_prec>("xllc") == 0.);
        REQUIRE(params.find<real_prec>("yllc") == 0.);
        REQUIRE(params.find<real_prec>("zllc") == 0.);
        REQUIRE(params.find<real_prec>("xobs") == 90.);
        REQUIRE(params.find<real_prec>("yobs") == 90.);
        REQUIRE(params.find<real_prec>("zobs") == 90.);
        REQUIRE(params.find<bool>("planepar") == true);
        REQUIRE(params.find<bool>("periodic") == true);
        REQUIRE(params.find<real_prec>("mass_factor") == 1.);
        REQUIRE(params.find<real_prec>("grad_psi_prior_factor") == 1.);
        REQUIRE(params.find<real_prec>("grad_psi_likeli_factor") == 1.);
        REQUIRE(params.find<bool>("grad_psi_prior_conjugate") == false);
        REQUIRE(params.find<bool>("grad_psi_likeli_conjugate") == false);
        REQUIRE(params.find<bool>("grad_psi_prior_times_i") == false);
        REQUIRE(params.find<bool>("grad_psi_likeli_times_i") == false);
        REQUIRE(params.find<bool>("div_dH_by_N") == false);
        REQUIRE(params.find<real_prec>("deltaQ_factor") == 1.0);
        REQUIRE(params.find<real_prec>("s_eps_total_fac") == 158.0);
        REQUIRE(params.find<int>("s_eps_total_Nx_norm") == 64);
        REQUIRE(params.find<real_prec>("s_eps_total_scaling") == 0.5);
      }
    }
  }
}
