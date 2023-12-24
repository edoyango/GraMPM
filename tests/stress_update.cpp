#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "grampm.hpp"
#include <algorithm>
#include <grampm_kernels.hpp>
#include <array>
#include <cmath>
#include <grampm_stress_update_functions.hpp>

TEST_CASE("Check correct Hooke's law") {
    
    const double dcell = 0.2;
    GraMPM::kernel_linear_bspline<double> knl(dcell);

    std::array<double, 3> bf {0., 0., 0.};
    GraMPM::MPM_system<double> p(1, bf, knl, std::array<double, 3>{0.}, std::array<double, 3>{0.99, 1.99, 2.99}, dcell);
    p.set_stress_update_function(
        GraMPM::stress_update::hookes_law<double>
    );

    p.p_strainratexx(0) = 1.;
    p.p_strainrateyy(0) = 2.;
    p.p_strainratezz(0) = 3.;
    p.p_strainratexy(0) = 4.;
    p.p_strainratexz(0) = 5.;
    p.p_strainrateyz(0) = 6.;

    p.m_E = 100.;
    p.m_v = 0.25;

    p.p_update_stress(0.1);

    REQUIRE(p.p_sigmaxx(0)==32.);
    REQUIRE(p.p_sigmayy(0)==40.);
    REQUIRE(p.p_sigmazz(0)==48.);
    REQUIRE(p.p_sigmaxy(0)==32.);
    REQUIRE(p.p_sigmaxz(0)==40.);
    REQUIRE(p.p_sigmayz(0)==48.);
}

TEST_CASE("Check correct Jaumann stress rate") {
    
    const double dcell = 0.2;
    // GraMPM::grid<double> g(0., 0., 0., 0.99, 1.99, 2.99, dcell);
    GraMPM::kernel_linear_bspline<double> knl(dcell);

    std::array<double, 3> bf {0., 0., 0.};
    GraMPM::MPM_system<double> p(1, bf, knl, std::array<double, 3>{0.}, std::array<double, 3>{0.99, 1.99, 2.99}, dcell);
    p.set_stress_update_function(
        GraMPM::stress_update::hookes_law<double>
    );

    p.p_spinratexy(0)= -0.8;
    p.p_spinratexz(0)= -0.9;
    p.p_spinrateyz(0)= -1.0;
    p.p_sigmaxx(0) = 100.;
    p.p_sigmayy(0) = 200.;
    p.p_sigmazz(0) = 300.;
    p.p_sigmaxy(0) = 400.;
    p.p_sigmaxz(0) = 500.;
    p.p_sigmayz(0) = 600.;

    CHECK(p.p_spinratexy(0)==-0.8);
    CHECK(p.p_spinratexz(0)==-0.9);
    CHECK(p.p_spinrateyz(0)==-1.0);
    CHECK(p.p_sigmaxx(0)==100.);
    CHECK(p.p_sigmayy(0)==200.);
    CHECK(p.p_sigmazz(0)==300.);
    CHECK(p.p_sigmaxy(0)==400.);
    CHECK(p.p_sigmaxz(0)==500.);
    CHECK(p.p_sigmayz(0)==600.);

    p.p_update_stress(0.1);

    REQUIRE(p.p_sigmaxx(0) == 254.);
    REQUIRE(p.p_sigmayy(0) == 256.);
    REQUIRE(p.p_sigmazz(0) == 90.);
    REQUIRE(p.p_sigmaxy(0) == 512.);
    REQUIRE(p.p_sigmaxz(0) == 526.);
    REQUIRE(p.p_sigmayz(0) == 534.);

}

TEST_CASE("Check DP elasto-plasticity") {
    const double dcell = 0.2;
    // GraMPM::grid<double> g(0., 0., 0., 0.99, 1.99, 2.99, dcell);
    GraMPM::kernel_linear_bspline<double> knl(dcell);

    std::array<double, 3> bf {0., 0., 0.};
    GraMPM::MPM_system<double> p(1, bf, knl, std::array<double, 3>{0.}, std::array<double, 3>{0.99, 1.99, 2.99}, dcell);
    p.set_stress_update_function(GraMPM::stress_update::drucker_prager_elastoplastic<double>);
    const double pi = std::acos(-1.);
    p.set_DP_params(pi/4., pi/36., 0.);

    double phi, psi, cohesion, alpha_phi, alpha_psi, k_c;
    p.DP_params(phi, psi, cohesion, alpha_phi, alpha_psi, k_c);
    REQUIRE(phi==pi/4.);
    REQUIRE(alpha_phi==2.*std::sin(phi)/(std::sqrt(3.)*(3.-std::sin(phi))));
    REQUIRE(psi==pi/36.);
    REQUIRE(alpha_psi==2.*std::sin(psi)/(std::sqrt(3.)*(3.-std::sin(phi))));
    REQUIRE(cohesion==0.);
    REQUIRE(k_c==0.);

    // check elasto-plastic compression
    p.p_strainratexx(0) = -1.;
    p.p_strainrateyy(0) = -2.;
    p.p_strainratezz(0) = -3.;
    p.p_strainratexy(0) = -4.;
    p.p_strainratexz(0) = -5.;
    p.p_strainrateyz(0) = -6.;

    p.m_E = 100.;
    p.m_v = 0.25;

    p.p_update_stress(0.1);

    REQUIRE(std::round(p.p_sigmaxx(0)*1e8)==-3952509695.);
    REQUIRE(std::round(p.p_sigmayy(0)*1e8)==-4496397315.);
    REQUIRE(std::round(p.p_sigmazz(0)*1e8)==-5040284935.);
    REQUIRE(std::round(p.p_sigmaxy(0)*1e8)==-2175550481.);
    REQUIRE(std::round(p.p_sigmaxz(0)*1e8)==-2719438102.);
    REQUIRE(std::round(p.p_sigmayz(0)*1e8)==-3263325722.);

    // check tensile correction
    p.p_sigmaxx(0) = 0.;
    p.p_sigmayy(0) = 0.;
    p.p_sigmazz(0) = 0.;
    p.p_sigmaxy(0) = 0.;
    p.p_sigmaxz(0) = 0.;
    p.p_sigmayz(0) = 0.;
    p.p_strainratexx(0) = 1.;
    p.p_strainrateyy(0) = 2.;
    p.p_strainratezz(0) = 3.;
    p.p_strainratexy(0) = 0.;
    p.p_strainratexz(0) = 0.;
    p.p_strainrateyz(0) = 0.;

    p.p_update_stress(0.1);
    REQUIRE(p.p_sigmaxx(0)==0.);
    REQUIRE(p.p_sigmayy(0)==0.);
    REQUIRE(p.p_sigmazz(0)==0.);
    REQUIRE(p.p_sigmaxy(0)==0.);
    REQUIRE(p.p_sigmaxz(0)==0.);
    REQUIRE(p.p_sigmayz(0)==0.);
}