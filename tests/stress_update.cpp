#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "grampm.hpp"
#include <algorithm>
#include <grampm_kernels.hpp>
#include <array>

TEST_CASE("Check correct Hooke's law") {
    
    const double dcell = 0.2;
    GraMPM::grid<double> g(0., 0., 0., 0.99, 1.99, 2.99, dcell);
    GraMPM::kernel_linear_bspline<double> knl(dcell);

    std::array<double, 3> bf {0., 0., 0.};
    GraMPM::particle_system<double> p(1, bf, g, knl);

    p.set_strainratexx(0, 1.);
    p.set_strainrateyy(0, 2.);
    p.set_strainratezz(0, 3.);
    p.set_strainratexy(0, 4.);
    p.set_strainratexz(0, 5.);
    p.set_strainrateyz(0, 6.);

    p.set_E(100.);
    p.set_v(0.25);

    p.update_stress(0.1);

    REQUIRE(p.sigmaxx(0)==32.);
    REQUIRE(p.sigmayy(0)==40.);
    REQUIRE(p.sigmazz(0)==48.);
    REQUIRE(p.sigmaxy(0)==32.);
    REQUIRE(p.sigmaxz(0)==40.);
    REQUIRE(p.sigmayz(0)==48.);
}

TEST_CASE("Check correct spin rate") {
    
    const double dcell = 0.2;
    GraMPM::grid<double> g(0., 0., 0., 0.99, 1.99, 2.99, dcell);
    GraMPM::kernel_linear_bspline<double> knl(dcell);

    std::array<double, 3> bf {0., 0., 0.};
    GraMPM::particle_system<double> p(1, bf, g, knl);

    p.set_spinratexy(0, -0.8);
    p.set_spinratexz(0, -0.9);
    p.set_spinrateyz(0, -1.0);
    p.set_sigmaxx(0, 100.);
    p.set_sigmayy(0, 200.);
    p.set_sigmazz(0, 300.);
    p.set_sigmaxy(0, 400.);
    p.set_sigmaxz(0, 500.);
    p.set_sigmayz(0, 600.);

    CHECK(p.spinratexy(0)==-0.8);
    CHECK(p.spinratexz(0)==-0.9);
    CHECK(p.spinrateyz(0)==-1.0);
    CHECK(p.sigmaxx(0)==100.);
    CHECK(p.sigmayy(0)==200.);
    CHECK(p.sigmazz(0)==300.);
    CHECK(p.sigmaxy(0)==400.);
    CHECK(p.sigmaxz(0)==500.);
    CHECK(p.sigmayz(0)==600.);

    p.update_stress(0.1);

    REQUIRE(p.sigmaxx(0) == 254.);
    REQUIRE(p.sigmayy(0) == 256.);
    REQUIRE(p.sigmazz(0) == 90.);
    REQUIRE(p.sigmaxy(0) == 512.);
    REQUIRE(p.sigmaxz(0) == 526.);
    REQUIRE(p.sigmayz(0) == 534.);

}