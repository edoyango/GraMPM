#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "grampm_kernels.hpp"
#include <cmath>

GraMPM::kernel_linear_bspline<double> knl(0.1);

TEST_CASE("Testing linear bspline returns correct values", "[knl]") {

    // check attributes are set correctly
    REQUIRE(knl.radius == 1.);
    REQUIRE(knl.dcell == 0.1);

    // check kernel function returns correct values
    // below round thing is because C++ is stupid
    REQUIRE(std::round(knl.w(0.05, 0.06, 0.07)*100) == 6);
    REQUIRE(knl.w(0.1, 0.06, 0.07)==0.);
    REQUIRE(knl.w(0.05, 0.1, 0.07)==0.);
    REQUIRE(knl.w(0.05, 0.06, 0.1)==0.);
    REQUIRE(knl.w(0.1, 0.1, 0.1)==0.);

    // check kernel gradients return correct values
    double dwdx, dwdy, dwdz;
    knl.dwdx(0.05, 0.06, 0.07, dwdx, dwdy, dwdz);
    REQUIRE(std::round(dwdx*10)==-12.);
    REQUIRE(std::round(dwdy*10)==-15.);
    REQUIRE(std::round(dwdz)==-2.);

    knl.dwdx(-0.05, -0.06, -0.07, dwdx, dwdy, dwdz);
    REQUIRE(std::round(dwdx*10)==12.);
    REQUIRE(std::round(dwdy*10)==15.);
    REQUIRE(std::round(dwdz)==2.);
}