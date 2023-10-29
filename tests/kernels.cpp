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
    REQUIRE(knl.w(0., 0., 0.)==1.);

    // check kernel gradients return correct values
    double dwdx, dwdy, dwdz;

    // checking gradient is 0 at 0
    // for 
    knl.dwdx(0., 0., 0., dwdx, dwdy, dwdz);
    REQUIRE(dwdx==0.);
    REQUIRE(dwdy==0.);
    REQUIRE(dwdz==0.);

    // checking gradient returns expected values in R+
    knl.dwdx(0.05, 0.06, 0.07, dwdx, dwdy, dwdz);
    REQUIRE(std::round(dwdx*10)==-12.);
    REQUIRE(std::round(dwdy*10)==-15.);
    REQUIRE(std::round(dwdz)==-2.);

    // checking symmetry about origin
    knl.dwdx(-0.05, -0.06, -0.07, dwdx, dwdy, dwdz);
    REQUIRE(std::round(dwdx*10)==12.);
    REQUIRE(std::round(dwdy*10)==15.);
    REQUIRE(std::round(dwdz)==2.);
}

GraMPM::kernel_cubic_bspline<double> knlc(0.1);

TEST_CASE("Testing cubic bspline returns correct values", "[knlc]") {

    // check attributes are set correctly
    REQUIRE(knlc.radius == 2.);
    REQUIRE(knlc.dcell == 0.1);

    // check kernel function returns correct values
    // below round thing is because C++ is stupid
    REQUIRE(std::round(1e7*knlc.w(0.05, 0.06, 0.07)) == std::round(1e7*0.0691787824));
    REQUIRE(std::round(1e6*knlc.w(0.15, 0.16, 0.17))==1.);
    REQUIRE(knlc.w(0.2, 0.1, 0.07)==0.);
    REQUIRE(knlc.w(0.05, 0.2, 0.07)==0.);
    REQUIRE(knlc.w(0.05, 0.06, 0.2)==0.);
    REQUIRE(knlc.w(0.2, 0.2, 0.2)==0.);
    REQUIRE(knlc.w(0., 0., 0.)==8./27.);

    // check kernel gradients return correct values
    double dwdx, dwdy, dwdz;

    // checking gradient is 0 at 0
    knl.dwdx(0., 0., 0., dwdx, dwdy, dwdz);
    REQUIRE(dwdx==0.);
    REQUIRE(dwdy==0.);
    REQUIRE(dwdz==0.);

    // checking gradient is 0 at 2*dcell
    knl.dwdx(0.2, 0.2, 0.2, dwdx, dwdy, dwdz);
    REQUIRE(dwdx==0.);
    REQUIRE(dwdy==0.);
    REQUIRE(dwdz==0.);

    // checking gradient returns expected values in R+
    knlc.dwdx(0.05, 0.06, 0.07, dwdx, dwdy, dwdz);
    REQUIRE(std::round(dwdx*1e7)==-9023319);
    REQUIRE(std::round(dwdy*1e7)==-11010771);
    REQUIRE(std::round(dwdz*1e7)==-13213181);

    // checking symmetry about origin
    knlc.dwdx(-0.05, -0.06, -0.07, dwdx, dwdy, dwdz);
    REQUIRE(std::round(dwdx*1e7)==9023319);
    REQUIRE(std::round(dwdy*1e7)==11010771);
    REQUIRE(std::round(dwdz*1e7)==13213181);
}