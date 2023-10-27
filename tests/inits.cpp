// tests to check that GraMPM and members have been initialized correctly

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "grampm.hpp"
#include <algorithm>
#include <grampm_kernels.hpp>

const double dcell = 0.1, mingrid[3] {-0.1, 0.05, 0.25}, maxgrid[3] {0.2, 0.4, 0.7};
const kernel_linear_bspline<double> knl(dcell);
GraMPM<double> myGraMPM(10, maxgrid, mingrid, dcell, knl);

TEST_CASE( "Grid is initialized properly", "[myGraMPM]" ) {
    REQUIRE( myGraMPM.g.ngridx == 5 );
    REQUIRE( myGraMPM.g.ngridy == 5 );
    REQUIRE( myGraMPM.g.ngridz == 6 );
}

TEST_CASE( "Particles is initalized properly", "[myGraMPM]") {
    
    for (int i = 0; i < 10; ++i) {
        REQUIRE(myGraMPM.p.getX(i)==0.);
        REQUIRE(myGraMPM.p.getY(i)==0.);
        REQUIRE(myGraMPM.p.getZ(i)==0.);
        REQUIRE(myGraMPM.p.getMass(i)==0.);
        REQUIRE(myGraMPM.p.getRavelledGId(i)==0);
    }
    REQUIRE(myGraMPM.p.getNtotal()==0);
}