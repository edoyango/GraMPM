#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "grampm.hpp"
#include <algorithm>

const double dcell = 0.1, mingrid[3] {0., 0., 0.}, maxgrid[3] {0.2, 0.3, 0.4};
const kernel_linear_bspline<double> knl(dcell);
GraMPM<double> myGraMPM(24, maxgrid, mingrid, dcell, knl);

TEST_CASE( "Particle position getters + setters work", "[myGraMPM]") {

    // set coordinates
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 4; ++k) {
                myGraMPM.p.setX(12*i+4*j+k, (i+0.5)*dcell);
                myGraMPM.p.setY(12*i+4*j+k, (j+0.5)*dcell);
                myGraMPM.p.setZ(12*i+4*j+k, (k+0.5)*dcell);
            }
        }
    }

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 4; ++k) {
                REQUIRE(myGraMPM.p.getX(12*i+4*j+k) == (i+0.5)*dcell);
                REQUIRE(myGraMPM.p.getY(12*i+4*j+k) == (j+0.5)*dcell);
                REQUIRE(myGraMPM.p.getZ(12*i+4*j+k) == (k+0.5)*dcell);
            }
        }
    }
}

TEST_CASE( "Particle property getters + setters work", "[myGraMPM]") {

    for (int i = 0; i < 24; ++i) {
        myGraMPM.p.setMass(i, 1000.*i);
        REQUIRE(myGraMPM.p.getMass(i) == 1000.*i);
    }
}

TEST_CASE( "Particle append works", "[myGraMPM]") {

     // set coordinates
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 4; ++k) {
                myGraMPM.p.appendPosition(
                    (i+0.5)*dcell,
                    (j+0.5)*dcell,
                    (k+0.5)*dcell);
            }
        }
    }

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 4; ++k) {
                particle<double> p = myGraMPM.p.getParticle(12*i+4*j+k);
                REQUIRE(p.x==(i+0.5)*dcell);
                REQUIRE(p.y==(j+0.5)*dcell);
                REQUIRE(p.z==(k+0.5)*dcell);
            }
        }
    }

    REQUIRE(myGraMPM.p.getNtotal()==24);
}