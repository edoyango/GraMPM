// Ensuring that functions to determine grid neighbours do their job

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "grampm.hpp"
#include <algorithm>
#include <grampm_kernels.hpp>

const double dcell = 0.1, mingrid[3] {-0.1, 0.05, 0.25}, maxgrid[3] {0.2, 0.4, 0.7};
const kernel_linear_bspline<double> knl(dcell);
GraMPM<double> myGraMPM(5, maxgrid, mingrid, dcell, knl);

TEST_CASE( "Correct ravelling of particles' grid indices", "[myGraMPM]") {
    myGraMPM.p.appendPosition(0.2, 0.4, 0.7);
    myGraMPM.calcPGridCells(0, 1);
    REQUIRE(myGraMPM.p.getRavelledGId(0)==90+18+4);
    myGraMPM.p.clear();
}

TEST_CASE( "Correct assignment and unravelling of particle/grid node neighbours and distances", "[myGraMPM]") {

    double dx = (maxgrid[0]-mingrid[0])/4;
    double dy = (maxgrid[1]-mingrid[1])/4;
    double dz = (maxgrid[2]-mingrid[2])/4;

    // initialize particles' position
    for (int i = 0; i < 5; ++i) {
        myGraMPM.p.appendPosition(
            mingrid[0] + i*dx,
            mingrid[1] + i*dy,
            mingrid[2] + i*dz
        );
    }

    // calculate particles' grid cell indices
    myGraMPM.calcPGridCells(0, 5);

    const int correctIdx[5] {0, 0, 1, 2, 3};
    const int correctIdy[5] {0, 0, 1, 2, 3};
    const int correctIdz[5] {0, 1, 2, 3, 4};

    for (int i = 0; i < 5; ++i) {
        int idx, idy, idz;
        myGraMPM.getPGridId(i, idx, idy, idz);
        REQUIRE(idx == correctIdx[i]);
        REQUIRE(idy == correctIdy[i]);
        REQUIRE(idz == correctIdz[i]);
    }

    // finding neighbour nodes (1 cell radius)

}