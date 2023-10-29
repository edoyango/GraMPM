// tests to check that GraMPM and members have been initialized correctly

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "grampm.hpp"
#include <algorithm>
#include <grampm_kernels.hpp>
#include <array>

const std::array<double, 3> mingridx_in {-0.1, 0.05, 0.15}, maxgridx_in {0.1, 0.3, 0.5};
const double dcell_in = 0.1;

// test array set interface
GraMPM::grid<double> grid(mingridx_in, maxgridx_in, dcell_in);

GraMPM::kernel_linear_bspline<double> knl(dcell_in);

GraMPM::particle_system<double> particles(5, grid, knl);

TEST_CASE("grid intialized correctly", "[grid]") {

    GraMPM::grid<double> grid1 {grid};

    REQUIRE(grid1.cell_size()==dcell_in);

    // test array get interface
    std::array<double, 3> mingridx_out {grid1.mingrid()}, maxgridx_out {grid1.maxgrid()};

    REQUIRE(mingridx_out[0]==mingridx_in[0]);
    REQUIRE(mingridx_out[1]==mingridx_in[1]);
    REQUIRE(mingridx_out[2]==mingridx_in[2]);

    REQUIRE(maxgridx_out[0]==maxgridx_in[0]);
    REQUIRE(maxgridx_out[1]==maxgridx_in[1]);
    REQUIRE(maxgridx_out[2]==maxgridx_in[2]);

    std::array<int, 3> ngridx_out = grid1.ngrid();
    REQUIRE(ngridx_out[0]==3);
    REQUIRE(ngridx_out[1]==4);
    REQUIRE(ngridx_out[2]==5);

    // test element-by-element get interface
    REQUIRE(grid1.mingridx()==mingridx_in[0]);
    REQUIRE(grid1.mingridy()==mingridx_in[1]);
    REQUIRE(grid1.mingridz()==mingridx_in[2]);
    
    REQUIRE(grid1.maxgridx()==maxgridx_in[0]);
    REQUIRE(grid1.maxgridy()==maxgridx_in[1]);
    REQUIRE(grid1.maxgridz()==maxgridx_in[2]);

    REQUIRE(grid1.ngridx()==3);
    REQUIRE(grid1.ngridy()==4);
    REQUIRE(grid1.ngridz()==5);

    for (int i = 0; i < grid1.ngridx()*grid1.ngridy()*grid1.ngridz(); ++i) {
        REQUIRE(grid1.mass(i)==0.);
    }

    // test element-by-element set interface
    GraMPM::grid<double> grid2(mingridx_in[0], mingridx_in[1], mingridx_in[2], maxgridx_in[0], maxgridx_in[1], 
        maxgridx_in[2], dcell_in);

    REQUIRE(grid2.cell_size()==dcell_in);

    // test array get interface
    mingridx_out = grid2.mingrid();
    REQUIRE(mingridx_out[0]==mingridx_in[0]);
    REQUIRE(mingridx_out[1]==mingridx_in[1]);
    REQUIRE(mingridx_out[2]==mingridx_in[2]);

    maxgridx_out = grid2.maxgrid();
    REQUIRE(maxgridx_out[0]==maxgridx_in[0]);
    REQUIRE(maxgridx_out[1]==maxgridx_in[1]);
    REQUIRE(maxgridx_out[2]==maxgridx_in[2]);

    ngridx_out = grid2.ngrid();
    REQUIRE(ngridx_out[0]==3);
    REQUIRE(ngridx_out[1]==4);
    REQUIRE(ngridx_out[2]==5);

    // test element-by-element get interface
    REQUIRE(grid2.mingridx()==mingridx_in[0]);
    REQUIRE(grid2.mingridy()==mingridx_in[1]);
    REQUIRE(grid2.mingridz()==mingridx_in[2]);
    
    REQUIRE(grid2.maxgridx()==maxgridx_in[0]);
    REQUIRE(grid2.maxgridy()==maxgridx_in[1]);
    REQUIRE(grid2.maxgridz()==maxgridx_in[2]);

    REQUIRE(grid2.ngridx()==3);
    REQUIRE(grid2.ngridy()==4);
    REQUIRE(grid2.ngridz()==5);
}

TEST_CASE("Particle initialized correctly", "[grid]") {

    // test that the grid attached to the particles instance is the same as the grid instance
    REQUIRE(particles.grid_address()==&grid);
    REQUIRE(particles.background_grid.ngridx()==3);
    REQUIRE(particles.background_grid.ngridy()==4);
    REQUIRE(particles.background_grid.ngridz()==5);

    // check that the particles instance has 5 particles, but all zeroed
    for (int i = 0; i < 5; ++i) {
        REQUIRE(particles.x(i)==0.);
        REQUIRE(particles.y(i)==0.);
        REQUIRE(particles.z(i)==0.);
        REQUIRE(particles.mass(i)==0.);
    }

    REQUIRE(particles.size()==5);
    REQUIRE(particles.capacity()==5);

    // check that manual setter functions work
    for (int i = 0; i < 5; ++i) {
        particles.set_x(i, -1.*i);
        REQUIRE(particles.x(i)==-1.*i);
        particles.set_y(i, -2.*i);
        REQUIRE(particles.y(i)==-2.*i);
        particles.set_z(i, -3.*i);
        REQUIRE(particles.z(i)==-3.*i);
        particles.set_mass(i, 30.*i);
        REQUIRE(particles.mass(i)==30.*i);
    }

    // check that aggregate interface works
    std::vector<GraMPM::particle<double>> pv;
    for (int i = 0; i < 5; ++i) {
        pv.push_back(
            GraMPM::particle<double>(i, 2.*i, 3.*i, 10.*i)
        );
    }
    GraMPM::particle_system<double> particles2(pv, grid, knl);

    for (int i = 0; i < 5; ++i) {
        REQUIRE(particles2.x(i)==i);
        REQUIRE(particles2.y(i)==2.*i);
        REQUIRE(particles2.z(i)==3.*i);
        REQUIRE(particles2.mass(i)==10.*i);
    }

    // check aggregate getters
    for (int i = 0; i < 5; ++i) {
        GraMPM::particle<double> p = particles2.at(i);
        REQUIRE(p.x==i);
        REQUIRE(p.y==2.*i);
        REQUIRE(p.z==3.*i);
        REQUIRE(p.mass==10.*i);
    }
}

TEST_CASE("Check clearing and resizing", "[particles]") {
    particles.clear();
    REQUIRE(particles.empty());
    particles.resize(3);
    REQUIRE(particles.size()==3);
    for (int i = 0; i < 3; ++i) {
        REQUIRE(particles.x(i)==0.);
        REQUIRE(particles.y(i)==0.);
        REQUIRE(particles.z(i)==0.);
        REQUIRE(particles.mass(i)==0.);
        REQUIRE(particles.ravelled_grid_idx(i)==0);
    }
}