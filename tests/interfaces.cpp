// tests to check that GraMPM and members have been initialized correctly

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "grampm.hpp"
#include <algorithm>
#include <grampm_kernels.hpp>
#include <array>

const std::array<double, 3> mingridx_in {-0.1, 0.05, 0.15}, maxgridx_in {0.1, 0.3, 0.5};
const double dcell_in = 0.1;
const std::array<double, 3> bf {1., 2., 3.};

// test array set interface
GraMPM::grid<double> grid(mingridx_in, maxgridx_in, dcell_in);

GraMPM::kernel_linear_bspline<double> knl(dcell_in);

GraMPM::particle_system<double> particles(5, bf, grid, knl);

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
        REQUIRE(grid1.momentumx(i)==0.);
        REQUIRE(grid1.momentumy(i)==0.);
        REQUIRE(grid1.momentumz(i)==0.);
        REQUIRE(grid1.forcex(i)==0.);
        REQUIRE(grid1.forcey(i)==0.);
        REQUIRE(grid1.forcez(i)==0.);
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

    // test that body force vector is set to the correct values
    std::array<double, 3> bf_out = particles.body_force();
    REQUIRE(bf_out[0] == 1.);
    REQUIRE(bf_out[1] == 2.);
    REQUIRE(bf_out[2] == 3.);

    // check that the particles instance has 5 particles, but all zeroed
    for (int i = 0; i < 5; ++i) {
        REQUIRE(particles.x(i)==0.);
        REQUIRE(particles.y(i)==0.);
        REQUIRE(particles.z(i)==0.);
        REQUIRE(particles.vx(i)==0.);
        REQUIRE(particles.vy(i)==0.);
        REQUIRE(particles.vz(i)==0.);
        REQUIRE(particles.mass(i)==0.);
        REQUIRE(particles.sigmaxx(i)==0.);
        REQUIRE(particles.sigmayy(i)==0.);
        REQUIRE(particles.sigmazz(i)==0.);
        REQUIRE(particles.sigmaxy(i)==0.);
        REQUIRE(particles.sigmaxz(i)==0.);
        REQUIRE(particles.sigmayz(i)==0.);
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
        particles.set_vx(i, -4.*i);
        REQUIRE(particles.vx(i)==-4.*i);
        particles.set_vy(i, -5.*i);
        REQUIRE(particles.vy(i)==-5.*i);
        particles.set_vz(i, -6.*i);
        REQUIRE(particles.vz(i)==-6.*i);
        particles.set_mass(i, 30.*i);
        REQUIRE(particles.mass(i)==30.*i);
        particles.set_rho(i, 40.*i);
        REQUIRE(particles.rho(i)==40.*i);
        particles.set_sigmaxx(i, -0.1*i);
        REQUIRE(particles.sigmaxx(i)==-0.1*i);
        particles.set_sigmayy(i, -0.2*i);
        REQUIRE(particles.sigmayy(i)==-0.2*i);
        particles.set_sigmazz(i, -0.3*i);
        REQUIRE(particles.sigmazz(i)==-0.3*i);
        particles.set_sigmaxy(i, -0.4*i);
        REQUIRE(particles.sigmaxy(i)==-0.4*i);
        particles.set_sigmaxz(i, -0.5*i);
        REQUIRE(particles.sigmaxz(i)==-0.5*i);
        particles.set_sigmayz(i, -0.6*i);
        REQUIRE(particles.sigmayz(i)==-0.6*i);
    }

    particles.set_body_force(2., 4., 6.);
    REQUIRE(particles.body_force(0)==2.);
    REQUIRE(particles.body_force(1)==4.);
    REQUIRE(particles.body_force(2)==6.);

    // check that aggregate interface works
    std::vector<GraMPM::particle<double>> pv;
    for (int i = 0; i < 5; ++i) {
        pv.push_back(
            GraMPM::particle<double>(i, 2.*i, 3.*i, 4.*i, 5.*i, 6.*i, 10.*i, 100.*i, -0.1*i, -0.2*i, -0.3*i, -0.4*i, 
                -0.5*i, -0.6*i)
        );
    }
    GraMPM::particle_system<double> particles2(pv, bf, grid, knl);

    for (int i = 0; i < 5; ++i) {
        REQUIRE(particles2.x(i)==i);
        REQUIRE(particles2.y(i)==2.*i);
        REQUIRE(particles2.z(i)==3.*i);
        REQUIRE(particles2.vx(i)==4.*i);
        REQUIRE(particles2.vy(i)==5.*i);
        REQUIRE(particles2.vz(i)==6.*i);
        REQUIRE(particles2.mass(i)==10.*i);
        REQUIRE(particles2.rho(i)==100.*i);
        REQUIRE(particles2.sigmaxx(i)==-0.1*i);
        REQUIRE(particles2.sigmayy(i)==-0.2*i);
        REQUIRE(particles2.sigmazz(i)==-0.3*i);
        REQUIRE(particles2.sigmaxy(i)==-0.4*i);
        REQUIRE(particles2.sigmaxz(i)==-0.5*i);
        REQUIRE(particles2.sigmayz(i)==-0.6*i);
    }

    // check aggregate getters
    for (int i = 0; i < 5; ++i) {
        GraMPM::particle<double> p = particles2.at(i);
        REQUIRE(p.x==i);
        REQUIRE(p.y==2.*i);
        REQUIRE(p.z==3.*i);
        REQUIRE(p.vx==4.*i);
        REQUIRE(p.vy==5.*i);
        REQUIRE(p.vz==6.*i);
        REQUIRE(p.mass==10.*i);
        REQUIRE(p.rho==100.*i);
        REQUIRE(p.sigmaxx==-0.1*i);
        REQUIRE(p.sigmayy==-0.2*i);
        REQUIRE(p.sigmazz==-0.3*i);
        REQUIRE(p.sigmaxy==-0.4*i);
        REQUIRE(p.sigmaxz==-0.5*i);
        REQUIRE(p.sigmayz==-0.6*i);
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
        REQUIRE(particles.vx(i)==0.);
        REQUIRE(particles.vy(i)==0.);
        REQUIRE(particles.vz(i)==0.);
        REQUIRE(particles.mass(i)==0.);
        REQUIRE(particles.rho(i)==0.);
        REQUIRE(particles.sigmaxx(i)==0.);
        REQUIRE(particles.sigmayy(i)==0.);
        REQUIRE(particles.sigmazz(i)==0.);
        REQUIRE(particles.sigmaxy(i)==0.);
        REQUIRE(particles.sigmaxz(i)==0.);
        REQUIRE(particles.sigmayz(i)==0.);
        REQUIRE(particles.ravelled_grid_idx(i)==0);
    }
}