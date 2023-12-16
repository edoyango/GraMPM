// Ensuring that functions to determine grid neighbours do their job

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "grampm.hpp"
#include <algorithm>
#include <grampm_kernels.hpp>
#include <array>

const double dcell = 0.1;
const std::array<double, 3> mingrid {-0.1, 0.05, 0.15}, maxgrid {0.1, 0.3, 0.5};
GraMPM::grid<double> g(mingrid, maxgrid, dcell);
GraMPM::kernel_linear_bspline<double> knl(dcell);
GraMPM::kernel_cubic_bspline<double> knlc(dcell);
GraMPM::particle_system<double> p(g, knl);

const double dx = (maxgrid[0]-mingrid[0])/4;
const double dy = (maxgrid[1]-mingrid[1])/4;
const double dz = (maxgrid[2]-mingrid[2])/4;

const int correct_idxx[5] {0, 0, 1, 1, 2};
const int correct_idxy[5] {0, 0, 1, 1, 2};
const int correct_idxz[5] {0, 0, 1, 2, 3};
const int correct_ravelled_idx[5] {0, 0, 1*4*5+1*5+1, 1*4*5+1*5+2, 2*4*5+2*5+3};

TEST_CASE( "Correct ravelling of particles' grid indices upon initialization", "[p]") {
    // create the particlesystem instance
    for (int i = 0; i < 5; ++i) {
        GraMPM::particle<double> pi(i*dx+mingrid[0], i*dy+mingrid[1], i*dz+mingrid[2], 0., 0., 0., 1., 0., 0., 0., 0., 
            0., 0., 0.);
        p.push_back(pi);
        // test that the assigned index is correct
        REQUIRE(p.ravelled_grid_idx(i)==correct_ravelled_idx[i]);
        std::array<int, 3> idx = p.grid_idx(i);
        REQUIRE(idx[0]==correct_idxx[i]);
        REQUIRE(idx[1]==correct_idxy[i]);
        REQUIRE(idx[2]==correct_idxz[i]);
    }
}

TEST_CASE( "Correct dynamic assignment and unravelling of particle grid indices", "[p]") {

    // calculate particles' grid cell indices
    p.clear();

    p.resize(5);

    for (int i = 0; i < 5; ++i) {
        p.set_x(i, i*dx+mingrid[0]);
        p.set_y(i, i*dy+mingrid[1]);
        p.set_z(i, i*dz+mingrid[2]);
        CHECK(p.ravelled_grid_idx(i)==0);
    }

    p.update_particle_to_cell_map(0, 5);

    for (int i = 0; i < 5; ++i) {
        std::array<int, 3> idx = p.grid_idx(i);
        REQUIRE(idx[0]==correct_idxx[i]);
        REQUIRE(idx[1]==correct_idxy[i]);
        REQUIRE(idx[2]==correct_idxz[i]);
    }

}

TEST_CASE("Correct determination of grid node neighbours (radius=1)", "[p]") {

    p.map_particles_to_grid();

    int n = 0;
    for (int di = 0; di <= 1; ++di) 
        for (int dj = 0; dj <= 1; ++dj)
            for (int dk = 0; dk <=1; ++dk) {
                for (int i = 0; i < 5; i+=4) {
                    REQUIRE(p.p2g_neighbour_node(i, n)==correct_ravelled_idx[i] + di*5*4 + dj*5 + dk);
                    std::array<int, 3> idx = p.grid_idx(i);
                    REQUIRE(p.p2g_neighbour_node_dx(i, n)==p.x(i)-((idx[0]+di)*dcell+mingrid[0]));
                    REQUIRE(p.p2g_neighbour_node_dy(i, n)==p.y(i)-((idx[1]+dj)*dcell+mingrid[1]));
                    REQUIRE(p.p2g_neighbour_node_dz(i, n)==p.z(i)-((idx[2]+dk)*dcell+mingrid[2]));
                    double w, dwdx, dwdy, dwdz;
                    knl.w_dwdx(
                        p.p2g_neighbour_node_dx(i, n), 
                        p.p2g_neighbour_node_dy(i, n), 
                        p.p2g_neighbour_node_dz(i, n),
                        w,
                        dwdx, dwdy, dwdz);
                    REQUIRE(p.p2g_neighbour_node_w(i, n)==w);
                    REQUIRE(p.p2g_neighbour_node_dwdx(i, n)==dwdx);
                    REQUIRE(p.p2g_neighbour_node_dwdy(i, n)==dwdy);
                    REQUIRE(p.p2g_neighbour_node_dwdz(i, n)==dwdz);
                }
                n++;
            }
}

TEST_CASE("Correct determination of grid node neighbours (radius=2)") {

    // initialize particle system with cublic spline kernel
    const std::array<double, 3> mingrid {-0.1, 0.05, 0.15}, maxgrid {0.19, 0.4, 0.6};
    GraMPM::grid<double> g(mingrid, maxgrid, dcell);
    GraMPM::kernel_cubic_bspline<double> knlc(dcell);
    GraMPM::particle_system<double> p(g, knlc);

    CHECK(p.background_grid.ngridx()==4);
    CHECK(p.background_grid.ngridy()==5);
    CHECK(p.background_grid.ngridz()==6);

    p.push_back(GraMPM::particle<double>(0.01, 0.16, 0.26, 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.));
    p.push_back(GraMPM::particle<double>(0.01, 0.3, 0.5, 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.));

    int correct_ravelled_idx[2] {37, 45};

    CHECK(p.ravelled_grid_idx(0)==correct_ravelled_idx[0]);
    CHECK(p.ravelled_grid_idx(1)==correct_ravelled_idx[1]);

    p.map_particles_to_grid();

    int n = 0;
    for (int di = -1; di <= 2; ++di) 
        for (int dj = -1; dj <=2; ++dj)
            for (int dk = -1; dk <=2; ++dk) {
                for (int i = 0; i < 2; ++i) {
                    REQUIRE(p.p2g_neighbour_node(i, n)==correct_ravelled_idx[i] + di*6*5 + dj*6 + dk);
                    std::array<int, 3> idx = p.grid_idx(i);
                    REQUIRE(p.p2g_neighbour_node_dx(i, n)==p.x(i)-((idx[0]+di)*dcell+mingrid[0]));
                    REQUIRE(p.p2g_neighbour_node_dy(i, n)==p.y(i)-((idx[1]+dj)*dcell+mingrid[1]));
                    REQUIRE(p.p2g_neighbour_node_dz(i, n)==p.z(i)-((idx[2]+dk)*dcell+mingrid[2]));
                    double w, dwdx, dwdy, dwdz;
                    knlc.w_dwdx(
                        p.p2g_neighbour_node_dx(i, n), 
                        p.p2g_neighbour_node_dy(i, n), 
                        p.p2g_neighbour_node_dz(i, n), 
                        w,
                        dwdx, dwdy, dwdz);
                    REQUIRE(p.p2g_neighbour_node_w(i, n)==w);
                    REQUIRE(p.p2g_neighbour_node_dwdx(i, n)==dwdx);
                    REQUIRE(p.p2g_neighbour_node_dwdy(i, n)==dwdy);
                    REQUIRE(p.p2g_neighbour_node_dwdz(i, n)==dwdz);
                }
                n++;
            }

}