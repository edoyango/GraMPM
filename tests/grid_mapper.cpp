// Ensuring that functions to determine grid neighbours do their job

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "grampm.hpp"
#include <algorithm>
#include <grampm_kernels.hpp>
#include <array>

const double dcell = 0.1;
const std::array<double, 3> mingrid {-0.1, 0.05, 0.15}, maxgrid {0.1, 0.3, 0.5};
GraMPM::Grid<double> g(mingrid, maxgrid, dcell);
GraMPM::kernel_linear_bspline<double> knl(dcell);
GraMPM::ParticleSystem<double> p(g, knl);


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
        GraMPM::particle<double> pi(i*dx+mingrid[0], i*dy+mingrid[1], i*dz+mingrid[2], 1.);
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


const int correct_neighbours0[8] {correct_ravelled_idx[0], 
                                  correct_ravelled_idx[0]+1, 
                                  correct_ravelled_idx[0]+5, 
                                  correct_ravelled_idx[0]+6,
                                  correct_ravelled_idx[0]+20,
                                  correct_ravelled_idx[0]+21,
                                  correct_ravelled_idx[0]+25, 
                                  correct_ravelled_idx[0]+26
                                 };
const double correct_dx0[8] {0., 0., 0., 0., -0.1, -0.1, -0.1, -0.1};
const double correct_dy0[8] {0., 0., -0.1, -0.1, 0., 0., -0.1, -0.1};
const double correct_dz0[8] {0., -0.1, 0., -0.1, 0., 0.-0.1, 0., -0.1};
const int correct_neighbours4[8] {correct_ravelled_idx[4], 
                                  correct_ravelled_idx[4]+1, 
                                  correct_ravelled_idx[4]+5, 
                                  correct_ravelled_idx[4]+6,
                                  correct_ravelled_idx[4]+20,
                                  correct_ravelled_idx[4]+21,
                                  correct_ravelled_idx[4]+25, 
                                  correct_ravelled_idx[4]+26
                                 };
const double correct_dx4[8] {0., 0., 0., 0., -0.1, -0.1, -0.1, -0.1};
const double correct_dy4[8] {0.05, 0.05, -0.05, -0.05, 0.05, 0.05, -0.05, -0.05};
const double correct_dz4[8] {0.05, -0.05, 0.05, -0.05, 0.05, -0.05, 0.05, -0.05};
TEST_CASE("Correct determination of grid node neighbours", "[p]") {

    p.get_neighbour_nodes();

    // for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 8; ++j) {
        REQUIRE(p.p2g_neighbour_node(0, j)==correct_neighbours0[j]);
        REQUIRE(std::abs(p.p2g_neighbour_node_dx(0, j)-correct_dx0[j])<1e-14);
        REQUIRE(std::abs(p.p2g_neighbour_node_dy(0, j)-correct_dy0[j])<1e-14);
        REQUIRE(std::abs(p.p2g_neighbour_node_dz(0, j)-correct_dz0[j])<1e-14);
    }
    for (int j = 0; j < 8; ++j) {
        REQUIRE(p.p2g_neighbour_node(4, j)==correct_neighbours4[j]);
        REQUIRE(std::abs(p.p2g_neighbour_node_dx(4, j)-correct_dx4[j])<1e-14);
        REQUIRE(std::abs(p.p2g_neighbour_node_dy(4, j)-correct_dy4[j])<1e-14);
        REQUIRE(std::abs(p.p2g_neighbour_node_dz(4, j)-correct_dz4[j])<1e-14);
    }
    // }
}