#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "grampm.hpp"
#include <algorithm>
#include <grampm_kernels.hpp>
#include <array>

void generate_particles(GraMPM::particle_system<double> &p) {
    for (int i = 0; i < 10; ++i) 
        for (int j = 0; j < 20; ++j) 
            for (int k = 0; k < 30; ++k) {
                GraMPM::particle<double> p_i(0.1*(i+0.5), 0.1*(j+0.5), 0.1*(k+0.5), 10.*(i+j+k));
                p.push_back(p_i);
            }
}

TEST_CASE("Map particles masses to grid (linear bspline)") {
    const double dcell = 0.2;
    GraMPM::grid<double> g(0., 0., 0., 0.99, 1.99, 2.99, dcell);

    CHECK(g.ngridx()==6);
    CHECK(g.ngridy()==11);
    CHECK(g.ngridz()==16);
    GraMPM::kernel_linear_bspline<double> knl(dcell);

    GraMPM::particle_system<double> p(g, knl);

    generate_particles(p);

    p.map_particles_to_grid();
    p.map_mass_to_grid();

    // check total mass conservation
    double psum = 0., gsum = 0.;
    for (int i = 0; i < p.size(); ++i)
        psum += p.mass(i);
    for (int i = 0; i < p.background_grid.ncells(); ++i)
        gsum += p.background_grid.mass(i);

    REQUIRE(psum==gsum);

    // check a few nodal values
    REQUIRE(p.background_grid.mass(1, 1, 1)==360.);
    REQUIRE(std::round(p.background_grid.mass(3, 4, 5))==1800.);
    REQUIRE(std::round(p.background_grid.mass(5, 10, 15)*10.)==5625.);
}

TEST_CASE("Map particles masses to grid (cubic bspline)") {
    const double dcell = 0.2;
    GraMPM::grid<double> g(-0.2, -0.2, -0.2, 1.19, 2.19, 3.19, dcell);

    CHECK(g.ngridx()==8);
    CHECK(g.ngridy()==13);
    CHECK(g.ngridz()==18);
    GraMPM::kernel_cubic_bspline<double> knl(dcell);

    GraMPM::particle_system<double> p(g, knl);

    generate_particles(p);

    p.map_particles_to_grid();
    p.map_mass_to_grid();

    // check total mass conservation
    double psum = 0., gsum = 0.;
    for (int i = 0; i < p.size(); ++i)
        psum += p.mass(i);
    for (int i = 0; i < p.background_grid.ncells(); ++i)
        gsum += p.background_grid.mass(i);

    // should be correct to 14 sigfigs
    REQUIRE(psum*1e7==std::round(gsum*1e7));

    // check a few nodal values
    REQUIRE(std::round(p.background_grid.mass(2, 2, 2)*1e10)==3426422542996);
    REQUIRE(std::round(p.background_grid.mass(4, 5, 6))==1800.);
    REQUIRE(std::round(p.background_grid.mass(6, 11, 16)*1e5)==55609375);
}