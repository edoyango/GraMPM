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

    double psum = 0., gsum = 0.;
    for (int i = 0; i < p.size(); ++i)
        psum += p.mass(i);
    for (int i = 0; i < p.background_grid.ncells(); ++i)
        gsum += p.background_grid.mass(i);

    REQUIRE(psum==gsum);
}