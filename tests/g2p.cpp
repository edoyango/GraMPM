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
                GraMPM::particle<double> p_i(0.1*(i+0.5), 0.1*(j+0.5), 0.1*(k+0.5), 
                                             -i, -j, -k,
                                             10.*(i+j+k+3.), 100.*(i+j+k+3.),
                                             0., 0., 0., 0., 0., 0.);
                p.push_back(p_i);
            }
}

TEST_CASE("Calculate particles' accelerations (linear bspline)") {
    
    const double dcell = 0.2;
    GraMPM::grid<double> g(0., 0., 0., 0.99, 1.99, 2.99, dcell);

    CHECK(g.ngridx()==6);
    CHECK(g.ngridy()==11);
    CHECK(g.ngridz()==16);
    GraMPM::kernel_linear_bspline<double> knl(dcell);

    GraMPM::particle_system<double> p(g, knl);

    generate_particles(p);

    // setup grid
    for (int i = 0 ; i < p.background_grid.ngridx(); ++i) 
        for (int j = 0; j < p.background_grid.ngridy(); ++j)
            for (int k = 0; k < p.background_grid.ngridz(); ++k) {
                p.background_grid.set_momentumx(i, j, k, 0.1*(i+j+k));
                p.background_grid.set_momentumy(i, j, k, 0.2*(i+j+k));
                p.background_grid.set_momentumz(i, j, k, 0.3*(i+j+k));
                p.background_grid.set_forcex(i, j, k, 0.4*(i+j+k));
                p.background_grid.set_forcey(i, j, k, 0.5*(i+j+k));
                p.background_grid.set_forcez(i, j, k, 0.6*(i+j+k));
            }

    p.map_particles_to_grid();
    p.map_mass_to_grid();
    p.map_acceleration_to_particles();

    // check conservation
    double psum = 0., gsum = 0.;
    for (int i = 0; i < p.size(); ++i)
        psum += p.ax(i)*p.mass(i);
    for (int i = 0; i < p.background_grid.ncells(); ++i)
        gsum += p.background_grid.forcex(i);

    REQUIRE(std::round(psum)==std::round(gsum));

    psum = 0., gsum = 0.;
    for (int i = 0; i < p.size(); ++i)
        psum += p.az(i)*p.mass(i);
    for (int i = 0; i < p.background_grid.ncells(); ++i)
        gsum += p.background_grid.forcez(i);

    REQUIRE(std::round(psum)==std::round(gsum));

    psum = 0., gsum = 0.;
    for (int i = 0; i < p.size(); ++i)
        psum += p.az(i)*p.mass(i);
    for (int i = 0; i < p.background_grid.ncells(); ++i)
        gsum += p.background_grid.forcez(i);

    REQUIRE(std::round(psum)==std::round(gsum));
}