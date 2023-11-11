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
    // sum particles' force (m*a) and momentum (dxdt*m)
    double psum[6] {0., 0., 0., 0., 0., 0.}, gsum[6] {0., 0., 0., 0., 0., 0.};
    for (int i = 0; i < p.size(); ++i) {
        psum[0] += p.ax(i)*p.mass(i);
        psum[1] += p.ay(i)*p.mass(i);
        psum[2] += p.az(i)*p.mass(i);
        psum[3] += p.dxdt(i)*p.mass(i);
        psum[4] += p.dydt(i)*p.mass(i);
        psum[5] += p.dzdt(i)*p.mass(i);
    }
    for (int i = 0; i < p.background_grid.ncells(); ++i) {
        gsum[0] += p.background_grid.forcex(i);
        gsum[1] += p.background_grid.forcey(i);
        gsum[2] += p.background_grid.forcez(i);
        gsum[3] += p.background_grid.momentumx(i);
        gsum[4] += p.background_grid.momentumy(i);
        gsum[5] += p.background_grid.momentumz(i);
    }

    // test
    for (int i = 0; i < 6; ++i)
        REQUIRE(std::round(psum[i]) == std::round(gsum[i]));
}

TEST_CASE("Calculate particles' strain/spin rates (linear bspline)") {
    
    const double dcell = 0.2;
    GraMPM::grid<double> g(0., 0., 0., 0.99, 1.99, 2.99, dcell);

    CHECK(g.ngridx()==6);
    CHECK(g.ngridy()==11);
    CHECK(g.ngridz()==16);
    GraMPM::kernel_linear_bspline<double> knl(dcell);

    GraMPM::particle_system<double> p(g, knl);

    generate_particles(p);

    // setup grid
    for (int i = 0; i < p.background_grid.ngridx(); ++i)
        for (int j = 0; j < p.background_grid.ngridy(); ++j) 
            for (int k = 0; k < p.background_grid.ngridz(); ++k) {
                const double x = i*0.2 + p.background_grid.mingridx();
                const double y = j*0.2 + p.background_grid.mingridy();
                const double z = k*0.2 + p.background_grid.mingridz();
                p.background_grid.set_momentumx(i, j, k, 0.1*(x-y-z));
                p.background_grid.set_momentumy(i, j, k, 0.2*(y-x-z));
                p.background_grid.set_momentumz(i, j, k, 0.3*(z-x-y));
                p.background_grid.set_mass(i, j, k, 1.);
            }
    
    p.map_particles_to_grid();
    p.map_strainrate_to_particles();

    REQUIRE(std::round(p.strainratexx(0)*100.)==10.);
    REQUIRE(std::round(p.strainrateyy(0)*100.)==20.);
    REQUIRE(std::round(p.strainratezz(0)*100.)==30.);
    REQUIRE(std::round(p.strainratexy(0)*100.)==-15.);
    REQUIRE(std::round(p.strainratexz(0)*100.)==-20.);
    REQUIRE(std::round(p.strainrateyz(0)*100.)==-25.);
    REQUIRE(std::round(p.spinratexy(0)*100.)==5.);
    REQUIRE(std::round(p.spinratexz(0)*100.)==10.);
    REQUIRE(std::round(p.spinrateyz(0)*100.)==5.);

    REQUIRE(std::round(p.strainratexx(p.size()-1)*100.)==10.);
    REQUIRE(std::round(p.strainrateyy(p.size()-1)*100.)==20.);
    REQUIRE(std::round(p.strainratezz(p.size()-1)*100.)==30.);
    REQUIRE(std::round(p.strainratexy(p.size()-1)*100.)==-15.);
    REQUIRE(std::round(p.strainratexz(p.size()-1)*100.)==-20.);
    REQUIRE(std::round(p.strainrateyz(p.size()-1)*100.)==-25.);
    REQUIRE(std::round(p.spinratexy(p.size()-1)*100.)==5.);
    REQUIRE(std::round(p.spinratexz(p.size()-1)*100.)==10.);
    REQUIRE(std::round(p.spinrateyz(p.size()-1)*100.)==5.);
}