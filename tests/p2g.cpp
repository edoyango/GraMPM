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
                                             10.*(i+j+k), 100.*(i+j+k+3.),
                                             0., 0., 0., 0., 0., 0.);
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

TEST_CASE("Map particles momentums to grid (linear bspline)") {
    const double dcell = 0.2;
    GraMPM::grid<double> g(0., 0., 0., 0.99, 1.99, 2.99, dcell);

    CHECK(g.ngridx()==6);
    CHECK(g.ngridy()==11);
    CHECK(g.ngridz()==16);
    GraMPM::kernel_linear_bspline<double> knl(dcell);

    GraMPM::particle_system<double> p(g, knl);

    generate_particles(p);

    p.map_particles_to_grid();
    p.map_momentum_to_grid();

    double psum = 0., gsum = 0.;
    for (int i = 0; i < p.size(); ++i)
        psum += p.mass(i)*p.vx(i);
    for (int i = 0; i < p.background_grid.ncells(); ++i)
        gsum += p.background_grid.momentumx(i);

    REQUIRE(psum==gsum);

    psum = 0., gsum = 0.;
    for (int i = 0; i < p.size(); ++i)
        psum += p.mass(i)*p.vy(i);
    for (int i = 0; i < p.background_grid.ncells(); ++i)
        gsum += p.background_grid.momentumy(i);

    REQUIRE(psum==gsum);

    psum = 0., gsum = 0.;
    for (int i = 0; i < p.size(); ++i)
        psum += p.mass(i)*p.vz(i);
    for (int i = 0; i < p.background_grid.ncells(); ++i)
        gsum += p.background_grid.momentumz(i);

    REQUIRE(psum==gsum);

    // check a few nodal values
    REQUIRE(std::round(p.background_grid.momentumx(1, 1, 1))==-600.);
    REQUIRE(std::round(p.background_grid.momentumx(3, 4, 5))==-9960.);
    REQUIRE(std::round(p.background_grid.momentumx(5, 10, 15)*100.)==-492375.);
    REQUIRE(std::round(p.background_grid.momentumy(1, 1, 1))==-600.);
    REQUIRE(std::round(p.background_grid.momentumy(3, 4, 5))==-13560.);
    REQUIRE(std::round(p.background_grid.momentumy(5, 10, 15)*100.)==-1054875.);
    REQUIRE(std::round(p.background_grid.momentumz(1, 1, 1))==-600.);
    REQUIRE(std::round(p.background_grid.momentumz(3, 4, 5))==-17160);
    REQUIRE(std::round(p.background_grid.momentumz(5, 10, 15)*100.)==-1617375.);
}

TEST_CASE("Calculate force on grid (linear bspline)") {
    const double dcell = 0.2;
    const std::array<double, 3> bf {1., 2., 3.};
    GraMPM::grid<double> g(0., 0., 0., 0.99, 1.99, 2.99, dcell);

    CHECK(g.ngridx()==6);
    CHECK(g.ngridy()==11);
    CHECK(g.ngridz()==16);
    GraMPM::kernel_linear_bspline<double> knl(dcell);

    GraMPM::particle_system<double> p(g, knl);

    generate_particles(p);
    
    p.set_body_force(bf);

    p.map_particles_to_grid();

    p.map_force_to_grid();

    double psum = 0., gsum = 0.;
    for (int i = 0; i < p.size(); ++i)
        psum += p.mass(i)*p.body_force(0);
    for (int i = 0; i < p.background_grid.ncells(); ++i)
        gsum += p.background_grid.forcex(i);

    REQUIRE(psum==gsum);

    psum = 0., gsum = 0.;
    for (int i = 0; i < p.size(); ++i)
        psum += p.mass(i)*p.body_force(1);
    for (int i = 0; i < p.background_grid.ncells(); ++i)
        gsum += p.background_grid.forcey(i);

    REQUIRE(psum==gsum);

    psum = 0., gsum = 0.;
    for (int i = 0; i < p.size(); ++i)
        psum += p.mass(i)*p.body_force(2);
    for (int i = 0; i < p.background_grid.ncells(); ++i)
        gsum += p.background_grid.forcez(i);

    REQUIRE(psum==gsum);

    // check a few nodal values
    REQUIRE(std::round(p.background_grid.forcex(1, 1, 1))==360.);
    REQUIRE(std::round(p.background_grid.forcex(3, 4, 5))==1800.);
    REQUIRE(std::round(p.background_grid.forcex(5, 10, 15)*10.)==5625.);
    REQUIRE(std::round(p.background_grid.forcey(1, 1, 1))==720.);
    REQUIRE(std::round(p.background_grid.forcey(3, 4, 5))==3600.);
    REQUIRE(std::round(p.background_grid.forcey(5, 10, 15)*10.)==11250.);
    REQUIRE(std::round(p.background_grid.forcez(1, 1, 1))==1080.);
    REQUIRE(std::round(p.background_grid.forcez(3, 4, 5))==5400.);
    REQUIRE(std::round(p.background_grid.forcez(5, 10, 15)*10.)==16875.);

    // try with non-zero stresses
    for (int i = 0; i < p.size(); ++i) {
        p.set_sigmaxx(i, p.x(i));
        p.set_sigmayy(i, p.y(i));
        p.set_sigmazz(i, p.z(i));
        p.set_sigmaxy(i, p.x(i)-p.y(i));
        p.set_sigmaxz(i, p.x(i)-p.z(i));
        p.set_sigmayz(i, p.y(i)-p.z(i));
    }

    p.map_force_to_grid();

    // check a few nodal values
    REQUIRE(std::round(p.background_grid.forcex(1, 1, 1)*1e10)==3596170004735.);
    REQUIRE(std::round(p.background_grid.forcex(3, 4, 5)*1e10)==17992941517835.);
    REQUIRE(std::round(p.background_grid.forcex(5, 10, 15)*1e10)==5644457328379.);
    REQUIRE(std::round(p.background_grid.forcey(1, 1, 1)*1e10)==7205551538826.);
    REQUIRE(std::round(p.background_grid.forcey(3, 4, 5)*1e10)==36007203137483.);
    REQUIRE(std::round(p.background_grid.forcey(5, 10, 15)*1e10)==11250948927821.);
    REQUIRE(std::round(p.background_grid.forcez(1, 1, 1)*1e10)==10814933072917.);
    REQUIRE(std::round(p.background_grid.forcez(3, 4, 5)*1e10)==54021315681200.);
    REQUIRE(std::round(p.background_grid.forcez(5, 10, 15)*1e10)==16876423394507.);

}