#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "grampm.hpp"
#include <algorithm>
#include <grampm_kernels.hpp>
#include <array>

// helper function to generate particles
void generate_particles(GraMPM::MPM_system<double> &p) {
    for (int i = 0; i < 10; ++i) 
        for (int j = 0; j < 20; ++j) 
            for (int k = 0; k < 30; ++k) {
                GraMPM::particle<double> p_i(0.1*(i+0.5), 0.1*(j+0.5), 0.1*(k+0.5), 
                                             -i, -j, -k,
                                             10.*(i+j+k), 100.*(i+j+k+3.),
                                             0., 0., 0., 0., 0., 0.);
                p.p_push_back(p_i);
            }
}

TEST_CASE("Map particles masses to grid (linear bspline)") {
    const double dcell = 0.2;
    std::array<double, 3> bf, mingrid {0., 0., 0.}, maxgrid {0.99, 1.99, 2.99};

    GraMPM::kernel_linear_bspline<double> knl(dcell);

    GraMPM::MPM_system<double> p(bf, knl, mingrid, maxgrid, dcell);

    CHECK(p.g_ngridx()==6);
    CHECK(p.g_ngridy()==11);
    CHECK(p.g_ngridz()==16);

    generate_particles(p);

    p.map_particles_to_grid();
    p.map_p2g_mass();

    // check total mass conservation
    // sum particles' mass
    double psum = 0., gsum = 0.;
    for (size_t i = 0; i < p.p_size(); ++i)
        psum += p.p_mass(i);
    // sum grid's mass
    for (size_t i = 0; i < p.g_size(); ++i)
        gsum += p.g_mass(i);

    REQUIRE(psum==gsum);

    // check a few nodal values
    REQUIRE(p.g_mass(1, 1, 1)==360.);
    REQUIRE(std::round(p.g_mass(3, 4, 5))==1800.);
    REQUIRE(std::round(p.g_mass(5, 10, 15)*10.)==5625.);
}

TEST_CASE("Map particles masses to grid (cubic bspline)") {
    const double dcell = 0.2;
    std::array<double, 3> bf, mingrid {-0.2, -0.2, -0.2}, maxgrid {1.19, 2.19, 3.19};

    GraMPM::kernel_cubic_bspline<double> knl(dcell);

    GraMPM::MPM_system<double> p(bf, knl, mingrid, maxgrid, dcell);
    
    CHECK(p.g_ngridx()==8);
    CHECK(p.g_ngridy()==13);
    CHECK(p.g_ngridz()==18);

    generate_particles(p);

    p.map_particles_to_grid();
    p.map_p2g_mass();

    // check total mass conservation
    // sum particles' mass
    double psum = 0., gsum = 0.;
    for (size_t i = 0; i < p.p_size(); ++i)
        psum += p.p_mass(i);
    // sum grid's mass
    for (size_t i = 0; i < p.g_size(); ++i)
        gsum += p.g_mass(i);

    // should be correct to 14 sigfigs
    REQUIRE(psum*1e7==std::round(gsum*1e7));

    // check a few nodal values
    REQUIRE(std::round(p.g_mass(2, 2, 2)*1e10)==3426422542996);
    REQUIRE(std::round(p.g_mass(4, 5, 6))==1800.);
    REQUIRE(std::round(p.g_mass(6, 11, 16)*1e5)==55609375);
}

TEST_CASE("Map particles momentums to grid (linear bspline)") {
    const double dcell = 0.2;
    std::array<double, 3> bf, mingrid {0., 0., 0.}, maxgrid {0.99, 1.99, 2.99};

    GraMPM::kernel_linear_bspline<double> knl(dcell);

    GraMPM::MPM_system<double> p(bf, knl, mingrid, maxgrid, dcell);

    CHECK(p.g_ngridx()==6);
    CHECK(p.g_ngridy()==11);
    CHECK(p.g_ngridz()==16);

    generate_particles(p);

    p.map_particles_to_grid();
    p.map_p2g_momentum();

    // check momentum conservation
    // sum particles' momentum
    double psum[3] {0., 0., 0.}, gsum[3] {0., 0., 0.};
    for (size_t i = 0; i < p.p_size(); ++i) {
        psum[0] += p.p_mass(i)*p.p_vx(i);
        psum[1] += p.p_mass(i)*p.p_vy(i);
        psum[2] += p.p_mass(i)*p.p_vz(i);
    }

    // sum grid's momentum
    for (size_t i = 0; i < p.g_size(); ++i) {
        gsum[0] += p.g_momentumx(i);
        gsum[1] += p.g_momentumy(i);
        gsum[2] += p.g_momentumz(i);
    }

    // test
    REQUIRE(psum[0]==gsum[0]);
    REQUIRE(psum[1]==gsum[1]);
    REQUIRE(psum[2]==gsum[2]);

    // check a few nodal values
    REQUIRE(std::round(p.g_momentumx(1, 1, 1))==-600.);
    REQUIRE(std::round(p.g_momentumx(3, 4, 5))==-9960.);
    REQUIRE(std::round(p.g_momentumx(5, 10, 15)*100.)==-492375.);
    REQUIRE(std::round(p.g_momentumy(1, 1, 1))==-600.);
    REQUIRE(std::round(p.g_momentumy(3, 4, 5))==-13560.);
    REQUIRE(std::round(p.g_momentumy(5, 10, 15)*100.)==-1054875.);
    REQUIRE(std::round(p.g_momentumz(1, 1, 1))==-600.);
    REQUIRE(std::round(p.g_momentumz(3, 4, 5))==-17160);
    REQUIRE(std::round(p.g_momentumz(5, 10, 15)*100.)==-1617375.);
}

TEST_CASE("Map particles momentums to grid (cubic bspline)") {
    const double dcell = 0.2;
    std::array<double, 3> bf, mingrid {-0.2, -0.2, -0.2}, maxgrid {1.19, 2.19, 3.19};

    GraMPM::kernel_cubic_bspline<double> knl(dcell);

    GraMPM::MPM_system<double> p(bf, knl, mingrid, maxgrid, dcell);

    CHECK(p.g_ngridx()==8);
    CHECK(p.g_ngridy()==13);
    CHECK(p.g_ngridz()==18);

    generate_particles(p);

    p.map_particles_to_grid();
    p.map_p2g_momentum();

    // check conservation
    double psum[3] {0., 0., 0.}, gsum[3] {0., 0., 0.};
    for (size_t i = 0; i < p.p_size(); ++i) {
        psum[0] += p.p_mass(i)*p.p_vx(i);
        psum[1] += p.p_mass(i)*p.p_vy(i);
        psum[2] += p.p_mass(i)*p.p_vz(i);
    }
    for (size_t i = 0; i < p.g_size(); ++i) {
        gsum[0] += p.g_momentumx(i);
        gsum[1] += p.g_momentumy(i);
        gsum[2] += p.g_momentumz(i);
    }

    REQUIRE(psum[0]*1e6==std::round(gsum[0]*1e6));
    REQUIRE(psum[1]*1e6==std::round(gsum[1]*1e6));
    REQUIRE(psum[2]*1e6==std::round(gsum[2]*1e6));

    // check a few nodal values
    REQUIRE(std::round(p.g_momentumx(2, 2, 2)*1e6)==-627705941.);
    REQUIRE(std::round(p.g_momentumx(4, 5, 6)*1e6)==-10006666667.);
    REQUIRE(std::round(p.g_momentumx(6, 11, 16)*1e6)==-4751120334.);
    REQUIRE(std::round(p.g_momentumy(2, 2, 2)*1e6)==-627705941.);
    REQUIRE(std::round(p.g_momentumy(4, 5, 6)*1e6)==-13606666667.);
    REQUIRE(std::round(p.g_momentumy(6, 11, 16)*1e6)==-10312057834.);
    REQUIRE(std::round(p.g_momentumz(2, 2, 2)*1e6)==-627705941.);
    REQUIRE(std::round(p.g_momentumz(4, 5, 6)*1e6)==-17206666667.);
    REQUIRE(std::round(p.g_momentumz(6, 11, 16)*1e6)==-15872995334.);
}

TEST_CASE("Calculate force on grid (linear bspline)") {
    const double dcell = 0.2;
    std::array<double, 3> bf {1., 2., 3.}, mingrid {0., 0., 0.}, maxgrid {0.99, 1.99, 2.99};

    GraMPM::kernel_linear_bspline<double> knl(dcell);

    GraMPM::MPM_system<double> p(bf, knl, mingrid, maxgrid, dcell);

    CHECK(p.g_ngridx()==6);
    CHECK(p.g_ngridy()==11);
    CHECK(p.g_ngridz()==16);

    generate_particles(p);

    p.map_particles_to_grid();

    p.map_p2g_force();

    // check conservation
    double psum[3] {0., 0., 0.}, gsum[3] {0., 0., 0.};
    for (size_t i = 0; i < p.p_size(); ++i) {
        psum[0] += p.p_mass(i)*p.body_force(0);
        psum[1] += p.p_mass(i)*p.body_force(1);
        psum[2] += p.p_mass(i)*p.body_force(2);
    }
    for (size_t i = 0; i < p.g_size(); ++i) {
        gsum[0] += p.g_forcex(i);
        gsum[1] += p.g_forcey(i);
        gsum[2] += p.g_forcez(i);
    }

    REQUIRE(psum[0]==gsum[0]);
    REQUIRE(psum[1]==gsum[1]);
    REQUIRE(psum[2]==gsum[2]);

    // check a few nodal values
    REQUIRE(std::round(p.g_forcex(1, 1, 1))==360.);
    REQUIRE(std::round(p.g_forcex(3, 4, 5))==1800.);
    REQUIRE(std::round(p.g_forcex(5, 10, 15)*10.)==5625.);
    REQUIRE(std::round(p.g_forcey(1, 1, 1))==720.);
    REQUIRE(std::round(p.g_forcey(3, 4, 5))==3600.);
    REQUIRE(std::round(p.g_forcey(5, 10, 15)*10.)==11250.);
    REQUIRE(std::round(p.g_forcez(1, 1, 1))==1080.);
    REQUIRE(std::round(p.g_forcez(3, 4, 5))==5400.);
    REQUIRE(std::round(p.g_forcez(5, 10, 15)*10.)==16875.);

    // try with non-zero stresses
    for (size_t i = 0; i < p.p_size(); ++i) {
        p.p_sigmaxx(i) = p.p_x(i);
        p.p_sigmayy(i) = p.p_y(i);
        p.p_sigmazz(i) = p.p_z(i);
        p.p_sigmaxy(i) = p.p_x(i)-p.p_y(i);
        p.p_sigmaxz(i) = p.p_x(i)-p.p_z(i);
        p.p_sigmayz(i) = p.p_y(i)-p.p_z(i);
    }

    p.map_p2g_force();

    // check a few nodal values
    REQUIRE(std::round(p.g_forcex(1, 1, 1)*1e10)==3596170004735.);
    REQUIRE(std::round(p.g_forcex(3, 4, 5)*1e10)==17992941517835.);
    REQUIRE(std::round(p.g_forcex(5, 10, 15)*1e10)==5644457328379.);
    REQUIRE(std::round(p.g_forcey(1, 1, 1)*1e10)==7205551538826.);
    REQUIRE(std::round(p.g_forcey(3, 4, 5)*1e10)==36007203137483.);
    REQUIRE(std::round(p.g_forcey(5, 10, 15)*1e10)==11250948927821.);
    REQUIRE(std::round(p.g_forcez(1, 1, 1)*1e10)==10814933072917.);
    REQUIRE(std::round(p.g_forcez(3, 4, 5)*1e10)==54021315681200.);
    REQUIRE(std::round(p.g_forcez(5, 10, 15)*1e10)==16876423394507.);

}

TEST_CASE("Calculate force on grid (cubic bspline)") {
    const double dcell = 0.2;
    std::array<double, 3> bf {1., 2., 3.}, mingrid {-0.2, -0.2, -0.2}, maxgrid {1.19, 2.19, 3.19};

    GraMPM::kernel_cubic_bspline<double> knl(dcell);

    GraMPM::MPM_system<double> p(bf, knl, mingrid, maxgrid, dcell);

    generate_particles(p);
    
    CHECK(p.g_ngridx()==8);
    CHECK(p.g_ngridy()==13);
    CHECK(p.g_ngridz()==18);

    p.map_particles_to_grid();

    p.map_p2g_force();

    // check conservation
    double psum[3] {0., 0., 0.}, gsum[3] {0., 0., 0.};
    for (size_t i = 0; i < p.p_size(); ++i) {
        psum[0] += p.p_mass(i)*p.body_force(0);
        psum[1] += p.p_mass(i)*p.body_force(1);
        psum[2] += p.p_mass(i)*p.body_force(2);
    }
    for (size_t i = 0; i < p.g_size(); ++i) {
        gsum[0] += p.g_forcex(i);
        gsum[1] += p.g_forcey(i);
        gsum[2] += p.g_forcez(i);
    }

    REQUIRE(std::round(psum[0]*1e6)==std::round(gsum[0]*1e6));
    REQUIRE(std::round(psum[1]*1e6)==std::round(gsum[1]*1e6));
    REQUIRE(std::round(psum[2]*1e6)==std::round(gsum[2]*1e6));

    // check a few nodal values
    REQUIRE(std::round(p.g_forcex(2, 2, 2)*1e6)==342642254.);
    REQUIRE(std::round(p.g_forcex(4, 5, 6)*1e6)==1800000000.);
    REQUIRE(std::round(p.g_forcex(6, 11, 16)*1e6)==556093750.);
    REQUIRE(std::round(p.g_forcey(2, 2, 2)*1e6)==685284509.);
    REQUIRE(std::round(p.g_forcey(4, 5, 6)*1e6)==3600000000.);
    REQUIRE(std::round(p.g_forcey(6, 11, 16)*1e6)==1112187500.);
    REQUIRE(std::round(p.g_forcez(2, 2, 2)*1e6)==1027926763.);
    REQUIRE(std::round(p.g_forcez(4, 5, 6)*1e6)==5400000000.);
    REQUIRE(std::round(p.g_forcez(6, 11, 16)*1e6)==1668281250.);

    // try with non-zero stresses
    for (size_t i = 0; i < p.p_size(); ++i) {
        p.p_sigmaxx(i) = p.p_x(i);
        p.p_sigmayy(i) = p.p_y(i);
        p.p_sigmazz(i) = p.p_z(i);
        p.p_sigmaxy(i) = p.p_x(i)-p.p_y(i);
        p.p_sigmaxz(i) = p.p_x(i)-p.p_z(i);
        p.p_sigmayz(i) = p.p_y(i)-p.p_z(i);
    }

    p.map_p2g_force();

    // check a few nodal values
    REQUIRE(std::round(p.g_forcex(2, 2, 2)*1e6)==342415998.);
    REQUIRE(std::round(p.g_forcex(4, 5, 6)*1e6)==1799294307.);
    REQUIRE(std::round(p.g_forcex(6, 11, 16)*1e6)==557428540.);
    REQUIRE(std::round(p.g_forcey(2, 2, 2)*1e6)==685779305.);
    REQUIRE(std::round(p.g_forcey(4, 5, 6)*1e6)==3600719937.);
    REQUIRE(std::round(p.g_forcey(6, 11, 16)*1e6)==1112282300.);
    REQUIRE(std::round(p.g_forcez(2, 2, 2)*1e6)==1029142612.);
    REQUIRE(std::round(p.g_forcez(4, 5, 6)*1e6)==5402130522.);
    REQUIRE(std::round(p.g_forcez(6, 11, 16)*1e6)==1668440054.);

}