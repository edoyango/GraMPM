#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "grampm.hpp"
#include <algorithm>
#include <grampm_kernels.hpp>
#include <array>

void generate_particles(GraMPM::MPM_system<double> &p) {
    for (int i = 0; i < 10; ++i) 
        for (int j = 0; j < 20; ++j) 
            for (int k = 0; k < 30; ++k) {
                GraMPM::particle<double> p_i(0.1*(i+0.5), 0.1*(j+0.5), 0.1*(k+0.5), 
                                             -i, -j, -k,
                                             10.*(i+j+k+3.), 100.*(i+j+k+3.),
                                             0., 0., 0., 0., 0., 0.);
                p.p_push_back(p_i);
            }
}

TEST_CASE("Calculate particles' accelerations (linear bspline)") {
    
    const double dcell = 0.2;

    std::array<double, 3> bf {0.}, mingrid {0.}, maxgrid {0.99, 1.99, 2.99};

    GraMPM::kernel_linear_bspline<double> knl(dcell);

    GraMPM::MPM_system<double> p(bf, knl, mingrid, maxgrid, dcell);

    CHECK(p.g_ngridx()==6);
    CHECK(p.g_ngridy()==11);
    CHECK(p.g_ngridz()==16);

    generate_particles(p);

    // setup grid
    for (int i = 0 ; i < p.g_ngridx(); ++i) 
        for (int j = 0; j < p.g_ngridy(); ++j)
            for (int k = 0; k < p.g_ngridz(); ++k) {
                p.g_momentumx(i, j, k) = 0.1*(i+j+k);
                p.g_momentumy(i, j, k) = 0.2*(i+j+k);
                p.g_momentumz(i, j, k) = 0.3*(i+j+k);
                p.g_forcex(i, j, k) = 0.4*(i+j+k);
                p.g_forcey(i, j, k) = 0.5*(i+j+k);
                p.g_forcez(i, j, k) = 0.6*(i+j+k);
            }

    p.map_particles_to_grid();
    p.map_p2g_mass();
    p.map_g2p_acceleration();

    // check conservation
    // sum particles' force (m*a) and momentum (dxdt*m)
    double psum[6] {0., 0., 0., 0., 0., 0.}, gsum[6] {0., 0., 0., 0., 0., 0.};
    for (int i = 0; i < p.p_size(); ++i) {
        psum[0] += p.p_ax(i)*p.p_mass(i);
        psum[1] += p.p_ay(i)*p.p_mass(i);
        psum[2] += p.p_az(i)*p.p_mass(i);
        psum[3] += p.p_dxdt(i)*p.p_mass(i);
        psum[4] += p.p_dydt(i)*p.p_mass(i);
        psum[5] += p.p_dzdt(i)*p.p_mass(i);
    }
    for (int i = 0; i < p.g_size(); ++i) {
        gsum[0] += p.g_forcex(i);
        gsum[1] += p.g_forcey(i);
        gsum[2] += p.g_forcez(i);
        gsum[3] += p.g_momentumx(i);
        gsum[4] += p.g_momentumy(i);
        gsum[5] += p.g_momentumz(i);
    }

    // test
    for (int i = 0; i < 6; ++i)
        REQUIRE(std::round(psum[i]) == std::round(gsum[i]));
}

TEST_CASE("Calculate particles' accelerations (cubic bspline)") {
    
    const double dcell = 0.2;
    // GraMPM::grid<double> g(-0.2, -0.2, -0.2, 1.19, 2.19, 3.19, dcell);
    std::array<double, 3> bf {0.}, mingrid {-0.2, -0.2, -0.2}, maxgrid {1.19, 2.19, 3.19};

    GraMPM::kernel_cubic_bspline<double> knl(dcell);

    GraMPM::MPM_system<double> p(bf, knl, mingrid, maxgrid, dcell);
    CHECK(p.g_ngridx()==8);
    CHECK(p.g_ngridy()==13);
    CHECK(p.g_ngridz()==18);

    generate_particles(p);

    // setup grid
    for (int i = 0 ; i < p.g_ngridx(); ++i) 
        for (int j = 0; j < p.g_ngridy(); ++j)
            for (int k = 0; k < p.g_ngridz(); ++k) {
                p.g_momentumx(i, j, k) = 0.1*(i+j+k);
                p.g_momentumy(i, j, k) = 0.2*(i+j+k);
                p.g_momentumz(i, j, k) = 0.3*(i+j+k);
                p.g_forcex(i, j, k) = 0.4*(i+j+k);
                p.g_forcey(i, j, k) = 0.5*(i+j+k);
                p.g_forcez(i, j, k) = 0.6*(i+j+k);
            }

    p.map_particles_to_grid();
    p.map_p2g_mass();
    p.map_g2p_acceleration();

    // check conservation
    // sum particles' force (m*a) and momentum (dxdt*m)
    double psum[6] {0., 0., 0., 0., 0., 0.}, gsum[6] {0., 0., 0., 0., 0., 0.};
    for (int i = 0; i < p.p_size(); ++i) {
        psum[0] += p.p_ax(i)*p.p_mass(i);
        psum[1] += p.p_ay(i)*p.p_mass(i);
        psum[2] += p.p_az(i)*p.p_mass(i);
        psum[3] += p.p_dxdt(i)*p.p_mass(i);
        psum[4] += p.p_dydt(i)*p.p_mass(i);
        psum[5] += p.p_dzdt(i)*p.p_mass(i);
    }
    for (int i = 0; i < p.g_size(); ++i) {
        gsum[0] += p.g_forcex(i);
        gsum[1] += p.g_forcey(i);
        gsum[2] += p.g_forcez(i);
        gsum[3] += p.g_momentumx(i);
        gsum[4] += p.g_momentumy(i);
        gsum[5] += p.g_momentumz(i);
    }

    // test
    for (int i = 0; i < 6; ++i)
        REQUIRE(std::round(psum[i]) == std::round(gsum[i]));
}

TEST_CASE("Calculate particles' strain/spin rates (linear bspline)") {
    
    const double dcell = 0.2;
    // GraMPM::grid<double> g(0., 0., 0., 0.99, 1.99, 2.99, dcell);
    std::array<double, 3> bf {0.}, mingrid {0., 0., 0.}, maxgrid {0.99, 1.99, 2.99};

    GraMPM::kernel_linear_bspline<double> knl(dcell);

    GraMPM::MPM_system<double> p(bf, knl, mingrid, maxgrid, dcell);

    CHECK(p.g_ngridx()==6);
    CHECK(p.g_ngridy()==11);
    CHECK(p.g_ngridz()==16);

    generate_particles(p);

    // setup grid
    for (int i = 0; i < p.g_ngridx(); ++i)
        for (int j = 0; j < p.g_ngridy(); ++j) 
            for (int k = 0; k < p.g_ngridz(); ++k) {
                const double x = i*0.2 + p.g_mingridx();
                const double y = j*0.2 + p.g_mingridy();
                const double z = k*0.2 + p.g_mingridz();
                p.g_momentumx(i, j, k) = 0.1*(x-y-z);
                p.g_momentumy(i, j, k) = 0.2*(y-x-z);
                p.g_momentumz(i, j, k) = 0.3*(z-x-y);
                p.g_mass(i, j, k) = 1.;
            }
    
    p.map_particles_to_grid();
    p.map_g2p_strainrate();

    REQUIRE(std::round(p.p_strainratexx(0)*100.)==10.);
    REQUIRE(std::round(p.p_strainrateyy(0)*100.)==20.);
    REQUIRE(std::round(p.p_strainratezz(0)*100.)==30.);
    REQUIRE(std::round(p.p_strainratexy(0)*100.)==-15.);
    REQUIRE(std::round(p.p_strainratexz(0)*100.)==-20.);
    REQUIRE(std::round(p.p_strainrateyz(0)*100.)==-25.);
    REQUIRE(std::round(p.p_spinratexy(0)*100.)==5.);
    REQUIRE(std::round(p.p_spinratexz(0)*100.)==10.);
    REQUIRE(std::round(p.p_spinrateyz(0)*100.)==5.);

    REQUIRE(std::round(p.p_strainratexx(p.p_size()-1)*100.)==10.);
    REQUIRE(std::round(p.p_strainrateyy(p.p_size()-1)*100.)==20.);
    REQUIRE(std::round(p.p_strainratezz(p.p_size()-1)*100.)==30.);
    REQUIRE(std::round(p.p_strainratexy(p.p_size()-1)*100.)==-15.);
    REQUIRE(std::round(p.p_strainratexz(p.p_size()-1)*100.)==-20.);
    REQUIRE(std::round(p.p_strainrateyz(p.p_size()-1)*100.)==-25.);
    REQUIRE(std::round(p.p_spinratexy(p.p_size()-1)*100.)==5.);
    REQUIRE(std::round(p.p_spinratexz(p.p_size()-1)*100.)==10.);
    REQUIRE(std::round(p.p_spinrateyz(p.p_size()-1)*100.)==5.);
}

TEST_CASE("Calculate particles' strain/spin rates (cubic bspline)") {
    
    const double dcell = 0.2;
    std::array<double, 3> bf {0.}, mingrid {-0.2, -0.2, -0.2}, maxgrid {1.19, 2.19, 3.19};

    GraMPM::kernel_cubic_bspline<double> knl(dcell);

    GraMPM::MPM_system<double> p(bf, knl, mingrid, maxgrid, dcell);

    CHECK(p.g_ngridx()==8);
    CHECK(p.g_ngridy()==13);
    CHECK(p.g_ngridz()==18);

    generate_particles(p);

    // setup grid
    for (int i = 0; i < p.g_ngridx(); ++i)
        for (int j = 0; j < p.g_ngridy(); ++j) 
            for (int k = 0; k < p.g_ngridz(); ++k) {
                const double x = i*0.2 + p.g_mingridx();
                const double y = j*0.2 + p.g_mingridy();
                const double z = k*0.2 + p.g_mingridz();
                p.g_momentumx(i, j, k) = 0.1*(x-y-z);
                p.g_momentumy(i, j, k) = 0.2*(y-x-z);
                p.g_momentumz(i, j, k) = 0.3*(z-x-y);
                p.g_mass(i, j, k) = 1.;
            }
    
    p.map_particles_to_grid();
    p.map_g2p_strainrate();

    REQUIRE(std::round(p.p_strainratexx(0)*100.)==10.);
    REQUIRE(std::round(p.p_strainrateyy(0)*100.)==20.);
    REQUIRE(std::round(p.p_strainratezz(0)*100.)==30.);
    REQUIRE(std::round(p.p_strainratexy(0)*100.)==-15.);
    REQUIRE(std::round(p.p_strainratexz(0)*100.)==-20.);
    REQUIRE(std::round(p.p_strainrateyz(0)*100.)==-25.);
    REQUIRE(std::round(p.p_spinratexy(0)*100.)==5.);
    REQUIRE(std::round(p.p_spinratexz(0)*100.)==10.);
    REQUIRE(std::round(p.p_spinrateyz(0)*100.)==5.);

    REQUIRE(std::round(p.p_strainratexx(p.p_size()-1)*100.)==10.);
    REQUIRE(std::round(p.p_strainrateyy(p.p_size()-1)*100.)==20.);
    REQUIRE(std::round(p.p_strainratezz(p.p_size()-1)*100.)==30.);
    REQUIRE(std::round(p.p_strainratexy(p.p_size()-1)*100.)==-15.);
    REQUIRE(std::round(p.p_strainratexz(p.p_size()-1)*100.)==-20.);
    REQUIRE(std::round(p.p_strainrateyz(p.p_size()-1)*100.)==-25.);
    REQUIRE(std::round(p.p_spinratexy(p.p_size()-1)*100.)==5.);
    REQUIRE(std::round(p.p_spinratexz(p.p_size()-1)*100.)==10.);
    REQUIRE(std::round(p.p_spinrateyz(p.p_size()-1)*100.)==5.);
}