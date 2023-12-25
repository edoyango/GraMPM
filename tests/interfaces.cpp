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

GraMPM::kernel_linear_bspline<double> knl(dcell_in);

GraMPM::MPM_system<double> myMPM(5, bf, knl, mingridx_in, maxgridx_in, dcell_in);

TEST_CASE("grid intialized correctly", "[grid]") {

    REQUIRE(myMPM.g_cell_size()==dcell_in);

    // test array get interface
    std::array<double, 3> mingridx_out {myMPM.g_mingrid()}, maxgridx_out {myMPM.g_maxgrid()};

    REQUIRE(mingridx_out[0]==mingridx_in[0]);
    REQUIRE(mingridx_out[1]==mingridx_in[1]);
    REQUIRE(mingridx_out[2]==mingridx_in[2]);

    REQUIRE(maxgridx_out[0]==maxgridx_in[0]);
    REQUIRE(maxgridx_out[1]==maxgridx_in[1]);
    REQUIRE(maxgridx_out[2]==maxgridx_in[2]);

    std::array<size_t, 3> ngridx_out = myMPM.g_ngrid();
    REQUIRE(ngridx_out[0]==3);
    REQUIRE(ngridx_out[1]==4);
    REQUIRE(ngridx_out[2]==5);

    // test element-by-element get interface
    REQUIRE(myMPM.g_mingridx()==mingridx_in[0]);
    REQUIRE(myMPM.g_mingridy()==mingridx_in[1]);
    REQUIRE(myMPM.g_mingridz()==mingridx_in[2]);
    
    REQUIRE(myMPM.g_maxgridx()==maxgridx_in[0]);
    REQUIRE(myMPM.g_maxgridy()==maxgridx_in[1]);
    REQUIRE(myMPM.g_maxgridz()==maxgridx_in[2]);

    REQUIRE(myMPM.g_ngridx()==3);
    REQUIRE(myMPM.g_ngridy()==4);
    REQUIRE(myMPM.g_ngridz()==5);
    REQUIRE(myMPM.g_size()==60);

    for (size_t i = 0; i < myMPM.g_size(); ++i) {
        REQUIRE(myMPM.g_mass(i)==0.);
        REQUIRE(myMPM.g_momentumx(i)==0.);
        REQUIRE(myMPM.g_momentumy(i)==0.);
        REQUIRE(myMPM.g_momentumz(i)==0.);
        REQUIRE(myMPM.g_forcex(i)==0.);
        REQUIRE(myMPM.g_forcey(i)==0.);
        REQUIRE(myMPM.g_forcez(i)==0.);
    }

    for (size_t i = 0; i < myMPM.g_ngridx(); ++i) 
        for (size_t j = 0; j < myMPM.g_ngridy(); ++j)
            for (size_t k = 0; k < myMPM.g_ngridz(); ++k) {
                myMPM.g_mass(i, j, k) = (i+j+k)*1.;
                REQUIRE(myMPM.g_mass(i, j, k)==(i+j+k)*1.);
                myMPM.g_momentumx(i, j, k) = (i+j+k)*2.;
                REQUIRE(myMPM.g_momentumx(i, j, k)==(i+j+k)*2.);
                myMPM.g_momentumy(i, j, k) = (i+j+k)*3.;
                REQUIRE(myMPM.g_momentumy(i, j, k)==(i+j+k)*3.);
                myMPM.g_momentumz(i, j, k) = (i+j+k)*4.;
                REQUIRE(myMPM.g_momentumz(i, j, k)==(i+j+k)*4.);
                myMPM.g_forcex(i, j, k) = (i+j+k)*5.;
                REQUIRE(myMPM.g_forcex(i, j, k)==(i+j+k)*5.);
                myMPM.g_forcey(i, j, k) = (i+j+k)*6.;
                REQUIRE(myMPM.g_forcey(i, j, k)==(i+j+k)*6.);
                myMPM.g_forcez(i, j, k) = (i+j+k)*7.;
                REQUIRE(myMPM.g_forcez(i, j, k)==(i+j+k)*7.);
            }
}

TEST_CASE("Particle initialized correctly", "[grid]") {

    // test that body force vector is set to the correct values
    REQUIRE(myMPM.body_force(0) == 1.);
    REQUIRE(myMPM.body_force(1) == 2.);
    REQUIRE(myMPM.body_force(2) == 3.);

    // check that the myMPM instance has 5 myMPM, but all zeroed
    for (int i = 0; i < 5; ++i) {
        REQUIRE(myMPM.p_x(i)==0.);
        REQUIRE(myMPM.p_y(i)==0.);
        REQUIRE(myMPM.p_z(i)==0.);
        REQUIRE(myMPM.p_vx(i)==0.);
        REQUIRE(myMPM.p_vy(i)==0.);
        REQUIRE(myMPM.p_vz(i)==0.);
        REQUIRE(myMPM.p_ax(i)==0.);
        REQUIRE(myMPM.p_ay(i)==0.);
        REQUIRE(myMPM.p_az(i)==0.);
        REQUIRE(myMPM.p_dxdt(i)==0.);
        REQUIRE(myMPM.p_dydt(i)==0.);
        REQUIRE(myMPM.p_dzdt(i)==0.);
        REQUIRE(myMPM.p_mass(i)==0.);
        REQUIRE(myMPM.p_sigmaxx(i)==0.);
        REQUIRE(myMPM.p_sigmayy(i)==0.);
        REQUIRE(myMPM.p_sigmazz(i)==0.);
        REQUIRE(myMPM.p_sigmaxy(i)==0.);
        REQUIRE(myMPM.p_sigmaxz(i)==0.);
        REQUIRE(myMPM.p_sigmayz(i)==0.);
        REQUIRE(myMPM.p_strainratexx(i)==0.);
        REQUIRE(myMPM.p_strainrateyy(i)==0.);
        REQUIRE(myMPM.p_strainratezz(i)==0.);
        REQUIRE(myMPM.p_strainratexy(i)==0.);
        REQUIRE(myMPM.p_strainratexz(i)==0.);
        REQUIRE(myMPM.p_strainrateyz(i)==0.);
        REQUIRE(myMPM.p_spinratexy(i)==0.);
        REQUIRE(myMPM.p_spinratexz(i)==0.);
        REQUIRE(myMPM.p_spinrateyz(i)==0.);
    }

    REQUIRE(myMPM.p_size()==5);

    // check that aggregate interface works
    std::vector<GraMPM::particle<double>> pv;
    for (int i = 0; i < 5; ++i) {
        pv.push_back(
            GraMPM::particle<double>(i, 2.*i, 3.*i, 4.*i, 5.*i, 6.*i, 10.*i, 100.*i, -0.1*i, -0.2*i, -0.3*i, -0.4*i, 
                -0.5*i, -0.6*i, 7.*i, 8.*i, 9.*i, 10.*i, 11.*i, 12.*i, -0.7*i, -0.8*i, -0.9*i, -1.0*i, -1.1*i, -1.2*i,
                -1.3*i, -1.4*i, -1.5*i)
        );
    }
    GraMPM::MPM_system<double> myMPM2(pv, bf, knl, mingridx_in, maxgridx_in, dcell_in);

    for (int i = 0; i < 5; ++i) {
        REQUIRE(myMPM2.p_x(i)==i);
        REQUIRE(myMPM2.p_y(i)==2.*i);
        REQUIRE(myMPM2.p_z(i)==3.*i);
        REQUIRE(myMPM2.p_vx(i)==4.*i);
        REQUIRE(myMPM2.p_vy(i)==5.*i);
        REQUIRE(myMPM2.p_vz(i)==6.*i);
        REQUIRE(myMPM2.p_ax(i)==7.*i);
        REQUIRE(myMPM2.p_ay(i)==8.*i);
        REQUIRE(myMPM2.p_az(i)==9.*i);
        REQUIRE(myMPM2.p_dxdt(i)==10.*i);
        REQUIRE(myMPM2.p_dydt(i)==11.*i);
        REQUIRE(myMPM2.p_dzdt(i)==12.*i);
        REQUIRE(myMPM2.p_mass(i)==10.*i);
        REQUIRE(myMPM2.p_rho(i)==100.*i);
        REQUIRE(myMPM2.p_sigmaxx(i)==-0.1*i);
        REQUIRE(myMPM2.p_sigmayy(i)==-0.2*i);
        REQUIRE(myMPM2.p_sigmazz(i)==-0.3*i);
        REQUIRE(myMPM2.p_sigmaxy(i)==-0.4*i);
        REQUIRE(myMPM2.p_sigmaxz(i)==-0.5*i);
        REQUIRE(myMPM2.p_sigmayz(i)==-0.6*i);
        REQUIRE(myMPM2.p_strainratexx(i)==-0.7*i);
        REQUIRE(myMPM2.p_strainrateyy(i)==-0.8*i);
        REQUIRE(myMPM2.p_strainratezz(i)==-0.9*i);
        REQUIRE(myMPM2.p_strainratexy(i)==-1.0*i);
        REQUIRE(myMPM2.p_strainratexz(i)==-1.1*i);
        REQUIRE(myMPM2.p_strainrateyz(i)==-1.2*i);
        REQUIRE(myMPM2.p_spinratexy(i)==-1.3*i);
        REQUIRE(myMPM2.p_spinratexz(i)==-1.4*i);
        REQUIRE(myMPM2.p_spinrateyz(i)==-1.5*i);
    }

    // check aggregate getters
    for (int i = 0; i < 5; ++i) {
        GraMPM::particle<double> p = myMPM2.p_at(i);
        REQUIRE(p.x[0]==i);
        REQUIRE(p.x[1]==2.*i);
        REQUIRE(p.x[2]==3.*i);
        REQUIRE(p.v[0]==4.*i);
        REQUIRE(p.v[1]==5.*i);
        REQUIRE(p.v[2]==6.*i);
        REQUIRE(p.a[0]==7.*i);
        REQUIRE(p.a[1]==8.*i);
        REQUIRE(p.a[2]==9.*i);
        REQUIRE(p.dxdt[0]==10.*i);
        REQUIRE(p.dxdt[1]==11.*i);
        REQUIRE(p.dxdt[2]==12.*i);
        REQUIRE(p.mass==10.*i);
        REQUIRE(p.rho==100.*i);
        REQUIRE(p.sigma[0]==-0.1*i);
        REQUIRE(p.sigma[1]==-0.2*i);
        REQUIRE(p.sigma[2]==-0.3*i);
        REQUIRE(p.sigma[3]==-0.4*i);
        REQUIRE(p.sigma[4]==-0.5*i);
        REQUIRE(p.sigma[5]==-0.6*i);
        REQUIRE(p.strainrate[0]==-0.7*i);
        REQUIRE(p.strainrate[1]==-0.8*i);
        REQUIRE(p.strainrate[2]==-0.9*i);
        REQUIRE(p.strainrate[3]==-1.0*i);
        REQUIRE(p.strainrate[4]==-1.1*i);
        REQUIRE(p.strainrate[5]==-1.2*i);
        REQUIRE(p.spinrate[0]==-1.3*i);
        REQUIRE(p.spinrate[1]==-1.4*i);
        REQUIRE(p.spinrate[2]==-1.5*i);
    }
}

TEST_CASE("IO", "[myMPM]") {
    // currently incomplete

    myMPM.save_to_file("testfile", 1);

    GraMPM::MPM_system<double> myMPM3("testfile0000001", bf, knl, mingridx_in, maxgridx_in, dcell_in);

    REQUIRE(myMPM3.p_size()==5);

    // check that manual setter functions work
    for (int i = 0; i < 1; ++i) {
        REQUIRE(myMPM3.p_x(i)==-1.*i);
        REQUIRE(myMPM3.p_y(i)==-2.*i);
        REQUIRE(myMPM3.p_z(i)==-3.*i);
        REQUIRE(myMPM3.p_vx(i)==-4.*i);
        REQUIRE(myMPM3.p_vy(i)==-5.*i);
        REQUIRE(myMPM3.p_vz(i)==-6.*i);
        REQUIRE(myMPM3.p_ax(i)==-7.*i);
        REQUIRE(myMPM3.p_ay(i)==-8.*i);
        REQUIRE(myMPM3.p_az(i)==-9.*i);
        REQUIRE(myMPM3.p_dxdt(i)==-10.*i);
        REQUIRE(myMPM3.p_dydt(i)==-11.*i);
        REQUIRE(myMPM3.p_dzdt(i)==-12.*i);
        REQUIRE(myMPM3.p_mass(i)==30.*i);
        REQUIRE(myMPM3.p_rho(i)==40.*i);
        REQUIRE(myMPM3.p_sigmaxx(i)==-0.1*i);
        REQUIRE(myMPM3.p_sigmayy(i)==-0.2*i);
        REQUIRE(myMPM3.p_sigmazz(i)==-0.3*i);
        REQUIRE(myMPM3.p_sigmaxy(i)==-0.4*i);
        REQUIRE(myMPM3.p_sigmaxz(i)==-0.5*i);
        REQUIRE(myMPM3.p_sigmayz(i)==-0.6*i);
        REQUIRE(myMPM3.p_strainratexx(i)==-0.7*i);
        REQUIRE(myMPM3.p_strainrateyy(i)==-0.8*i);
        REQUIRE(myMPM3.p_strainratezz(i)==-0.9*i);
        REQUIRE(myMPM3.p_strainratexy(i)==-1.0*i);
        REQUIRE(myMPM3.p_strainratexz(i)==-1.1*i);
        REQUIRE(myMPM3.p_strainrateyz(i)==-1.2*i);
        REQUIRE(myMPM3.p_spinratexy(i)==-1.3*i);
        REQUIRE(myMPM3.p_spinratexz(i)==-1.4*i);
        REQUIRE(myMPM3.p_spinrateyz(i)==-1.5*i);
    }
}

TEST_CASE("Check clearing and resizing", "[myMPM]") {
    myMPM.p_clear();
    REQUIRE(myMPM.p_empty());
    myMPM.p_resize(3);
    REQUIRE(myMPM.p_size()==3);
    for (int i = 0; i < 3; ++i) {
        REQUIRE(myMPM.p_x(i)==0.);
        REQUIRE(myMPM.p_y(i)==0.);
        REQUIRE(myMPM.p_z(i)==0.);
        REQUIRE(myMPM.p_vx(i)==0.);
        REQUIRE(myMPM.p_vy(i)==0.);
        REQUIRE(myMPM.p_vz(i)==0.);
        REQUIRE(myMPM.p_ax(i)==0.);
        REQUIRE(myMPM.p_ay(i)==0.);
        REQUIRE(myMPM.p_az(i)==0.);
        REQUIRE(myMPM.p_dxdt(i)==0.);
        REQUIRE(myMPM.p_dydt(i)==0.);
        REQUIRE(myMPM.p_dzdt(i)==0.);
        REQUIRE(myMPM.p_mass(i)==0.);
        REQUIRE(myMPM.p_rho(i)==0.);
        REQUIRE(myMPM.p_sigmaxx(i)==0.);
        REQUIRE(myMPM.p_sigmayy(i)==0.);
        REQUIRE(myMPM.p_sigmazz(i)==0.);
        REQUIRE(myMPM.p_sigmaxy(i)==0.);
        REQUIRE(myMPM.p_sigmaxz(i)==0.);
        REQUIRE(myMPM.p_sigmayz(i)==0.);
        REQUIRE(myMPM.p_strainratexx(i)==0.);
        REQUIRE(myMPM.p_strainrateyy(i)==0.);
        REQUIRE(myMPM.p_strainratezz(i)==0.);
        REQUIRE(myMPM.p_strainratexy(i)==0.);
        REQUIRE(myMPM.p_strainratexz(i)==0.);
        REQUIRE(myMPM.p_strainrateyz(i)==0.);
        REQUIRE(myMPM.p_spinratexy(i)==0.);
        REQUIRE(myMPM.p_spinratexz(i)==0.);
        REQUIRE(myMPM.p_spinrateyz(i)==0.);
        REQUIRE(myMPM.p_grid_idx(i)==0);
    }
}