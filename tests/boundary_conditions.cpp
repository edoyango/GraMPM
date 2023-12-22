#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include <catch.hpp>
#include <grampm.hpp>
#include <grampm_kernels.hpp>

void apply_lower_momentum(GraMPM::MPM_system<double> &self, const int &timestep, const double &dt) {
    for (int i = 0; i < self.g_ngridx; ++i)
        for (int j = 0; j < self.g_ngridy; ++j) {
            self.g_set_momentumx(i, j, 0, -1.);
            self.g_set_momentumy(i, j, 0, -1.);
            self.g_set_momentumz(i, j, 0, -1.);
        }
}

void apply_lower_force(GraMPM::MPM_system<double> &self, const int &timestep, const double &dt) {
    for (int i = 0; i < self.g_ngridx; ++i)
        for (int j = 0; j < self.g_ngridy; ++j) {
            self.g_set_forcex(i, j, 0, -2.);
            self.g_set_forcey(i, j, 0, -2.);
            self.g_set_forcez(i, j, 0, -2.);
        }
}

void apply_west_momentum(GraMPM::MPM_system<double> &self, const int &timestep, const double &dt) {
    for (int j = 0; j < self.g_ngridy; ++j)
        for (int k = 0; k < self.g_ngridz; ++k) {
            self.g_set_momentumx(0, j, k, -3.);
            self.g_set_momentumy(0, j, k, -3.);
            self.g_set_momentumz(0, j, k, -3.);
        }
}

void apply_west_force(GraMPM::MPM_system<double> &self, const int &timestep, const double &dt) {
    for (int j = 0; j < self.g_ngridy; ++j)
        for (int k = 0; k < self.g_ngridz; ++k) {
            self.g_set_forcex(0, j, k, -4.);
            self.g_set_forcey(0, j, k, -4.);
            self.g_set_forcez(0, j, k, -4.);
        }
}

TEST_CASE("Check that grid momentum and force boundary conditions have been applied correctly") {

    // check that the function works when initialized with the grid
    const double dcell = 0.2;
    // GraMPM::grid<double> g(0., 0., 0., 0.99, 1.99, 2.99, dcell, apply_lower_momentum, apply_lower_force);
    std::array<double, 3> bf, mingrid {0., 0., 0.}, maxgrid {0.99, 1.99, 2.99};
    GraMPM::kernel_cubic_bspline<double> knl(dcell);
    GraMPM::MPM_system<double> myMPM(bf, knl, mingrid, maxgrid, dcell);
    myMPM.g_set_momentum_boundary_function(apply_lower_momentum);
    myMPM.g_set_force_boundary_function(apply_lower_force);

    for (int i = 0; i < myMPM.g_size; ++i) {
        myMPM.g_momentumx[i] = 1.;
        myMPM.g_momentumy[i] = 1.;
        myMPM.g_momentumz[i] = 1.;
        myMPM.g_forcex[i] = 2.;
        myMPM.g_forcey[i] = 2.;
        myMPM.g_forcez[i] = 2.;
    }

    myMPM.g_apply_momentum_boundary_conditions(1, 0.);

    for (int i = 0; i < myMPM.g_ngridx; ++i)
        for (int j = 0; j < myMPM.g_ngridy; ++j)
            for (int k = 0; k < myMPM.g_ngridz; ++k) {
                if (k == 0) {
                    REQUIRE(myMPM.g_get_momentumx(i, j, k)==-1.);
                    REQUIRE(myMPM.g_get_momentumy(i, j, k)==-1.);
                    REQUIRE(myMPM.g_get_momentumz(i, j, k)==-1.);
                } else {
                    REQUIRE(myMPM.g_get_momentumx(i, j, k)==1.);
                    REQUIRE(myMPM.g_get_momentumy(i, j, k)==1.);
                    REQUIRE(myMPM.g_get_momentumz(i, j, k)==1.);
                }
            }

    myMPM.g_apply_force_boundary_conditions(1, 0.);

    for (int i = 0; i < myMPM.g_ngridx; ++i)
        for (int j = 0; j < myMPM.g_ngridy; ++j)
            for (int k = 0; k < myMPM.g_ngridz; ++k) {
                if (k == 0) {
                    REQUIRE(myMPM.g_get_forcex(i, j, k)==-2.);
                    REQUIRE(myMPM.g_get_forcey(i, j, k)==-2.);
                    REQUIRE(myMPM.g_get_forcez(i, j, k)==-2.);
                } else {
                    REQUIRE(myMPM.g_get_forcex(i, j, k)==2.);
                    REQUIRE(myMPM.g_get_forcey(i, j, k)==2.);
                    REQUIRE(myMPM.g_get_forcez(i, j, k)==2.);
                }
            }

    // setting new momentum boundary condition function
    myMPM.g_set_momentum_boundary_function(apply_west_momentum);
    myMPM.g_apply_momentum_boundary_conditions(1, 0.);

    for (int i = 0; i < myMPM.g_ngridx; ++i)
        for (int j = 0; j < myMPM.g_ngridy; ++j)
            for (int k = 0; k < myMPM.g_ngridz; ++k) {
                if (i == 0) {
                    REQUIRE(myMPM.g_get_momentumx(i, j, k)==-3.);
                    REQUIRE(myMPM.g_get_momentumy(i, j, k)==-3.);
                    REQUIRE(myMPM.g_get_momentumz(i, j, k)==-3.);
                } else if (k == 0) {
                    REQUIRE(myMPM.g_get_momentumx(i, j, k)==-1.);
                    REQUIRE(myMPM.g_get_momentumy(i, j, k)==-1.);
                    REQUIRE(myMPM.g_get_momentumz(i, j, k)==-1.);
                } else {
                    REQUIRE(myMPM.g_get_momentumx(i, j, k)==1.);
                    REQUIRE(myMPM.g_get_momentumy(i, j, k)==1.);
                    REQUIRE(myMPM.g_get_momentumz(i, j, k)==1.);
                }
            }

    // setting new force boundary condition function
    myMPM.g_set_force_boundary_function(apply_west_force);
    myMPM.g_apply_force_boundary_conditions(1, 0.);

    for (int i = 0; i < myMPM.g_ngridx; ++i)
        for (int j = 0; j < myMPM.g_ngridy; ++j)
            for (int k = 0; k < myMPM.g_ngridz; ++k) {
                if (i == 0) {
                    REQUIRE(myMPM.g_get_forcex(i, j, k)==-4.);
                    REQUIRE(myMPM.g_get_forcey(i, j, k)==-4.);
                    REQUIRE(myMPM.g_get_forcez(i, j, k)==-4.);
                } else if (k == 0) {
                    REQUIRE(myMPM.g_get_forcex(i, j, k)==-2.);
                    REQUIRE(myMPM.g_get_forcey(i, j, k)==-2.);
                    REQUIRE(myMPM.g_get_forcez(i, j, k)==-2.);
                } else {
                    REQUIRE(myMPM.g_get_forcex(i, j, k)==2.);
                    REQUIRE(myMPM.g_get_forcey(i, j, k)==2.);
                    REQUIRE(myMPM.g_get_forcez(i, j, k)==2.);
                }
            }

}