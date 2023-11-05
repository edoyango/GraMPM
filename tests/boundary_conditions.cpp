#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include <catch.hpp>
#include <grampm.hpp>
#include <grampm_kernels.hpp>

void apply_lower_momentum(GraMPM::grid<double> &self, const int &timestep, const double &dt) {
    for (int i = 0; i < self.ngridx(); ++i)
        for (int j = 0; j < self.ngridy(); ++j) {
            self.set_momentumx(i, j, 0, -1.);
            self.set_momentumy(i, j, 0, -1.);
            self.set_momentumz(i, j, 0, -1.);
        }
}

void apply_lower_force(GraMPM::grid<double> &self, const int &timestep, const double &dt) {
    for (int i = 0; i < self.ngridx(); ++i)
        for (int j = 0; j < self.ngridy(); ++j) {
            self.set_forcex(i, j, 0, -2.);
            self.set_forcey(i, j, 0, -2.);
            self.set_forcez(i, j, 0, -2.);
        }
}

void apply_west_momentum(GraMPM::grid<double> &self, const int &timestep, const double &dt) {
    for (int j = 0; j < self.ngridy(); ++j)
        for (int k = 0; k < self.ngridz(); ++k) {
            self.set_momentumx(0, j, k, -3.);
            self.set_momentumy(0, j, k, -3.);
            self.set_momentumz(0, j, k, -3.);
        }
}

void apply_west_force(GraMPM::grid<double> &self, const int &timestep, const double &dt) {
    for (int j = 0; j < self.ngridy(); ++j)
        for (int k = 0; k < self.ngridz(); ++k) {
            self.set_forcex(0, j, k, -4.);
            self.set_forcey(0, j, k, -4.);
            self.set_forcez(0, j, k, -4.);
        }
}

TEST_CASE("Check that grid momentum and force boundary conditions have been applied correctly") {

    // check that the function works when initialized with the grid
    const double dcell = 0.2;
    GraMPM::grid<double> g(0., 0., 0., 0.99, 1.99, 2.99, dcell, apply_lower_momentum, apply_lower_force);

    for (int i = 0; i < g.ncells(); ++i) {
        g.set_momentumx(i, 1.);
        g.set_momentumy(i, 1.);
        g.set_momentumz(i, 1.);
        g.set_forcex(i, 2.);
        g.set_forcey(i, 2.);
        g.set_forcez(i, 2.);
    }

    g.apply_momentum_boundary_conditions(1, 0.);

    for (int i = 0; i < g.ngridx(); ++i)
        for (int j = 0; j < g.ngridy(); ++j)
            for (int k = 0; k < g.ngridz(); ++k) {
                if (k == 0) {
                    REQUIRE(g.momentumx(i, j, k)==-1.);
                    REQUIRE(g.momentumy(i, j, k)==-1.);
                    REQUIRE(g.momentumz(i, j, k)==-1.);
                } else {
                    REQUIRE(g.momentumx(i, j, k)==1.);
                    REQUIRE(g.momentumy(i, j, k)==1.);
                    REQUIRE(g.momentumz(i, j, k)==1.);
                }
            }

    g.apply_force_boundary_conditions(1, 0.);

    for (int i = 0; i < g.ngridx(); ++i)
        for (int j = 0; j < g.ngridy(); ++j)
            for (int k = 0; k < g.ngridz(); ++k) {
                if (k == 0) {
                    REQUIRE(g.forcex(i, j, k)==-2.);
                    REQUIRE(g.forcey(i, j, k)==-2.);
                    REQUIRE(g.forcez(i, j, k)==-2.);
                } else {
                    REQUIRE(g.forcex(i, j, k)==2.);
                    REQUIRE(g.forcey(i, j, k)==2.);
                    REQUIRE(g.forcez(i, j, k)==2.);
                }
            }

    // setting new momentum boundary condition function
    g.set_momentum_boundary_function(apply_west_momentum);
    g.apply_momentum_boundary_conditions(1, 0.);

    for (int i = 0; i < g.ngridx(); ++i)
        for (int j = 0; j < g.ngridy(); ++j)
            for (int k = 0; k < g.ngridz(); ++k) {
                if (i == 0) {
                    REQUIRE(g.momentumx(i, j, k)==-3.);
                    REQUIRE(g.momentumy(i, j, k)==-3.);
                    REQUIRE(g.momentumz(i, j, k)==-3.);
                } else if (k == 0) {
                    REQUIRE(g.momentumx(i, j, k)==-1.);
                    REQUIRE(g.momentumy(i, j, k)==-1.);
                    REQUIRE(g.momentumz(i, j, k)==-1.);
                } else {
                    REQUIRE(g.momentumx(i, j, k)==1.);
                    REQUIRE(g.momentumy(i, j, k)==1.);
                    REQUIRE(g.momentumz(i, j, k)==1.);
                }
            }

    // setting new force boundary condition function
    g.set_force_boundary_function(apply_west_force);
    g.apply_force_boundary_conditions(1, 0.);

    for (int i = 0; i < g.ngridx(); ++i)
        for (int j = 0; j < g.ngridy(); ++j)
            for (int k = 0; k < g.ngridz(); ++k) {
                if (i == 0) {
                    REQUIRE(g.forcex(i, j, k)==-4.);
                    REQUIRE(g.forcey(i, j, k)==-4.);
                    REQUIRE(g.forcez(i, j, k)==-4.);
                } else if (k == 0) {
                    REQUIRE(g.forcex(i, j, k)==-2.);
                    REQUIRE(g.forcey(i, j, k)==-2.);
                    REQUIRE(g.forcez(i, j, k)==-2.);
                } else {
                    REQUIRE(g.forcex(i, j, k)==2.);
                    REQUIRE(g.forcey(i, j, k)==2.);
                    REQUIRE(g.forcez(i, j, k)==2.);
                }
            }

}