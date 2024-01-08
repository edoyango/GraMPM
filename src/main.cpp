#include <grampm.hpp>
#include <grampm_kernels.hpp>
#include <array>
#include <grampm_integrators.hpp>
#include <cmath>
#include <grampm_stress_update_functions.hpp>

static constexpr double nbuff = 2.;

static void momentum_boundary(GraMPM::MPM_system<double> &self, const size_t &timestep, const double &dt) {

    // floor/ceiling
    #pragma omp for collapse(2)
    for (size_t i = 0; i < self.g_ngridx(); ++i) {
        for (size_t j = 0; j < self.g_ngridy(); ++j) {
            for (size_t k = 0; k < nbuff+1; ++k) {
                self.g_momentumx(i, j, k) = 0.;
                self.g_momentumy(i, j, k) = 0.;
                self.g_momentumz(i, j, k) = 0.;
            }
            for (size_t k = self.g_ngridz()-nbuff-1; k < self.g_ngridz(); ++k) {
                self.g_momentumx(i, j, k) = 0.;
                self.g_momentumy(i, j, k) = 0.;
                self.g_momentumz(i, j, k) = 0.;
            }
        }
    }

    // east/west wall
    #pragma omp for collapse(2)
    for (size_t j = 0; j < self.g_ngridy(); ++j) {
        for (size_t k = 0; k < self.g_ngridz(); ++k) {
            for (size_t i = 0; i < nbuff+1; ++i) {
                self.g_momentumx(i, j, k) = 0.;
                self.g_momentumy(i, j, k) = 0.;
                self.g_momentumz(i, j, k) = 0.;
            }
            for (size_t i = self.g_ngridx()-nbuff-1; i < self.g_ngridx(); ++i) {
                self.g_momentumx(i, j, k) = 0.;
                self.g_momentumy(i, j, k) = 0.;
                self.g_momentumz(i, j, k) = 0.;
            }
        }
    }

    // north/south wall
    #pragma omp for collapse(2)
    for (size_t i = 0; i < self.g_ngridx(); ++i) {
        for (size_t k = 0; k < self.g_ngridz(); ++k) {
            for (size_t j = 0; j < nbuff+1; ++j) {
                self.g_momentumy(i, j, k) = 0.;
            }
            for (size_t j = self.g_ngridy()-nbuff-1; j < self.g_ngridy(); ++j) {
                self.g_momentumy(i, j, k) = 0.;
            }
        }
    }
}

static void force_boundary(GraMPM::MPM_system<double> &self, const size_t &timestep, const double &dt) {
    // floor/ceiling
    #pragma omp for collapse(2)
    for (size_t i = 0; i < self.g_ngridx(); ++i) {
        for (size_t j = 0; j < self.g_ngridy(); ++j) {
            for (size_t k = 0; k < nbuff+1; ++k) {
                self.g_forcex(i, j, k) = 0.;
                self.g_forcey(i, j, k) = 0.;
                self.g_forcez(i, j, k) = 0.;
            }
            for (size_t k = self.g_ngridz()-nbuff-1; k < self.g_ngridz(); ++k) {
                self.g_forcex(i, j, k) = 0.;
                self.g_forcey(i, j, k) = 0.;
                self.g_forcez(i, j, k) = 0.;
            }
        }
    }

    // east/west wall
    #pragma omp for collapse(2)
    for (size_t j = 0; j < self.g_ngridy(); ++j) {
        for (size_t k = 0; k < self.g_ngridz(); ++k) {
            for (size_t i = 0; i < nbuff+1; ++i) {
                self.g_forcex(i, j, k) = 0.;
                self.g_forcey(i, j, k) = 0.;
                self.g_forcez(i, j, k) = 0.;
            }
            for (size_t i = self.g_ngridx()-nbuff-1; i < self.g_ngridx(); ++i) {
                self.g_forcex(i, j, k) = 0.;
                self.g_forcey(i, j, k) = 0.;
                self.g_forcez(i, j, k) = 0.;
            }
        }
    }

    // north/south wall
    #pragma omp for collapse(2)
    for (size_t i = 0; i < self.g_ngridx(); ++i) {
        for (size_t k = 0; k < self.g_ngridz(); ++k) {
            for (size_t j = 0; j < nbuff+1; ++j) {
                self.g_forcey(i, j, k) = 0.;
            }
            for (size_t j = self.g_ngridy()-nbuff-1; j < self.g_ngridy(); ++j) {
                self.g_forcey(i, j, k) = 0.;
            }
        }
    }
}

int main() {

    const double dcell=1.;
    const std::array<double, 3> mingrid {-nbuff*dcell, -nbuff*dcell, -nbuff*dcell}, 
        maxgrid {74.99+nbuff*dcell, 5.99+nbuff*dcell, 39.99+nbuff*dcell}, gf {0., 0., -9.81};
    // GraMPM::kernel_linear_bspline<double> knl(dcell);
    GraMPM::kernel_cubic_bspline<double> knl(dcell);
    GraMPM::MPM_system<double> myMPM(gf, knl, mingrid, maxgrid, dcell);
    // myMPM.set_stress_update_function(GraMPM::stress_update::drucker_prager_elastoplastic<double>);
    myMPM.set_stress_update_function(GraMPM::stress_update::eos<double>);
    myMPM.g_set_momentum_boundary_function(momentum_boundary);
    myMPM.g_set_force_boundary_function(force_boundary);

    const double rho_ini = 1000.;

    // generate particles
    for (int i = 0; i < 50; ++i) {
        for (int j = 0; j < 12; ++j) {
            for (int k = 0; k < 50; ++k) {
                GraMPM::particle<double> p;
                p.x[0] = (i+0.5)*dcell/2.;
                p.x[1] = (j+0.5)*dcell/2.;
                p.x[2] = (k+0.5)*dcell/2.;
                p.rho = rho_ini;
                p.mass = rho_ini*dcell*dcell*dcell/8.;
                myMPM.p_push_back(p);
            }
        }
    }

    const double E = 0.86e6, v = 0.3;
    const double K = E/(3.*(1.-2.*v)), G = E/(2.*(1.+v));
    // const double c = std::sqrt((K+4./3.*G)/rho_ini);
    const double c = 1500.;
    const double dt = dcell/c;

    const double pi = std::acos(-1.);
    myMPM.set_stress_update_param("E", E);
    myMPM.set_stress_update_param("v", v);
    myMPM.set_stress_update_param("phi", pi/9.);
    myMPM.set_stress_update_param("psi", 0.);
    myMPM.set_stress_update_param("cohesion", 0.);
    myMPM.set_stress_update_param("c", c);
    myMPM.set_stress_update_param("reference-density", rho_ini);
    myMPM.set_stress_update_param("kinematic viscosity", 0.01);

    myMPM.save_to_file("outputdata/p_", 0);

    GraMPM::integrators::MUSL<double>(myMPM, dt, 20000, 100, 100);

    return 0;
}