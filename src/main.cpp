#include <grampm.hpp>
#include <grampm_kernels.hpp>
#include <array>
#include <grampm_integrators.hpp>
#include <cmath>
#include <grampm_stress_update_functions.hpp>

static void momentum_boundary(GraMPM::MPM_system<double> &self, const size_t &timestep, const double &dt) {

    // floor/ceiling
    const size_t zup = self.g_ngridz()-2;
    #pragma omp for collapse(2)
    for (size_t i = 0; i < self.g_ngridx(); ++i) {
        for (size_t j = 0; j < self.g_ngridy(); ++j) {
            self.g_momentumx(i, j, 0) = 0.;
            self.g_momentumx(i, j, 1) = 0.;
            self.g_momentumx(i, j, 2) = 0.;
            self.g_momentumy(i, j, 0) = 0.;
            self.g_momentumy(i, j, 1) = 0.;
            self.g_momentumy(i, j, 2) = 0.;
            self.g_momentumz(i, j, 0) = 0.;
            self.g_momentumz(i, j, 1) = 0.;
            self.g_momentumz(i, j, 2) = 0.;
            self.g_momentumx(i, j, zup-2) = 0.;
            self.g_momentumx(i, j, zup-1) = 0.;
            self.g_momentumx(i, j, zup) = 0.;
            self.g_momentumy(i, j, zup-2) = 0.;
            self.g_momentumy(i, j, zup-1) = 0.;
            self.g_momentumy(i, j, zup) = 0.;
            self.g_momentumz(i, j, zup-2) = 0.;
            self.g_momentumz(i, j, zup-1) = 0.;
            self.g_momentumz(i, j, zup) = 0.;
        }
    }

    // east/west wall
    const size_t xup = self.g_ngridx()-1;
    #pragma omp for collapse(2)
    for (size_t j = 0; j < self.g_ngridy(); ++j) {
        for (size_t k = 0; k < self.g_ngridz(); ++k) {
            self.g_momentumx(0, j, k) = 0.;
            self.g_momentumx(1, j, k) = 0.;
            self.g_momentumx(2, j, k) = 0.;
            self.g_momentumy(0, j, k) = 0.;
            self.g_momentumy(1, j, k) = 0.;
            self.g_momentumy(2, j, k) = 0.;
            self.g_momentumz(0, j, k) = 0.;
            self.g_momentumz(1, j, k) = 0.;
            self.g_momentumz(2, j, k) = 0.;
            self.g_momentumx(xup-2, j, k) = 0.;
            self.g_momentumx(xup-1, j, k) = 0.;
            self.g_momentumx(xup, j, k) = 0.;
            self.g_momentumy(xup-2, j, k) = 0.;
            self.g_momentumy(xup-1, j, k) = 0.;
            self.g_momentumy(xup, j, k) = 0.;
            self.g_momentumz(xup-2, j, k) = 0.;
            self.g_momentumz(xup-1, j, k) = 0.;
            self.g_momentumz(xup, j, k) = 0.;
        }
    }

    // north/south wall
    const size_t yup = self.g_ngridy()-1;
    #pragma omp for collapse(2)
    for (size_t i = 0; i < self.g_ngridx(); ++i) {
        for (size_t k = 0; k < self.g_ngridz(); ++k) {
            // self.g_momentumx(i, 0, k) = 0.;
            // self.g_momentumx(i, 1, k) = 0.;
            // self.g_momentumx(i, 2, k) = 0.;
            self.g_momentumy(i, 0, k) = 0.;
            self.g_momentumy(i, 1, k) = 0.;
            self.g_momentumy(i, 2, k) = 0.;
            // self.g_momentumz(i, 0, k) = 0.;
            // self.g_momentumz(i, 1, k) = 0.;
            // self.g_momentumz(i, 2, k) = 0.;
            // self.g_momentumx(i, yup-2, k) = 0.;
            // self.g_momentumx(i, yup-1, k) = 0.;
            // self.g_momentumx(i, yup, k) = 0.;
            self.g_momentumy(i, yup-2, k) = 0.;
            self.g_momentumy(i, yup-1, k) = 0.;
            self.g_momentumy(i, yup, k) = 0.;
            // self.g_momentumz(i, yup-2, k) = 0.;
            // self.g_momentumz(i, yup-1, k) = 0.;
            // self.g_momentumz(i, yup, k) = 0.;
        }
    }
}

static void force_boundary(GraMPM::MPM_system<double> &self, const size_t &timestep, const double &dt) {
    // floor/ceiling
    const size_t zup = self.g_ngridz()-1;
    #pragma omp for collapse(2)
    for (size_t i = 0; i < self.g_ngridx(); ++i) {
        for (size_t j = 0; j < self.g_ngridy(); ++j) {
            self.g_forcex(i, j, 0) = 0.;
            self.g_forcex(i, j, 1) = 0.;
            self.g_forcex(i, j, 2) = 0.;
            self.g_forcey(i, j, 0) = 0.;
            self.g_forcey(i, j, 1) = 0.;
            self.g_forcey(i, j, 2) = 0.;
            self.g_forcez(i, j, 0) = 0.;
            self.g_forcez(i, j, 1) = 0.;
            self.g_forcez(i, j, 2) = 0.;
            self.g_forcex(i, j, zup-2) = 0.;
            self.g_forcex(i, j, zup-1) = 0.;
            self.g_forcex(i, j, zup) = 0.;
            self.g_forcey(i, j, zup-2) = 0.;
            self.g_forcey(i, j, zup-1) = 0.;
            self.g_forcey(i, j, zup) = 0.;
            self.g_forcez(i, j, zup-2) = 0.;
            self.g_forcez(i, j, zup-1) = 0.;
            self.g_forcez(i, j, zup) = 0.;
        }
    }

    // east/west wall
    const size_t xup = self.g_ngridx()-1;
    #pragma omp for collapse(2)
    for (size_t j = 0; j < self.g_ngridy(); ++j) {
        for (size_t k = 0; k < self.g_ngridz(); ++k) {
            self.g_forcex(0, j, k) = 0.;
            self.g_forcex(1, j, k) = 0.;
            self.g_forcex(2, j, k) = 0.;
            self.g_forcey(0, j, k) = 0.;
            self.g_forcey(1, j, k) = 0.;
            self.g_forcey(2, j, k) = 0.;
            self.g_forcez(0, j, k) = 0.;
            self.g_forcez(1, j, k) = 0.;
            self.g_forcez(2, j, k) = 0.;
            self.g_forcex(xup-2, j, k) = 0.;
            self.g_forcex(xup-1, j, k) = 0.;
            self.g_forcex(xup, j, k) = 0.;
            self.g_forcey(xup-2, j, k) = 0.;
            self.g_forcey(xup-1, j, k) = 0.;
            self.g_forcey(xup, j, k) = 0.;
            self.g_forcez(xup-2, j, k) = 0.;
            self.g_forcez(xup-1, j, k) = 0.;
            self.g_forcez(xup, j, k) = 0.;
        }
    }

    // north/south wall
    const size_t yup = self.g_ngridy()-1;
    #pragma omp for collapse(2)
    for (size_t i = 0; i < self.g_ngridx(); ++i) {
        for (size_t k = 0; k < self.g_ngridz(); ++k) {
            // self.g_forcex(i, 0, k) = 0.;
            // self.g_forcex(i, 1, k) = 0.;
            // self.g_forcex(i, 2, k) = 0.;
            self.g_forcey(i, 0, k) = 0.;
            self.g_forcey(i, 1, k) = 0.;
            self.g_forcey(i, 2, k) = 0.;
            // self.g_forcez(i, 0, k) = 0.;
            // self.g_forcez(i, 1, k) = 0.;
            // self.g_forcez(i, 2, k) = 0.;
            // self.g_forcex(i, yup-2, k) = 0.;
            // self.g_forcex(i, yup-1, k) = 0.;
            // self.g_forcex(i, yup, k) = 0.;
            self.g_forcey(i, yup-2, k) = 0.;
            self.g_forcey(i, yup-1, k) = 0.;
            self.g_forcey(i, yup, k) = 0.;
            // self.g_forcez(i, yup-2, k) = 0.;
            // self.g_forcez(i, yup-1, k) = 0.;
            // self.g_forcez(i, yup, k) = 0.;
        }
    }
}

int main() {

    const double dcell=1.;
    const std::array<double, 3> mingrid {-2.*dcell, -2.*dcell, -2.*dcell}, 
        maxgrid {74.99+2.*dcell, 5.99+2.*dcell, 39.99+2.*dcell}, gf {0., 0., -9.81};
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

    myMPM.save_to_file("outputdata/p_", 0);

    GraMPM::integrators::MUSL<double>(myMPM, dt, 20000, 100, 100);

    return 0;
}