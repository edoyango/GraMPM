#include <grampm.hpp>
#include <grampm_kernels.hpp>
#include <array>
#include <grampm_integrators.hpp>
#include <cmath>
#include <grampm_stress_update_functions.hpp>

static void momentum_boundary(GraMPM::MPM_system<double> &self, const size_t &timestep, const double &dt) {

    // floor
    #pragma omp for collapse(2)
    for (size_t i = 0; i < self.g_ngridx(); ++i) {
        for (size_t j = 0; j < self.g_ngridy(); ++j) {
            self.g_momentumx(i, j, 0) = 0.;
            self.g_momentumx(i, j, 1) = 0.;
            self.g_momentumy(i, j, 0) = 0.;
            self.g_momentumy(i, j, 1) = 0.;
            self.g_momentumz(i, j, 0) = 0.;
            self.g_momentumz(i, j, 1) = 0.;
        }
    }

    // east/west wall
    const size_t xup = self.g_ngridx()-1;
    #pragma omp for collapse(2)
    for (size_t j = 0; j < self.g_ngridy(); ++j) {
        for (size_t k = 0; k < self.g_ngridz(); ++k) {
            self.g_momentumx(0, j, k) = 0.;
            self.g_momentumx(1, j, k) = 0.;
            self.g_momentumx(xup-1, j, k) = 0.;
            self.g_momentumx(xup, j, k) = 0.;
        }
    }

    // north/south wall
    const size_t yup = self.g_ngridy()-1;
    #pragma omp for collapse(2)
    for (size_t i = 0; i < self.g_ngridx(); ++i) {
        for (size_t k = 0; k < self.g_ngridz(); ++k) {
            self.g_momentumy(i, 0, k) = 0.;
            self.g_momentumy(i, 1, k) = 0.;
            self.g_momentumy(i, yup-1, k) = 0.;
            self.g_momentumy(i, yup, k) = 0.;
        }
    }
}

static void force_boundary(GraMPM::MPM_system<double> &self, const size_t &timestep, const double &dt) {
    // floor
    #pragma omp for collapse(2)
    for (size_t i = 0; i < self.g_ngridx(); ++i) {
        for (size_t j = 0; j < self.g_ngridy(); ++j) {
            self.g_forcex(i, j, 0) = 0.;
            self.g_forcex(i, j, 1) = 0.;
            self.g_forcey(i, j, 0) = 0.;
            self.g_forcey(i, j, 1) = 0.;
            self.g_forcez(i, j, 0) = 0.;
            self.g_forcez(i, j, 1) = 0.;
        }
    }

    // east/west wall
    const size_t xup = self.g_ngridx()-1;
    #pragma omp for collapse(2)
    for (size_t j = 0; j < self.g_ngridy(); ++j) {
        for (size_t k = 0; k < self.g_ngridz(); ++k) {
            self.g_forcex(0, j, k) = 0.;
            self.g_forcex(1, j, k) = 0.;
            self.g_forcex(xup-1, j, k) = 0.;
            self.g_forcex(xup, j, k) = 0.;
        }
    }

    // north/south wall
    const size_t yup = self.g_ngridy()-1;
    #pragma omp for collapse(2)
    for (size_t i = 0; i < self.g_ngridx(); ++i) {
        for (size_t k = 0; k < self.g_ngridz(); ++k) {
            self.g_forcey(i, 0, k) = 0.;
            self.g_forcey(i, 1, k) = 0.;
            self.g_forcey(i, yup-1, k) = 0.;
            self.g_forcey(i, yup, k) = 0.;
        }
    }
}

int main() {

    const double dcell=0.002;
    const std::array<double, 3> mingrid {-dcell, -dcell, -dcell}, 
        maxgrid {0.299+dcell, 0.011+dcell, 0.049+dcell}, gf {0., 0., -9.81};
    // GraMPM::kernel_linear_bspline<double> knl(dcell);
    GraMPM::kernel_cubic_bspline<double> knl(dcell);
    GraMPM::MPM_system<double> myMPM(gf, knl, mingrid, maxgrid, dcell);
    myMPM.set_stress_update_function(GraMPM::stress_update::drucker_prager_elastoplastic<double>);
    myMPM.g_set_momentum_boundary_function(momentum_boundary);
    myMPM.g_set_force_boundary_function(force_boundary);

    const double rho_ini = 1650.;

    // generate particles
    for (int i = 0; i < 100; ++i) {
        for (int j = 0; j < 12; ++j) {
            for (int k = 0; k < 50; ++k) {
                GraMPM::particle<double> p;
                p.x = (i+0.5)*dcell/2.;
                p.y = (j+0.5)*dcell/2.;
                p.z = (k+0.5)*dcell/2.;
                p.rho = rho_ini;
                p.mass = rho_ini*dcell*dcell*dcell/8.;
                myMPM.p_push_back(p);
            }
        }
    }

    const double E = 0.86e6, v = 0.3;
    const double K = E/(3.*(1.-2.*v)), G = E/(2.*(1.+v));
    const double c = std::sqrt((K+4./3.*G)/rho_ini);
    const double dt = dcell/c;

    myMPM.p_E() = E;
    myMPM.p_v() = v;
    const double pi = std::acos(-1.);
    myMPM.set_DP_params(pi/9., 0., 0.);

    myMPM.save_to_file("outputdata/p_", 0);

    GraMPM::integrators::MUSL<double>(myMPM, dt, 500, 500, 600);

    return 0;
}