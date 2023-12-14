#include <grampm.hpp>
#include <grampm_kernels.hpp>
#include <array>
#include <grampm_integrators.hpp>
#include <cmath>
#include <grampm_stress_update_functions.hpp>

static void momentum_boundary(GraMPM::grid<double> &self, const int &timestep, const double &dt) {

    // floor
    for (int i = 0; i < self.ngridx(); ++i)
        for (int j = 0; j < self.ngridy(); ++j) {
            self.set_momentumx(i, j, 0, 0.);
            self.set_momentumx(i, j, 1, 0.);
            self.set_momentumy(i, j, 0, 0.);
            self.set_momentumy(i, j, 1, 0.);
            self.set_momentumz(i, j, 0, 0.);
            self.set_momentumz(i, j, 1, 0.);
        }

    // east/west wall
    const int xup = self.ngridx()-1;
    for (int j = 0; j < self.ngridy(); ++j)
        for (int k = 0; k < self.ngridz(); ++k) {
            self.set_momentumx(0, j, k, 0.);
            self.set_momentumx(1, j, k, 0.);
            self.set_momentumx(xup-1, j, k, 0.);
            self.set_momentumx(xup, j, k, 0.);
        }

    // north/south wall
    const int yup = self.ngridy()-1;
    for (int i = 0; i < self.ngridx(); ++i) 
        for (int k = 0; k < self.ngridz(); ++k) {
            self.set_momentumy(i, 0, k, 0.);
            self.set_momentumy(i, 1, k, 0.);
            self.set_momentumy(i, yup-1, k, 0.);
            self.set_momentumy(i, yup, k, 0.);
        }
}

static void force_boundary(GraMPM::grid<double> &self, const int &timestep, const double &dt) {
    // floor
    for (int i = 0; i < self.ngridx(); ++i)
        for (int j = 0; j < self.ngridy(); ++j) {
            self.set_forcex(i, j, 0, 0.);
            self.set_forcex(i, j, 1, 0.);
            self.set_forcey(i, j, 0, 0.);
            self.set_forcey(i, j, 1, 0.);
            self.set_forcez(i, j, 0, 0.);
            self.set_forcez(i, j, 1, 0.);
        }

    // east/west wall
    const int xup = self.ngridx()-1;
    for (int j = 0; j < self.ngridy(); ++j)
        for (int k = 0; k < self.ngridz(); ++k) {
            self.set_forcex(0, j, k, 0.);
            self.set_forcex(1, j, k, 0.);
            self.set_forcex(xup-1, j, k, 0.);
            self.set_forcex(xup, j, k, 0.);
        }

    // north/south wall
    const int yup = self.ngridy()-1;
    for (int i = 0; i < self.ngridx(); ++i) 
        for (int k = 0; k < self.ngridz(); ++k) {
            self.set_forcey(i, 0, k, 0.);
            self.set_forcey(i, 1, k, 0.);
            self.set_forcey(i, yup-1, k, 0.);
            self.set_forcey(i, yup, k, 0.);
        }
}

int main() {

    const int maxn = 10;
    const double dcell=0.002;
    const std::array<double, 3> mingrid {-dcell, -dcell, -dcell}, 
        maxgrid {0.299+dcell, 0.011+dcell, 0.049+dcell}, gf {0., 0., -9.81};
    // GraMPM::kernel_linear_bspline<double> knl(dcell);
    GraMPM::kernel_cubic_bspline<double> knl(dcell);
    GraMPM::grid<double> g(mingrid, maxgrid, dcell, momentum_boundary, force_boundary);
    GraMPM::particle_system<double> ps(g, knl);
    ps.set_stress_update_function(GraMPM::stress_update::drucker_prager_elastoplastic<double>);

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
                ps.push_back(p);
            }
        }
    }

    const double E = 0.86e6, v = 0.3;
    const double K = E/(3.*(1.-2.*v)), G = E/(2.*(1.+v));
    const double c = std::sqrt((K+4./3.*G)/rho_ini);
    const double dt = dcell/c;

    ps.set_body_force(gf);
    ps.set_E(E);
    ps.set_v(v);
    const double pi = std::acos(-1.);
    ps.set_DP_params(pi/9., 0., 0.);

    ps.save_to_file("outputdata/p_", 0);

    GraMPM::integrators::MUSL<double>(ps, dt, 6000, 50, 50);

    return 0;
}