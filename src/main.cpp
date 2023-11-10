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
            self.set_momentumy(i, j, 0, 0.);
            self.set_momentumz(i, j, 0, 0.);
        }

    // east/west wall
    const int xup = self.ngridx()-1;
    for (int j = 0; j < self.ngridy(); ++j)
        for (int k = 0; k < self.ngridz(); ++k) {
            self.set_momentumx(0, j, k, 0.);
            self.set_momentumx(xup, j, k, 0.);
        }

    // north/south wall
    const int yup = self.ngridy()-1;
    for (int i = 0; i < self.ngridx(); ++i) 
        for (int k = 0; k < self.ngridz(); ++k) {
            self.set_momentumy(i, 0, k, 0.);
            self.set_momentumy(i, yup, k, 0.);
        }
}

static void force_boundary(GraMPM::grid<double> &self, const int &timestep, const double &dt) {
    // floor
    for (int i = 0; i < self.ngridx(); ++i)
        for (int j = 0; j < self.ngridy(); ++j) {
            self.set_forcex(i, j, 0, 0.);
            self.set_forcey(i, j, 0, 0.);
            self.set_forcez(i, j, 0, 0.);
        }

    // east/west wall
    const int xup = self.ngridx()-1;
    for (int j = 0; j < self.ngridy(); ++j)
        for (int k = 0; k < self.ngridz(); ++k) {
            self.set_forcex(0, j, k, 0.);
            self.set_forcex(xup, j, k, 0.);
        }

    // north/south wall
    const int yup = self.ngridy()-1;
    for (int i = 0; i < self.ngridx(); ++i) 
        for (int k = 0; k < self.ngridz(); ++k) {
            self.set_forcey(i, 0, k, 0.);
            self.set_forcey(i, yup, k, 0.);
        }
}

int main() {

    const int maxn = 10;
    const std::array<double, 3> mingrid {0., 0., 0.}, maxgrid {0.299, 0.019, 0.049}, gf {0., 0., -9.81};
    const double dcell=0.004;
    GraMPM::kernel_linear_bspline<double> knl(dcell);
    GraMPM::grid<double> g(mingrid, maxgrid, dcell, momentum_boundary, force_boundary);
    GraMPM::particle_system<double> ps(g, knl);
    ps.set_stress_update_function(GraMPM::stress_update::hookes_law<double>);

    const double rho_ini = 1650.;

    // generate particles
    for (int i = 0; i < 50; ++i) {
        for (int j = 0; j < 10; ++j) {
            for (int k = 0; k < 25; ++k) {
                GraMPM::particle<double> p;
                p.x = mingrid[0] + (i+0.5)*dcell/2.;
                p.y = mingrid[1] + (j+0.5)*dcell/2.;
                p.z = mingrid[2] + (k+0.5)*dcell/2.;
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
    ps.set_E(0.86e6);
    ps.set_v(0.3);

    ps.save_to_file("outputdata/p_", 0);

    GraMPM::integrators::MUSL<double>(ps, dt, 100, 1, 1);

    return 0;
}