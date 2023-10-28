#include "grampm.hpp"
#include <grampm_kernels.hpp>
#include <array>

int main() {

    const int maxn = 10;
    const std::array<double, 3> mingrid {0., 0., 0.}, maxgrid {0.3, 0.1, 0.05};
    const double dcell=0.004;
    GraMPM::kernel_linear_bspline<double> knl(dcell);
    GraMPM::grid<double> g(mingrid, maxgrid, dcell);
    GraMPM::particle_system<double> p(g, knl);
    return 0;
}