#include "grampm.hpp"
#include <grampm_kernels.hpp>

int main() {

    const int maxn = 10;
    const double mingrid[3] {0., 0., 0.}, maxgrid[3] {0.3, 0.1, 0.05}, dcell=0.004;
    const kernel_linear_bspline<double> knl(dcell);
    GraMPM<double> myMPM(10, maxgrid, mingrid, dcell, knl);
    return 0;
}