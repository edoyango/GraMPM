#ifndef GRAMPM_KERNELS
#define GRAMPM_KERNELS

#include <cmath>

// fast (?) sign function from SO https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T>
static int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

namespace GraMPM {

    template<typename F>
    class kernel_base {
        public:
            kernel_base(const F dc, const F r)
                : dcell(dc)
                , radius(r)
            {
            }
            const double radius, dcell;
            const F w(const F &dx, const F &dy, const F &dz) const {
                return w1(dx)*w1(dy)*w1(dz);
            }
            virtual void dwdx(const F &dx, const F &dy, const F &dz, F &dwdx, F &dwdy, F &dwdz) const 
            { 
                dwdx = 0.;
                dwdy = 0.;
                dwdz = 0.;
            }
        protected:
            virtual double w1(const double &dr) const {return 0.;}
    };

    template<typename F>
    class kernel_linear_bspline : public kernel_base<F> {
        protected:
            F w1(const F &dr) const override {
                return std::max(0., 1.-std::abs(dr)/kernel_base<F>::dcell);
            }
        public:
            kernel_linear_bspline(const F dc)
                : kernel_base<F>(dc, 1.)
            {
            }

            void dwdx(const F &dx, const F &dy, const F &dz, F &dwdx, F &dwdy, F &dwdz) const override {
                dwdx = -sgn(dx)/kernel_base<F>::dcell*w1(dy)*w1(dz);
                dwdy = -sgn(dy)/kernel_base<F>::dcell*w1(dx)*w1(dz);
                dwdz = -sgn(dz)/kernel_base<F>::dcell*w1(dx)*w1(dy);
            }
    };

}
#endif