#ifndef GRAMPM_KERNELS
#define GRAMPM_KERNELS

#include <cmath>

// fast (?) sign function from SO https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T>
static int sgn(T val) {
    return (T(0.) < val) - (val < T(0.));
}

namespace GraMPM {

    template<typename F>
    class kernel_base {
        public:
            kernel_base(const F dc, const F r)
                : radius(r)
                , dcell(dc)
            {
            }
            const F radius, dcell;
            void w_dwdx(const F &dx, const F&dy, const F &dz, F &w, F &dwdx, F &dwdy, F &dwdz) const {
                const F qx = std::abs(dx)/dcell, qy = std::abs(dy)/dcell, qz = std::abs(dz)/dcell;
                const F w1x = w1(qx), w1y = w1(qy), w1z = w1(qz);
                w = w1x*w1y*w1z;
                dwdx = w1y*w1z*dw1dq(qx)*dqdr(dx);
                dwdy = w1x*w1z*dw1dq(qy)*dqdr(dy);
                dwdz = w1x*w1y*dw1dq(qz)*dqdr(dz);
            }
        protected:
            virtual F w1(const F &q) const {return 0.;}
            virtual F dw1dq(const F &q) const {return 0.;}
            F dqdr(const F &dr) const {return sgn(dr)/dcell;}
    };

    template<typename F>
    class kernel_linear_bspline : public kernel_base<F> {
        protected:
            F w1(const F &q) const override {
                return std::max(0., 1.-q);
            }

            F dw1dq (const F &q) const override {
                return -1.;
            }
        public:
            kernel_linear_bspline(const F dc)
                : kernel_base<F>(dc, 1.)
            {
            }
    };

    template<typename F>
    class kernel_cubic_bspline : public kernel_base<F> {
        protected:
            F w1(const F &q) const override {
                const F dim2q {std::max(0., 2.-q)}, dim1q {std::max(0., 1.-q)};
                return 2./3.*(0.25*dim2q*dim2q*dim2q-dim1q*dim1q*dim1q);
            }

            F dw1dq (const F &q) const override {
                const F dim2q {std::max(0., 2.-q)}, dim1q {std::max(0., 1.-q)};
                return -2.*(0.25*dim2q*dim2q-dim1q*dim1q);
            }
        public:
            kernel_cubic_bspline(const F dc)
                : kernel_base<F>(dc, 2.)
            {
            }
    };

}
#endif