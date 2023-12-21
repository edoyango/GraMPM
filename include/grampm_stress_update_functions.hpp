#ifndef GRAMPM_STRESS_UPDATE
#define GRAMPM_STRESS_UPDATE

#include <grampm.hpp>

namespace GraMPM {
    namespace stress_update {

        template<typename F>
        void hookes_law(GraMPM::particle_system<F> &self, const F &dt) {

            const F D0 {self.m_E/((1.+self.m_v)*(1.-2.*self.m_v))};

            // DE*dstrain
            #pragma omp for
            for (int i = 0; i < self.m_size; ++i) {
                // jaumann stress rate (applied first, as it depends on current stress state)
                self.m_sigmaxx[i] -= dt*2.*(self.m_spinratexy[i]*self.m_sigmaxy[i] + self.m_spinratexz[i]*self.m_sigmaxz[i]);
                self.m_sigmayy[i] -= dt*2.*(-self.m_spinratexy[i]*self.m_sigmaxy[i] + self.m_spinrateyz[i]*self.m_sigmayz[i]);
                self.m_sigmazz[i] += dt*2.*(self.m_spinratexz[i]*self.m_sigmaxz[i] + self.m_spinrateyz[i]*self.m_sigmayz[i]);
                self.m_sigmaxy[i] += dt*(self.m_sigmaxx[i]*self.m_spinratexy[i] - self.m_sigmaxz[i]*self.m_spinrateyz[i] -
                    self.m_spinratexy[i]*self.m_sigmayy[i] - self.m_spinratexz[i]*self.m_sigmayz[i]);
                self.m_sigmaxz[i] += dt*(self.m_sigmaxx[i]*self.m_spinratexz[i] + self.m_sigmaxy[i]*self.m_spinrateyz[i] -
                    self.m_spinratexy[i]*self.m_sigmayz[i] - self.m_spinratexz[i]*self.m_sigmazz[i]);
                self.m_sigmayz[i] += dt*(self.m_sigmaxy[i]*self.m_spinratexz[i] + self.m_sigmayy[i]*self.m_spinrateyz[i] +
                    self.m_spinratexy[i]*self.m_sigmaxz[i] - self.m_spinrateyz[i]*self.m_sigmazz[i]);

                // update stress state with elastic increment
                self.m_sigmaxx[i] += dt*D0*((1.-self.m_v)*self.m_strainratexx[i] + self.m_v*self.m_strainrateyy[i] + self.m_v*self.m_strainratezz[i]);
                self.m_sigmayy[i] += dt*D0*(self.m_v*self.m_strainratexx[i] + (1.-self.m_v)*self.m_strainrateyy[i] + self.m_v*self.m_strainratezz[i]);
                self.m_sigmazz[i] += dt*D0*(self.m_v*self.m_strainratexx[i] + self.m_v*self.m_strainrateyy[i] + (1.-self.m_v)*self.m_strainratezz[i]);
                self.m_sigmaxy[i] += dt*D0*self.m_strainratexy[i]*(1.-2.*self.m_v);
                self.m_sigmaxz[i] += dt*D0*self.m_strainratexz[i]*(1.-2.*self.m_v);
                self.m_sigmayz[i] += dt*D0*self.m_strainrateyz[i]*(1.-2.*self.m_v);
            }

        }

        template<typename F>
        void drucker_prager_elastoplastic(GraMPM::particle_system<F> &self, const F &dt) {

            // performing initial elastic predictor step
            hookes_law(self, dt);

            // setting up for plastic corrector step (duplicated code, would like to remove this somehow)
            const F D0 {self.m_E/((1.+self.m_v)*(1.-2.*self.m_v))};

            F phi, psi, coh, alpha_phi, alpha_psi, k_c;
            self.DP_params(phi, psi, coh, alpha_phi, alpha_psi, k_c);

            // begin plastic correction
            #pragma omp for
            for (int i = 0; i < self.m_size; ++i) {
                // calculating invariants and deviatoric stress tensor
                F I1 = self.m_sigmaxx[i] + self.m_sigmayy[i] + self.m_sigmazz[i];
                const F s[6] {
                    self.m_sigmaxx[i] - I1/3.,
                    self.m_sigmayy[i] - I1/3.,
                    self.m_sigmazz[i] - I1/3.,
                    self.m_sigmaxy[i],
                    self.m_sigmaxz[i],
                    self.m_sigmayz[i],
                };
                const F J2 = 0.5*(s[0]*s[0] + s[1]*s[1] + s[2]*s[2] + 2.*(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]));
                // tensile correction 1
                if (J2 == 0. && I1 > k_c/alpha_phi) {
                    self.m_sigmaxx[i] = k_c/alpha_phi/3.;
                    self.m_sigmayy[i] = k_c/alpha_phi/3.;
                    self.m_sigmazz[i] = k_c/alpha_phi/3.;
                    I1 = k_c/alpha_phi;
                }

                // calculate yield function value
                // max used to reduce branching
                const F f {alpha_phi*I1 + std::sqrt(J2) - k_c};
                
                if (f > 1.e-13) {
                    const F snorm = 2.*std::sqrt(J2);
                    const F shat[6] {
                        s[0]/snorm,
                        s[1]/snorm,
                        s[2]/snorm,
                        s[3]/snorm,
                        s[4]/snorm,
                        s[5]/snorm,
                    };
                    const F dfdsig[6] {
                        alpha_phi + shat[0],
                        alpha_phi + shat[1],
                        alpha_phi + shat[2],
                        shat[3],
                        shat[4],
                        shat[5],
                    };
                    const F dgdsig[6] {
                        alpha_psi + shat[0],
                        alpha_psi + shat[1],
                        alpha_psi + shat[2],
                        shat[3],
                        shat[4],
                        shat[5],
                    };
                    const F dlambda {f/(
                        dfdsig[0]*D0*((1.-self.m_v)*dgdsig[0] + self.m_v*dgdsig[1] + self.m_v*dgdsig[2]) +
                        dfdsig[1]*D0*(self.m_v*dgdsig[0] + (1.-self.m_v)*dgdsig[1] + self.m_v*dgdsig[2]) +
                        dfdsig[2]*D0*(self.m_v*dgdsig[0] + self.m_v*dgdsig[1] + (1.-self.m_v)*dgdsig[2]) +
                        2.*dfdsig[3]*D0*dgdsig[3]*(1.-2.*self.m_v) +
                        2.*dfdsig[4]*D0*dgdsig[4]*(1.-2.*self.m_v) +
                        2.*dfdsig[5]*D0*dgdsig[5]*(1.-2.*self.m_v)
                    )};
                    
                    self.m_sigmaxx[i] -= dlambda*D0*((1.-self.m_v)*dgdsig[0] + self.m_v*dgdsig[1] + self.m_v*dgdsig[2]);
                    self.m_sigmayy[i] -= dlambda*D0*(self.m_v*dgdsig[0] + (1.-self.m_v)*dgdsig[1] + self.m_v*dgdsig[2]);
                    self.m_sigmazz[i] -= dlambda*D0*(self.m_v*dgdsig[0] + self.m_v*dgdsig[1] + (1.-self.m_v)*dgdsig[2]);
                    self.m_sigmaxy[i] -= dlambda*D0*dgdsig[3]*(1.-2.*self.m_v);
                    self.m_sigmaxz[i] -= dlambda*D0*dgdsig[4]*(1.-2.*self.m_v);
                    self.m_sigmayz[i] -= dlambda*D0*dgdsig[5]*(1.-2.*self.m_v);

                    I1 = self.m_sigmaxx[i] + self.m_sigmayy[i] + self.m_sigmazz[i];
                    if (I1 > k_c/alpha_phi) {
                        self.m_sigmaxx[i] = k_c/alpha_phi/3.;
                        self.m_sigmayy[i] = k_c/alpha_phi/3.;
                        self.m_sigmazz[i] = k_c/alpha_phi/3.;
                        self.m_sigmaxy[i] = 0.;
                        self.m_sigmaxz[i] = 0.;
                        self.m_sigmayz[i] = 0.;
                    }
                }

            }
        }
    }
}
#endif