#ifndef GRAMPM_STRESS_UPDATE
#define GRAMPM_STRESS_UPDATE

#include <grampm.hpp>

namespace GraMPM {
    namespace stress_update {

        template<typename F>
        void hookes_law(GraMPM::MPM_system<F> &self, const F &dt) {

            const F D0 {self.m_E/((1.+self.m_v)*(1.-2.*self.m_v))};

            std::vector<F> dsigmaxx(self.p_size), dsigmayy(self.p_size), dsigmazz(self.p_size), dsigmaxy(self.p_size), 
                dsigmaxz(self.p_size), dsigmayz(self.p_size);

            // DE*dstrain
            #pragma omp for
            for (int i = 0; i < self.p_size; ++i) {
                // jaumann stress rate (applied first, as it depends on current stress state)
                dsigmaxx[i] -= 2.*(self.p_spinratexy[i]*self.p_sigmaxy[i] + self.p_spinratexz[i]*self.p_sigmaxz[i]);
                dsigmayy[i] -= 2.*(-self.p_spinratexy[i]*self.p_sigmaxy[i] + self.p_spinrateyz[i]*self.p_sigmayz[i]);
                dsigmazz[i] += 2.*(self.p_spinratexz[i]*self.p_sigmaxz[i] + self.p_spinrateyz[i]*self.p_sigmayz[i]);
                dsigmaxy[i] += self.p_sigmaxx[i]*self.p_spinratexy[i] - self.p_sigmaxz[i]*self.p_spinrateyz[i] -
                    self.p_spinratexy[i]*self.p_sigmayy[i] - self.p_spinratexz[i]*self.p_sigmayz[i];
                dsigmaxz[i] += self.p_sigmaxx[i]*self.p_spinratexz[i] + self.p_sigmaxy[i]*self.p_spinrateyz[i] -
                    self.p_spinratexy[i]*self.p_sigmayz[i] - self.p_spinratexz[i]*self.p_sigmazz[i];
                dsigmayz[i] += self.p_sigmaxy[i]*self.p_spinratexz[i] + self.p_sigmayy[i]*self.p_spinrateyz[i] +
                    self.p_spinratexy[i]*self.p_sigmaxz[i] - self.p_spinrateyz[i]*self.p_sigmazz[i];

                // update stress state with elastic increment
                dsigmaxx[i] += D0*((1.-self.m_v)*self.p_strainratexx[i] + self.m_v*self.p_strainrateyy[i] + self.m_v*self.p_strainratezz[i]);
                dsigmayy[i] += D0*(self.m_v*self.p_strainratexx[i] + (1.-self.m_v)*self.p_strainrateyy[i] + self.m_v*self.p_strainratezz[i]);
                dsigmazz[i] += D0*(self.m_v*self.p_strainratexx[i] + self.m_v*self.p_strainrateyy[i] + (1.-self.m_v)*self.p_strainratezz[i]);
                dsigmaxy[i] += D0*self.p_strainratexy[i]*(1.-2.*self.m_v);
                dsigmaxz[i] += D0*self.p_strainratexz[i]*(1.-2.*self.m_v);
                dsigmayz[i] += D0*self.p_strainrateyz[i]*(1.-2.*self.m_v);

                self.p_sigmaxx[i] += dt*dsigmaxx[i];
                self.p_sigmayy[i] += dt*dsigmayy[i];
                self.p_sigmazz[i] += dt*dsigmazz[i];
                self.p_sigmaxy[i] += dt*dsigmaxy[i];
                self.p_sigmaxz[i] += dt*dsigmaxz[i];
                self.p_sigmayz[i] += dt*dsigmayz[i];
            }

        }

        template<typename F>
        void drucker_prager_elastoplastic(GraMPM::MPM_system<F> &self, const F &dt) {

            // performing initial elastic predictor step
            hookes_law(self, dt);

            // setting up for plastic corrector step (duplicated code, would like to remove this somehow)
            const F D0 {self.m_E/((1.+self.m_v)*(1.-2.*self.m_v))};

            F phi, psi, coh, alpha_phi, alpha_psi, k_c;
            self.DP_params(phi, psi, coh, alpha_phi, alpha_psi, k_c);

            // begin plastic correction
            #pragma omp for
            for (int i = 0; i < self.p_size; ++i) {
                // calculating invariants and deviatoric stress tensor
                F I1 = self.p_sigmaxx[i] + self.p_sigmayy[i] + self.p_sigmazz[i];
                const F s[6] {
                    self.p_sigmaxx[i] - I1/3.,
                    self.p_sigmayy[i] - I1/3.,
                    self.p_sigmazz[i] - I1/3.,
                    self.p_sigmaxy[i],
                    self.p_sigmaxz[i],
                    self.p_sigmayz[i],
                };
                const F J2 = 0.5*(s[0]*s[0] + s[1]*s[1] + s[2]*s[2] + 2.*(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]));
                // tensile correction 1
                if (J2 == 0. && I1 > k_c/alpha_phi) {
                    self.p_sigmaxx[i] = k_c/alpha_phi/3.;
                    self.p_sigmayy[i] = k_c/alpha_phi/3.;
                    self.p_sigmazz[i] = k_c/alpha_phi/3.;
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
                    
                    self.p_sigmaxx[i] -= dlambda*D0*((1.-self.m_v)*dgdsig[0] + self.m_v*dgdsig[1] + self.m_v*dgdsig[2]);
                    self.p_sigmayy[i] -= dlambda*D0*(self.m_v*dgdsig[0] + (1.-self.m_v)*dgdsig[1] + self.m_v*dgdsig[2]);
                    self.p_sigmazz[i] -= dlambda*D0*(self.m_v*dgdsig[0] + self.m_v*dgdsig[1] + (1.-self.m_v)*dgdsig[2]);
                    self.p_sigmaxy[i] -= dlambda*D0*dgdsig[3]*(1.-2.*self.m_v);
                    self.p_sigmaxz[i] -= dlambda*D0*dgdsig[4]*(1.-2.*self.m_v);
                    self.p_sigmayz[i] -= dlambda*D0*dgdsig[5]*(1.-2.*self.m_v);

                    I1 = self.p_sigmaxx[i] + self.p_sigmayy[i] + self.p_sigmazz[i];
                    if (I1 > k_c/alpha_phi) {
                        self.p_sigmaxx[i] = k_c/alpha_phi/3.;
                        self.p_sigmayy[i] = k_c/alpha_phi/3.;
                        self.p_sigmazz[i] = k_c/alpha_phi/3.;
                        self.p_sigmaxy[i] = 0.;
                        self.p_sigmaxz[i] = 0.;
                        self.p_sigmayz[i] = 0.;
                    }
                }

            }
        }
    }
}
#endif