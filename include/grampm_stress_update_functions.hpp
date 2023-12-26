#ifndef GRAMPM_STRESS_UPDATE
#define GRAMPM_STRESS_UPDATE

#include <grampm.hpp>

namespace GraMPM {
    namespace stress_update {

        template<typename F>
        void hookes_law(GraMPM::MPM_system<F> &self, const F &dt) {

            const F E = self.get_stress_update_param("E"), v = self.get_stress_update_param("v");

            const F D0 {E/((1.+v)*(1.-2.*v))};

            // DE*dstrain
            #pragma omp for
            for (size_t i = 0; i < self.p_size(); ++i) {

                // update stress state with elastic increment
                double dsigmaxx = D0*((1.-v)*self.p_strainratexx(i) + v*self.p_strainrateyy(i) + v*self.p_strainratezz(i));
                double dsigmayy = D0*(v*self.p_strainratexx(i) + (1.-v)*self.p_strainrateyy(i) + v*self.p_strainratezz(i));
                double dsigmazz = D0*(v*self.p_strainratexx(i) + v*self.p_strainrateyy(i) + (1.-v)*self.p_strainratezz(i));
                double dsigmaxy = D0*self.p_strainratexy(i)*(1.-2.*v);
                double dsigmaxz = D0*self.p_strainratexz(i)*(1.-2.*v);
                double dsigmayz = D0*self.p_strainrateyz(i)*(1.-2.*v);

                // jaumann stress rate
                dsigmaxx -= 2.*(self.p_spinratexy(i)*self.p_sigmaxy(i) + self.p_spinratexz(i)*self.p_sigmaxz(i));
                dsigmayy -= 2.*(-self.p_spinratexy(i)*self.p_sigmaxy(i) + self.p_spinrateyz(i)*self.p_sigmayz(i));
                dsigmazz += 2.*(self.p_spinratexz(i)*self.p_sigmaxz(i) + self.p_spinrateyz(i)*self.p_sigmayz(i));
                dsigmaxy += self.p_sigmaxx(i)*self.p_spinratexy(i) - self.p_sigmaxz(i)*self.p_spinrateyz(i) -
                    self.p_spinratexy(i)*self.p_sigmayy(i) - self.p_spinratexz(i)*self.p_sigmayz(i);
                dsigmaxz += self.p_sigmaxx(i)*self.p_spinratexz(i) + self.p_sigmaxy(i)*self.p_spinrateyz(i) -
                    self.p_spinratexy(i)*self.p_sigmayz(i) - self.p_spinratexz(i)*self.p_sigmazz(i);
                dsigmayz += self.p_sigmaxy(i)*self.p_spinratexz(i) + self.p_sigmayy(i)*self.p_spinrateyz(i) +
                    self.p_spinratexy(i)*self.p_sigmaxz(i) - self.p_spinrateyz(i)*self.p_sigmazz(i);

                self.p_sigmaxx(i) += dt*dsigmaxx;
                self.p_sigmayy(i) += dt*dsigmayy;
                self.p_sigmazz(i) += dt*dsigmazz;
                self.p_sigmaxy(i) += dt*dsigmaxy;
                self.p_sigmaxz(i) += dt*dsigmaxz;
                self.p_sigmayz(i) += dt*dsigmayz;
            }

        }

        template<typename F>
        void drucker_prager_elastoplastic(GraMPM::MPM_system<F> &self, const F &dt) {

            // performing initial elastic predictor step
            hookes_law(self, dt);

            // setting up for plastic corrector step (duplicated code, would like to remove this somehow)
            const F E = self.get_stress_update_param("E"), v = self.get_stress_update_param("v");
            const F D0 {E/((1.+v)*(1.-2.*v))};

            const F phi = self.get_stress_update_param("phi"), psi = self.get_stress_update_param("psi"), 
                coh = self.get_stress_update_param("cohesion");
            const F alpha_phi = 2.*std::sin(phi)/(std::sqrt(3.)*(3.-std::sin(phi))), 
                alpha_psi = 2.*std::sin(psi)/(std::sqrt(3.)*(3.-std::sin(phi))), 
                k_c = 6.*coh*std::cos(phi)/(std::sqrt(3.)*(3.-std::sin(phi)));

            // begin plastic correction
            #pragma omp for
            for (size_t i = 0; i < self.p_size(); ++i) {
                // calculating invariants and deviatoric stress tensor
                F I1 = self.p_sigmaxx(i) + self.p_sigmayy(i) + self.p_sigmazz(i);
                const F s[6] {
                    self.p_sigmaxx(i) - I1/3.,
                    self.p_sigmayy(i) - I1/3.,
                    self.p_sigmazz(i) - I1/3.,
                    self.p_sigmaxy(i),
                    self.p_sigmaxz(i),
                    self.p_sigmayz(i),
                };
                const F J2 = 0.5*(s[0]*s[0] + s[1]*s[1] + s[2]*s[2] + 2.*(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]));
                // tensile correction 1
                if (J2 == 0. && I1 > k_c/alpha_phi) {
                    self.p_sigmaxx(i) = k_c/alpha_phi/3.;
                    self.p_sigmayy(i) = k_c/alpha_phi/3.;
                    self.p_sigmazz(i) = k_c/alpha_phi/3.;
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
                        dfdsig[0]*D0*((1.-v)*dgdsig[0] + v*dgdsig[1] + v*dgdsig[2]) +
                        dfdsig[1]*D0*(v*dgdsig[0] + (1.-v)*dgdsig[1] + v*dgdsig[2]) +
                        dfdsig[2]*D0*(v*dgdsig[0] + v*dgdsig[1] + (1.-v)*dgdsig[2]) +
                        2.*dfdsig[3]*D0*dgdsig[3]*(1.-2.*v) +
                        2.*dfdsig[4]*D0*dgdsig[4]*(1.-2.*v) +
                        2.*dfdsig[5]*D0*dgdsig[5]*(1.-2.*v)
                    )};
                    
                    self.p_sigmaxx(i) -= dlambda*D0*((1.-v)*dgdsig[0] + v*dgdsig[1] + v*dgdsig[2]);
                    self.p_sigmayy(i) -= dlambda*D0*(v*dgdsig[0] + (1.-v)*dgdsig[1] + v*dgdsig[2]);
                    self.p_sigmazz(i) -= dlambda*D0*(v*dgdsig[0] + v*dgdsig[1] + (1.-v)*dgdsig[2]);
                    self.p_sigmaxy(i) -= dlambda*D0*dgdsig[3]*(1.-2.*v);
                    self.p_sigmaxz(i) -= dlambda*D0*dgdsig[4]*(1.-2.*v);
                    self.p_sigmayz(i) -= dlambda*D0*dgdsig[5]*(1.-2.*v);

                    I1 = self.p_sigmaxx(i) + self.p_sigmayy(i) + self.p_sigmazz(i);
                    if (I1 > k_c/alpha_phi) {
                        self.p_sigmaxx(i) = k_c/alpha_phi/3.;
                        self.p_sigmayy(i) = k_c/alpha_phi/3.;
                        self.p_sigmazz(i) = k_c/alpha_phi/3.;
                        self.p_sigmaxy(i) = 0.;
                        self.p_sigmaxz(i) = 0.;
                        self.p_sigmayz(i) = 0.;
                    }
                }

            }
        }

        template<typename F>
        void eos(GraMPM::MPM_system<F> &self, const F &dt) {
            const F c = self.get_stress_update_param("c"), refrho = self.get_stress_update_param("reference-density"),
                mu = self.get_stress_update_param("kinematic viscosity");
            #pragma omp for
            for (size_t i = 0; i < self.p_size(); ++i) {
                F p = c*c*(self.p_rho(i) - refrho);
                const F divv = self.p_strainratexx(i) + self.p_strainrateyy(i) + self.p_strainratezz(i);
                if (divv < 0.) p += 4.*self.p_rho(i)*self.g_cell_size()*self.g_cell_size()*divv*divv + 
                    1.*self.p_rho(i)*self.g_cell_size()*c*std::abs(divv);
                self.p_sigmaxx(i) = -p+2.*mu;
                self.p_sigmayy(i) = -p+2.*mu;
                self.p_sigmazz(i) = -p+2.*mu;
                self.p_sigmaxy(i) = 2.*mu;
                self.p_sigmaxz(i) = 2.*mu;
                self.p_sigmayz(i) = 2.*mu;
            }
        }
    }
}
#endif