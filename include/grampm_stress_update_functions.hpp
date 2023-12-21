#ifndef GRAMPM_STRESS_UPDATE
#define GRAMPM_STRESS_UPDATE

#include <grampm.hpp>

namespace GraMPM {
    namespace stress_update {

        template<typename F>
        void hookes_law(GraMPM::particle_system<F> &self, const F &dt) {

            const F E {self.E()}, v {self.v()};
            const int size {static_cast<int>(self.size())};
            std::vector<F> &sigmaxx {*(self.sigmaxx())}, &sigmayy {*(self.sigmayy())}, 
                &sigmazz {*(self.sigmazz())}, &sigmaxy {*(self.sigmaxy())}, 
                &sigmaxz {*(self.sigmaxz())}, &sigmayz {*(self.sigmayz())};
            std::vector<F> &strainratexx {*(self.strainratexx())}, &strainrateyy {*(self.strainrateyy())}, 
                &strainratezz {*(self.strainratezz())}, &strainratexy {*(self.strainratexy())}, 
                &strainratexz {*(self.strainratexz())}, &strainrateyz {*(self.strainrateyz())};
            std::vector<F> &spinratexy {*(self.spinratexy())}, &spinratexz {*(self.spinratexz())},
                &spinrateyz {*(self.spinrateyz())};

            const F D0 {E/((1.+v)*(1.-2.*v))};

            std::vector<F> dsigmaxx(size), dsigmayy(size), dsigmazz(size), dsigmaxy(size), dsigmaxz(size), 
                dsigmayz(size);

            // DE*dstrain
            #pragma omp for
            for (int i = 0; i < size; ++i) {
                dsigmaxx[i] = D0*((1.-v)*strainratexx[i] + v*strainrateyy[i] + v*strainratezz[i]);
                dsigmayy[i] = D0*(v*strainratexx[i] + (1.-v)*strainrateyy[i] + v*strainratezz[i]);
                dsigmazz[i] = D0*(v*strainratexx[i] + v*strainrateyy[i] + (1.-v)*strainratezz[i]);
                dsigmaxy[i] = D0*strainratexy[i]*(1.-2.*v);
                dsigmaxz[i] = D0*strainratexz[i]*(1.-2.*v);
                dsigmayz[i] = D0*strainrateyz[i]*(1.-2.*v);
            // }
            
            // jaumann stress rate
            // for (int i = 0; i < size; ++i) {
                dsigmaxx[i] -= 2.*(spinratexy[i]*sigmaxy[i] + spinratexz[i]*sigmaxz[i]);
                dsigmayy[i] -= 2.*(-spinratexy[i]*sigmaxy[i] + spinrateyz[i]*sigmayz[i]);
                dsigmazz[i] += 2.*(spinratexz[i]*sigmaxz[i] + spinrateyz[i]*sigmayz[i]);
                dsigmaxy[i] += sigmaxx[i]*spinratexy[i] - sigmaxz[i]*spinrateyz[i] -
                    spinratexy[i]*sigmayy[i] - spinratexz[i]*sigmayz[i];
                dsigmaxz[i] += sigmaxx[i]*spinratexz[i] + sigmaxy[i]*spinrateyz[i] -
                    spinratexy[i]*sigmayz[i] - spinratexz[i]*sigmazz[i];
                dsigmayz[i] += sigmaxy[i]*spinratexz[i] + sigmayy[i]*spinrateyz[i] +
                    spinratexy[i]*sigmaxz[i] - spinrateyz[i]*sigmazz[i];
            // }

            // update original stress states
            // for (int i = 0; i < size; ++i) {
                sigmaxx[i] += dt*dsigmaxx[i];
                sigmayy[i] += dt*dsigmayy[i];
                sigmazz[i] += dt*dsigmazz[i];
                sigmaxy[i] += dt*dsigmaxy[i];
                sigmaxz[i] += dt*dsigmaxz[i];
                sigmayz[i] += dt*dsigmayz[i];
            }

        }

        template<typename F>
        void drucker_prager_elastoplastic(GraMPM::particle_system<F> &self, const F &dt) {

            // performing initial elastic predictor step
            hookes_law(self, dt);

            // setting up for plastic corrector step (duplicated code, would like to remove this somehow)
            const F E {self.E()}, v {self.v()};
            const int size {static_cast<int>(self.size())};
            std::vector<F> &sigmaxx {*(self.sigmaxx())}, &sigmayy {*(self.sigmayy())}, 
                &sigmazz {*(self.sigmazz())}, &sigmaxy {*(self.sigmaxy())}, 
                &sigmaxz {*(self.sigmaxz())}, &sigmayz {*(self.sigmayz())};
            std::vector<F> &strainratexx {*(self.strainratexx())}, &strainrateyy {*(self.strainrateyy())}, 
                &strainratezz {*(self.strainratezz())}, &strainratexy {*(self.strainratexy())}, 
                &strainratexz {*(self.strainratexz())}, &strainrateyz {*(self.strainrateyz())};
            std::vector<F> &spinratexy {*(self.spinratexy())}, &spinratexz {*(self.spinratexz())},
                &spinrateyz {*(self.spinrateyz())};

            const F D0 {E/((1.+v)*(1.-2.*v))};

            F phi, psi, coh, alpha_phi, alpha_psi, k_c;
            self.DP_params(phi, psi, coh, alpha_phi, alpha_psi, k_c);

            // begin plastic correction
            #pragma omp for
            for (int i = 0; i < size; ++i) {
                // calculating invariants and deviatoric stress tensor
                F I1 = sigmaxx[i] + sigmayy[i] + sigmazz[i];
                const F s[6] {
                    sigmaxx[i] - I1/3.,
                    sigmayy[i] - I1/3.,
                    sigmazz[i] - I1/3.,
                    sigmaxy[i],
                    sigmaxz[i],
                    sigmayz[i],
                };
                const F J2 = 0.5*(s[0]*s[0] + s[1]*s[1] + s[2]*s[2] + 2.*(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]));
                // tensile correction 1
                if (J2 == 0. && I1 > k_c/alpha_phi) {
                    sigmaxx[i] = k_c/alpha_phi/3.;
                    sigmayy[i] = k_c/alpha_phi/3.;
                    sigmazz[i] = k_c/alpha_phi/3.;
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
                    
                    sigmaxx[i] -= dlambda*D0*((1.-v)*dgdsig[0] + v*dgdsig[1] + v*dgdsig[2]);
                    sigmayy[i] -= dlambda*D0*(v*dgdsig[0] + (1.-v)*dgdsig[1] + v*dgdsig[2]);
                    sigmazz[i] -= dlambda*D0*(v*dgdsig[0] + v*dgdsig[1] + (1.-v)*dgdsig[2]);
                    sigmaxy[i] -= dlambda*D0*dgdsig[3]*(1.-2.*v);
                    sigmaxz[i] -= dlambda*D0*dgdsig[4]*(1.-2.*v);
                    sigmayz[i] -= dlambda*D0*dgdsig[5]*(1.-2.*v);

                    I1 = sigmaxx[i] + sigmayy[i] + sigmazz[i];
                    if (I1 > k_c/alpha_phi) {
                        sigmaxx[i] = k_c/alpha_phi/3.;
                        sigmayy[i] = k_c/alpha_phi/3.;
                        sigmazz[i] = k_c/alpha_phi/3.;
                        sigmaxy[i] = 0.;
                        sigmaxz[i] = 0.;
                        sigmayz[i] = 0.;
                    }
                    I1 = sigmaxx[i] + sigmayy[i] + sigmazz[i];
                    const F stest[6] {
                        sigmaxx[i] - I1/3.,
                        sigmayy[i] - I1/3.,
                        sigmazz[i] - I1/3.,
                        sigmaxy[i],
                        sigmaxz[i],
                        sigmayz[i],
                    };
                }

            }
        }
    }
}
#endif