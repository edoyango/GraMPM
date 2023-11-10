#ifndef GRAMPM_STRESS_UPDATE
#define GRAMPM_STRESS_UPDATE

#include <grampm.hpp>

namespace GraMPM {
    namespace stress_update {
        template<typename F>
        void hookes_law(GraMPM::particle_system<F> &self, const F &dt) {

            const F E {self.E()}, v {self.v()};
            const int size {self.size()};
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
            for (int i = 0; i < size; ++i) {
                dsigmaxx[i] = D0*((1.-v)*strainratexx[i] + v*strainrateyy[i] + v*strainratezz[i]);
                dsigmayy[i] = D0*(v*strainratexx[i] + (1.-v)*strainrateyy[i] + v*strainratezz[i]);
                dsigmazz[i] = D0*(v*strainratexx[i] + v*strainrateyy[i] + (1.-v)*strainratezz[i]);
                dsigmaxy[i] = D0*strainratexy[i]*(1.-2.*v);
                dsigmaxz[i] = D0*strainratexz[i]*(1.-2.*v);
                dsigmayz[i] = D0*strainrateyz[i]*(1.-2.*v);
            }
            
            // jaumann stress rate
            for (int i = 0; i < size; ++i) {
                dsigmaxx[i] -= 2.*(spinratexy[i]*sigmaxy[i] + spinratexz[i]*sigmaxz[i]);
                dsigmayy[i] -= 2.*(-spinratexy[i]*sigmaxy[i] + spinrateyz[i]*sigmayz[i]);
                dsigmazz[i] += 2.*(spinratexz[i]*sigmaxz[i] + spinrateyz[i]*sigmayz[i]);
                dsigmaxy[i] += sigmaxx[i]*spinratexy[i] - sigmaxz[i]*spinrateyz[i] -
                    spinratexy[i]*sigmayy[i] - spinratexz[i]*sigmayz[i];
                dsigmaxz[i] += sigmaxx[i]*spinratexz[i] + sigmaxy[i]*spinrateyz[i] -
                    spinratexy[i]*sigmayz[i] - spinratexz[i]*sigmazz[i];
                dsigmayz[i] += sigmaxy[i]*spinratexz[i] + sigmayy[i]*spinrateyz[i] +
                    spinratexy[i]*sigmaxz[i] - spinrateyz[i]*sigmazz[i];
            }

            // update original stress states
            for (int i = 0; i < size; ++i) {
                sigmaxx[i] += dt*dsigmaxx[i];
                sigmayy[i] += dt*dsigmayy[i];
                sigmazz[i] += dt*dsigmazz[i];
                sigmaxy[i] += dt*dsigmaxy[i];
                sigmaxz[i] += dt*dsigmaxz[i];
                sigmayz[i] += dt*dsigmayz[i];
            }

        }
    }
}
#endif