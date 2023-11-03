#ifndef GRAMPM_particle_ipp
#define GRAMPM_particle_ipp

namespace GraMPM {

    template<typename F>
    particle<F>::particle(const F &inx, const F &iny, const F &inz, const F &invx, const F &invy, const F &invz, 
        const F &inmass, const F &inrho, const F &insigmaxx, const F &insigmayy, const F &insigmazz, const F &insigmaxy,
        const F &insigmaxz, const F &insigmayz)
        : x {inx}
        , y {iny}
        , z {inz}
        , vx {invx}
        , vy {invy}
        , vz {invz}
        , mass {inmass}
        , rho {inrho}
        , sigmaxx {insigmaxx}
        , sigmayy {insigmayy}
        , sigmazz {insigmazz}
        , sigmaxy {insigmaxy}
        , sigmaxz {insigmaxz}
        , sigmayz {insigmayz}
    {
    }
}

#endif