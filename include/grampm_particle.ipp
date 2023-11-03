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
        , ax {0.}
        , ay {0.}
        , az {0.}
        , dvx {0.}
        , dvy {0.}
        , dvz {0.}
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

    template<typename F>
    particle<F>::particle(const F &inx, const F &iny, const F &inz, const F &invx, const F &invy, const F &invz, 
        const F &inmass, const F &inrho, const F &insigmaxx, const F &insigmayy, const F &insigmazz, const F &insigmaxy,
        const F &insigmaxz, const F &insigmayz, const F &inax, const F &inay, const F &inaz, const F &indvx, 
                const F & indvy, const F &indvz)
        : x {inx}
        , y {iny}
        , z {inz}
        , vx {invx}
        , vy {invy}
        , vz {invz}
        , ax {inax}
        , ay {inay}
        , az {inaz}
        , dvx {indvx}
        , dvy {indvy}
        , dvz {indvz}
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