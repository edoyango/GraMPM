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
        , dxdt {0.}
        , dydt {0.}
        , dzdt {0.}
        , mass {inmass}
        , rho {inrho}
        , sigmaxx {insigmaxx}
        , sigmayy {insigmayy}
        , sigmazz {insigmazz}
        , sigmaxy {insigmaxy}
        , sigmaxz {insigmaxz}
        , sigmayz {insigmayz}
        , strainratexx {0.}
        , strainrateyy {0.}
        , strainratezz {0.}
        , strainratexy {0.}
        , strainratexz {0.}
        , strainrateyz {0.}
        , spinratexy {0.}
        , spinratexz {0.}
        , spinrateyz {0.}
    {
    }

    template<typename F>
    particle<F>::particle(const F &inx, const F &iny, const F &inz, const F &invx, const F &invy, const F &invz, 
        const F &inmass, const F &inrho, const F &insigmaxx, const F &insigmayy, const F &insigmazz, const F &insigmaxy,
        const F &insigmaxz, const F &insigmayz, const F &inax, const F &inay, const F &inaz, const F &indxdt, 
        const F &indydt, const F &indzdt, const F &instrainratexx, const F &instrainrateyy, const F &instrainratezz, 
        const F &instrainratexy, const F &instrainratexz, const F &instrainrateyz, const F &inspinratexy,
        const F &inspinratexz, const F &inspinrateyz)
        : x {inx}
        , y {iny}
        , z {inz}
        , vx {invx}
        , vy {invy}
        , vz {invz}
        , ax {inax}
        , ay {inay}
        , az {inaz}
        , dxdt {indxdt}
        , dydt {indydt}
        , dzdt {indzdt}
        , mass {inmass}
        , rho {inrho}
        , sigmaxx {insigmaxx}
        , sigmayy {insigmayy}
        , sigmazz {insigmazz}
        , sigmaxy {insigmaxy}
        , sigmaxz {insigmaxz}
        , sigmayz {insigmayz}
        , strainratexx {instrainratexx}
        , strainrateyy {instrainrateyy}
        , strainratezz {instrainratezz}
        , strainratexy {instrainratexy}
        , strainratexz {instrainratexz}
        , strainrateyz {instrainrateyz}
        , spinratexy {inspinratexy}
        , spinratexz {inspinratexz}
        , spinrateyz {inspinrateyz}
    {
    }

    template<typename F> particle<F>::particle() {};
}

#endif